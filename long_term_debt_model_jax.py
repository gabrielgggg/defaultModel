"""
Sovereign default model with long-term debt (taste shocks / discrete choice).

JAX port of long_term_debt_model.py (itself a port of _fortran/src/).
Hot path is pure jax.numpy + jax.jit so the same code can run on CPU or GPU.

Install (CPU)::
    pip install -U "jax[cpu]"

GPU: follow https://docs.jax.dev/en/latest/installation.html for CUDA wheels.
"""

from __future__ import annotations

import math
from functools import partial

import jax
import jax.numpy as jnp
import matplotlib.pyplot as plt
import numpy as np
from jax import lax, random
from jax.scipy.special import logsumexp
from jax.scipy.stats import norm as jax_norm
from matplotlib.backends.backend_pdf import PdfPages

# Prefer float64 for parity with Fortran / Numba port
jax.config.update("jax_enable_x64", True)

PLOT_PDF_PATH = "long_term_debt_model_jax_plots.pdf"

# ---------------------------------------------------------------------------
# Parameters (match defMod.f90)
# ---------------------------------------------------------------------------
CRRA = 2.0
LBD0 = -0.48
LBD1 = 0.525
RHO_Y = 0.95
SIGMA_Y = 0.005
RF = 1.04**0.25 - 1.0
BETA = 0.9775
GAMM = 1.0 / (2.0 * 4.0)

MACAULAY = 5.0 * 4.0
DELTA = (1.0 + RF) / MACAULAY - RF
KAPPA = DELTA + RF

RHO_D = 5.0e-4
RHO_B = 1.0e-5

MAX_ITER = 1000
SIM_SZ = 100_000
Y_SZ = 31
B_SZ = 600
TOL_ERR_V = 1.0e-6
TOL_ERR_Q = 1.0e-6
B_MIN = 0.0
B_MAX = 0.75
VERY_NEGATIVE = -1.0e6

# Set True for a fast end-to-end smoke test (coarse grids)
DEBUG_SMALL = False
if DEBUG_SMALL:
    Y_SZ = 7
    B_SZ = 50
    SIM_SZ = 5_000
    MAX_ITER = 200

RNG_SEED = 1989


# ---------------------------------------------------------------------------
# Utility and grids
# ---------------------------------------------------------------------------
def u_fun(cons: jnp.ndarray | float) -> jnp.ndarray | float:
    cons = jnp.asarray(cons, dtype=jnp.float64)
    return (cons ** (1.0 - CRRA) - 1.0) / (1.0 - CRRA)


def h_fun(y_val: jnp.ndarray | float) -> jnp.ndarray | float:
    y_val = jnp.asarray(y_val, dtype=jnp.float64)
    return y_val - jnp.maximum(0.0, LBD0 * y_val + LBD1 * y_val**2)


def discretize_ar1(
    mean: float,
    rho: float,
    sigma: float,
    n: int,
    n_sds: float,
) -> tuple[jnp.ndarray, jnp.ndarray, jnp.ndarray]:
    """Tauchen-style AR(1) discretization (mid-point bins), as in NL.f90."""
    if n == 1:
        grid = jnp.array([mean], dtype=jnp.float64)
        tran = jnp.array([[1.0]], dtype=jnp.float64)
        stationary = jnp.array([1.0], dtype=jnp.float64)
        return grid, tran, stationary

    half_width = n_sds * sigma / math.sqrt(1.0 - rho * rho)
    grid = jnp.linspace(mean - half_width, mean + half_width, n, dtype=jnp.float64)

    # Midpoints between grid points: length n-1
    mid = 0.5 * (grid[:-1] + grid[1:])
    cond_mean = (1.0 - rho) * mean + rho * grid[:, None]  # (n, 1)

    # P(j=0): CDF at first midpoint
    # P(j=n-1): 1 - CDF at last midpoint
    # P interior: CDF(hi) - CDF(lo)
    z_first = (mid[0] - cond_mean[:, 0]) / sigma
    z_last = (mid[-1] - cond_mean[:, 0]) / sigma
    p0 = jax_norm.cdf(z_first)
    p_last = 1.0 - jax_norm.cdf(z_last)

    # Interior columns j = 1 .. n-2
    if n > 2:
        # mid[j-1] is lo, mid[j] is hi for state j
        z_mids = (mid[None, :] - cond_mean) / sigma  # (n, n-1)
        cdf_mids = jax_norm.cdf(z_mids)
        p_int = cdf_mids[:, 1:] - cdf_mids[:, :-1]  # (n, n-2)
        tran = jnp.concatenate(
            [p0[:, None], p_int, p_last[:, None]], axis=1
        )
    else:
        tran = jnp.stack([p0, p_last], axis=1)

    # Stationary distribution by power iteration
    stationary = jnp.zeros(n, dtype=jnp.float64).at[max(0, n // 2)].set(1.0)

    def body(carry):
        st, _err = carry
        stt0 = st @ tran
        err = jnp.max(jnp.abs(stt0 - st))
        return stt0, err

    def cond(carry):
        _st, err = carry
        return err > 1.0e-12

    stationary, _ = lax.while_loop(cond, body, (stationary, jnp.array(1.0)))
    return grid, tran, stationary


def prepare_shocks_and_grids() -> tuple[jnp.ndarray, jnp.ndarray, jnp.ndarray, jnp.ndarray]:
    b_grid = jnp.linspace(B_MIN, B_MAX, B_SZ, dtype=jnp.float64)
    y_grid, y_pi, y_stat = discretize_ar1(0.0, RHO_Y, SIGMA_Y, Y_SZ, 3.0)
    # Level adjustment as in defMod.f90 prepareShocksAndGrids
    y_grid = jnp.exp(y_grid - 0.5 * SIGMA_Y**2 / (1.0 - RHO_Y**2))
    return y_grid, y_pi, y_stat, b_grid


# ---------------------------------------------------------------------------
# VFI kernels (vectorized + jitted)
# ---------------------------------------------------------------------------
@jax.jit
def update_default_value(
    uhy: jnp.ndarray,
    y_pi: jnp.ndarray,
    v0: jnp.ndarray,
    vd0: jnp.ndarray,
    beta: float,
    gamm: float,
) -> jnp.ndarray:
    cont = y_pi @ (gamm * v0[:, 0] + (1.0 - gamm) * vd0)
    return uhy + beta * cont


@partial(jax.jit, static_argnames=("crra",))
def update_values(
    y_grid: jnp.ndarray,
    b_grid: jnp.ndarray,
    y_pi: jnp.ndarray,
    v0: jnp.ndarray,
    vd1: jnp.ndarray,
    q0: jnp.ndarray,
    beta: float,
    delta: float,
    kappa: float,
    rho_b: float,
    rho_d: float,
    crra: float,
    very_negative: float,
) -> tuple[jnp.ndarray, jnp.ndarray, jnp.ndarray, jnp.ndarray]:
    """Return vr, d_pol, v1, b_pol."""
    # E[V0(y', b') | y] for each b'  → (Y, B)
    ev = y_pi @ v0

    y = y_grid[:, None, None]  # (Y, 1, 1)
    b = b_grid[None, :, None]  # (1, B, 1)
    bp = b_grid[None, None, :]  # (1, 1, B)
    q = q0[:, None, :]  # (Y, 1, B)

    cc = y - kappa * b + q * (bp - (1.0 - delta) * b)
    u_cc = (cc ** (1.0 - crra) - 1.0) / (1.0 - crra)
    w = jnp.where(cc > 0.0, u_cc + beta * ev[:, None, :], very_negative)

    # Taste shocks over b' (rho_b): log-sum-exp
    wbar = jnp.max(w, axis=-1, keepdims=True)
    log_ex = (w - wbar) / rho_b
    log_sum = logsumexp(log_ex, axis=-1, keepdims=True)
    b_pol = jnp.exp(log_ex - log_sum)
    vr = jnp.squeeze(wbar, axis=-1) + rho_b * jnp.squeeze(log_sum, axis=-1)

    # Default vs repay (rho_d)
    vd = vd1[:, None]  # (Y, 1)
    wbar2 = jnp.maximum(vd, vr)
    e_d = jnp.exp((vd - wbar2) / rho_d)
    e_r = jnp.exp((vr - wbar2) / rho_d)
    the_sum2 = e_d + e_r
    d_pol = e_d / the_sum2
    v1 = wbar2 + rho_d * jnp.log(the_sum2)
    return vr, d_pol, v1, b_pol


@jax.jit
def update_prices(
    y_pi: jnp.ndarray,
    d_pol: jnp.ndarray,
    b_pol: jnp.ndarray,
    q0: jnp.ndarray,
    kappa: float,
    delta: float,
    rf: float,
) -> jnp.ndarray:
    # qbar[yp, bp] = sum_bpp b_pol[yp, bp, bpp] * q0[yp, bpp]
    qbar = jnp.sum(b_pol * q0[:, None, :], axis=-1)
    payoff = (1.0 - d_pol) * (kappa + (1.0 - delta) * qbar)
    return (y_pi @ payoff) / (1.0 + rf)


@jax.jit
def expected_b_prime(b_grid: jnp.ndarray, b_pol: jnp.ndarray) -> jnp.ndarray:
    return b_pol @ b_grid


@partial(jax.jit, static_argnames=("crra",))
def vfi_step(
    uhy: jnp.ndarray,
    y_grid: jnp.ndarray,
    b_grid: jnp.ndarray,
    y_pi: jnp.ndarray,
    v0: jnp.ndarray,
    vd0: jnp.ndarray,
    q0: jnp.ndarray,
    beta: float,
    gamm: float,
    delta: float,
    kappa: float,
    rho_b: float,
    rho_d: float,
    crra: float,
    very_negative: float,
    rf: float,
) -> tuple[
    jnp.ndarray,
    jnp.ndarray,
    jnp.ndarray,
    jnp.ndarray,
    jnp.ndarray,
    jnp.ndarray,
    jnp.ndarray,
    jnp.ndarray,
]:
    vd1 = update_default_value(uhy, y_pi, v0, vd0, beta, gamm)
    vr, d_pol, v1, b_pol = update_values(
        y_grid,
        b_grid,
        y_pi,
        v0,
        vd1,
        q0,
        beta,
        delta,
        kappa,
        rho_b,
        rho_d,
        crra,
        very_negative,
    )
    q1 = update_prices(y_pi, d_pol, b_pol, q0, kappa, delta, rf)
    err_v = jnp.maximum(
        jnp.max(jnp.abs(v1 - v0)),
        jnp.max(jnp.abs(vd1 - vd0)),
    )
    err_q = jnp.max(jnp.abs(q1 - q0))
    return v1, vd1, q1, vr, d_pol, b_pol, err_v, err_q


# ---------------------------------------------------------------------------
# Simulation (lax.scan + jax.random)
# ---------------------------------------------------------------------------
def _draw_discrete(key: jax.Array, pmf: jnp.ndarray) -> tuple[jax.Array, jnp.ndarray]:
    """Inverse-CDF draw from a PMF; returns (new_key, index)."""
    key, sub = random.split(key)
    u = random.uniform(sub, shape=())
    cdf = jnp.cumsum(pmf)
    # searchsorted: first index where cdf >= u
    ix = jnp.searchsorted(cdf, u, side="left")
    ix = jnp.clip(ix, 0, pmf.shape[0] - 1)
    return key, ix


def simulate(
    y_pi: jnp.ndarray,
    y_grid: jnp.ndarray,
    b_grid: jnp.ndarray,
    d_pol: jnp.ndarray,
    b_pol: jnp.ndarray,
    q1: jnp.ndarray,
    gamm: float,
    kappa: float,
    delta: float,
    very_negative: float,
    lbd0: float,
    lbd1: float,
    sim_sz: int,
    seed: int,
) -> tuple[
    jnp.ndarray,
    jnp.ndarray,
    jnp.ndarray,
    jnp.ndarray,
    jnp.ndarray,
    jnp.ndarray,
    jnp.ndarray,
    jnp.ndarray,
]:
    y_sz = y_grid.shape[0]
    key = random.PRNGKey(seed)

    # CEILING(ySz/2) in Fortran (1-based) → 0-based mid index
    y0 = (y_sz + 1) // 2 - 1
    b0 = 0
    bpr0 = 0
    d0 = 0

    # t=0 row (written but mostly placeholders like Numba)
    init_carry = (
        jnp.int32(y0),
        jnp.int32(b0),
        jnp.int32(bpr0),
        jnp.int32(d0),
        key,
    )

    def step(carry, _t):
        y_prev, b_prev, b_pr_prev, d_prev, key = carry

        key, k_reentry = random.split(key)
        u_reentry = random.uniform(k_reentry)

        # If last period default: re-enter with prob gamm (zero debt), else stay out
        reenter = u_reentry <= gamm
        d_t = jnp.where(d_prev == 1, jnp.where(reenter, jnp.int32(0), jnp.int32(1)), jnp.int32(0))
        b_t = jnp.where(
            d_prev == 1,
            jnp.where(reenter, jnp.int32(0), b_prev),
            b_pr_prev,
        )

        key, y_t = _draw_discrete(key, y_pi[y_prev, :])

        # Fresh default decision if currently in market
        key, k_def = random.split(key)
        u_def = random.uniform(k_def)
        default_draw = u_def <= d_pol[y_t, b_t]
        d_t = jnp.where((d_t == 0) & default_draw, jnp.int32(1), d_t)

        # If in market: choose b', spreads, consumption
        key, b_pr_mkt = _draw_discrete(key, b_pol[y_t, b_t, :])
        q = q1[y_t, b_pr_mkt]
        sp_mkt = kappa * (1.0 / q - 1.0)
        gdp_mkt = y_grid[y_t]
        c_mkt = (
            gdp_mkt
            - kappa * b_grid[b_t]
            + q * (b_grid[b_pr_mkt] - (1.0 - delta) * b_grid[b_t])
        )
        tb_mkt = gdp_mkt - c_mkt

        # Default/exclusion: penalty on output
        yv = y_grid[y_t]
        pen = jnp.maximum(0.0, lbd0 * yv + lbd1 * yv * yv)
        gdp_def = yv - pen

        in_mkt = d_t == 0
        b_pr_t = jnp.where(in_mkt, b_pr_mkt, b_t)
        sp_t = jnp.where(in_mkt, sp_mkt, very_negative)
        gdp_t = jnp.where(in_mkt, gdp_mkt, gdp_def)
        c_t = jnp.where(in_mkt, c_mkt, gdp_def)
        tb_t = jnp.where(in_mkt, tb_mkt, 0.0)

        new_carry = (y_t, b_t, b_pr_t, d_t, key)
        out = (y_t, b_t, b_pr_t, d_t, sp_t, c_t, gdp_t, tb_t)
        return new_carry, out

    # scan produces t=1 .. sim_sz-1; prepend t=0
    (_final, outs) = lax.scan(step, init_carry, jnp.arange(1, sim_sz))
    y_rest, b_rest, bpr_rest, d_rest, sp_rest, c_rest, gdp_rest, tb_rest = outs

    y_sim = jnp.concatenate([jnp.array([y0], dtype=jnp.int32), y_rest])
    b_sim = jnp.concatenate([jnp.array([b0], dtype=jnp.int32), b_rest])
    b_pr_sim = jnp.concatenate([jnp.array([bpr0], dtype=jnp.int32), bpr_rest])
    d_sim = jnp.concatenate([jnp.array([d0], dtype=jnp.int32), d_rest])
    sp_sim = jnp.concatenate(
        [jnp.array([very_negative], dtype=jnp.float64), sp_rest]
    )
    c_sim = jnp.concatenate([jnp.array([0.0], dtype=jnp.float64), c_rest])
    gdp_sim = jnp.concatenate(
        [jnp.asarray([y_grid[y0]], dtype=jnp.float64), gdp_rest]
    )
    tb_sim = jnp.concatenate([jnp.array([0.0], dtype=jnp.float64), tb_rest])
    return y_sim, b_sim, b_pr_sim, d_sim, sp_sim, c_sim, gdp_sim, tb_sim


# JIT simulate with static sim_sz for fixed scan length
simulate = jax.jit(simulate, static_argnames=("sim_sz", "seed"))


# ---------------------------------------------------------------------------
# Solve
# ---------------------------------------------------------------------------
def solve_model():
    print("Start!")
    print(f"JAX devices: {jax.devices()}  backend={jax.default_backend()}")
    y_grid, y_pi, y_stat, b_grid = prepare_shocks_and_grids()
    hy = h_fun(y_grid)
    uhy = u_fun(hy)

    # Initial value: autarky-like consumption y - kappa b, floored
    cons0 = jnp.maximum(y_grid[:, None] - KAPPA * b_grid[None, :], 0.01)
    v0 = u_fun(cons0)
    vd0 = uhy.copy()
    q0 = jnp.ones((Y_SZ, B_SZ), dtype=jnp.float64)

    print("Compiling JAX kernels (first iteration may be slow)...")
    it = 1
    err_v = 1.0
    err_q = 1.0
    vr = d_pol = v1 = b_pol = vd1 = q1 = None

    while it <= MAX_ITER and (err_v > TOL_ERR_V or err_q > TOL_ERR_Q):
        v1, vd1, q1, vr, d_pol, b_pol, err_v_j, err_q_j = vfi_step(
            uhy,
            y_grid,
            b_grid,
            y_pi,
            v0,
            vd0,
            q0,
            BETA,
            GAMM,
            DELTA,
            KAPPA,
            RHO_B,
            RHO_D,
            CRRA,
            VERY_NEGATIVE,
            RF,
        )
        # Block until ready for host-side print / convergence check
        err_v = float(err_v_j)
        err_q = float(err_q_j)
        if it % 10 == 0:
            print(f"{it:5d}{err_v:20.5e}{err_q:20.5e}")
        it += 1
        v0 = v1
        vd0 = vd1
        q0 = q1

    print(
        f"Converged (or max iter) after {it - 1} iterations: "
        f"errV={err_v:.5e}, errQ={err_q:.5e}"
    )
    eb_pr = expected_b_prime(b_grid, b_pol)

    return {
        "y_grid": y_grid,
        "y_pi": y_pi,
        "y_stat": y_stat,
        "b_grid": b_grid,
        "V": v1,
        "Vr": vr,
        "Vd": vd1,
        "q": q1,
        "dPol": d_pol,
        "bPol": b_pol,
        "EbPr": eb_pr,
    }


# ---------------------------------------------------------------------------
# Moments (match moments.m / loadSimulation.m)
# ---------------------------------------------------------------------------
def compute_valid_mask(d_sim: np.ndarray, k: int = 20, n: int = 20) -> np.ndarray:
    """valid[ix] if no default in d_sim[ix-n : ix] inclusive (MATLAB ix-N:ix)."""
    d_sim = np.asarray(d_sim)
    sz = d_sim.shape[0]
    valid = np.zeros(sz, dtype=bool)
    for ix in range(k + n, sz):
        if np.sum(d_sim[ix - n : ix + 1]) == 0:
            valid[ix] = True
    return valid


def print_moments(
    b_grid,
    b_sim,
    b_pr_sim,
    d_sim,
    sp_sim,
    c_sim,
    gdp_sim,
    tb_sim,
) -> np.ndarray:
    b_grid = np.asarray(b_grid)
    b_sim = np.asarray(b_sim)
    d_sim = np.asarray(d_sim)
    sp_sim = np.asarray(sp_sim)
    c_sim = np.asarray(c_sim)
    gdp_sim = np.asarray(gdp_sim)
    tb_sim = np.asarray(tb_sim)

    sp_ann = (1.0 + sp_sim) ** 4 - 1.0
    valid = compute_valid_mask(d_sim)

    b_sim_val = b_grid[b_sim[valid]]
    gdp_val = gdp_sim[valid]
    loggdp = np.log(gdp_val)
    logc = np.log(c_sim[valid])
    tby = tb_sim[valid] / gdp_val
    sp_val = sp_ann[valid]

    print(f"Mean Debt to GDP   {100.0 * np.mean(b_sim_val / gdp_val / 4.0):10.2f}")
    print(f"Mean Spread        {100.0 * np.mean(sp_val):10.2f}")
    print(f"Std Spread         {100.0 * np.std(sp_val, ddof=1):10.2f}")
    print(f"Std log C          {100.0 * np.std(logc, ddof=1):10.2f}")
    print(f"Std log GDP        {100.0 * np.std(loggdp, ddof=1):10.2f}")
    print(f"Corr Sp, GDP       {100.0 * np.corrcoef(sp_val, loggdp)[0, 1]:10.2f}")
    print(f"Corr TB/GDP, GDP   {100.0 * np.corrcoef(tby, loggdp)[0, 1]:10.2f}")
    return valid


# ---------------------------------------------------------------------------
# Plots (match simplePlots.m)
# ---------------------------------------------------------------------------
def make_plots(
    eq: dict,
    sim: dict,
    valid: np.ndarray,
    pdf_path: str = PLOT_PDF_PATH,
) -> None:
    y_grid = np.asarray(eq["y_grid"])
    b_grid = np.asarray(eq["b_grid"])
    q = np.asarray(eq["q"])
    d_pol = np.asarray(eq["dPol"])
    b_pol = np.asarray(eq["bPol"])
    eb_pr = np.asarray(eq["EbPr"])

    y_sz = y_grid.shape[0]
    sp_sim = np.asarray(sim["sp_sim"])
    b_sim = np.asarray(sim["b_sim"])
    sp_ann = (1.0 + sp_sim) ** 4 - 1.0

    ry_ix = int(round(y_sz / 2.0)) - 1
    if DEBUG_SMALL:
        step = max(1, y_sz // 3)
        test_ys = [max(0, ry_ix - step), ry_ix, min(y_sz - 1, ry_ix + step)]
        b_fix_ix = min(149, B_SZ - 1)
    else:
        test_ys = [ry_ix - 10, ry_ix, ry_ix + 10]
        b_fix_ix = 149

    plt.rcParams.update({"font.size": 14})
    figures: list = []

    def y_label(yi: int) -> str:
        return f"y = {y_grid[yi]:.4f}"

    # 1. Spread histogram
    fig, ax = plt.subplots()
    ax.hist(sp_ann[valid], bins=50, density=True)
    ax.axvline(np.mean(sp_ann[valid]), linewidth=1)
    ax.set_xlim(0.0, 0.1)
    ax.set_xlabel("Spread")
    ax.set_ylabel("Frequency")
    ax.set_title("Spreads")
    fig.tight_layout()
    figures.append(fig)

    # 2. Debt histogram
    fig, ax = plt.subplots()
    debt = b_grid[b_sim[valid]]
    ax.hist(debt, bins=50, density=True)
    ax.axvline(np.mean(debt), linewidth=1)
    ax.set_xlabel("Debt")
    ax.set_ylabel("Frequency")
    ax.set_title("Debt")
    fig.tight_layout()
    figures.append(fig)

    # 3. q schedule
    fig, ax = plt.subplots()
    for yi in test_ys:
        ax.plot(b_grid, q[yi, :], label=y_label(yi))
    ax.set_title("q")
    ax.set_xlabel("Debt Next Period (B')")
    ax.legend()
    fig.tight_layout()
    figures.append(fig)

    # 4. Annualized spread curves
    fig, ax = plt.subplots()
    for yi in test_ys:
        q_safe = np.maximum(q[yi, :], 1.0e-14)
        sp = KAPPA * (1.0 / q_safe - 1.0)
        ax.plot(b_grid, (1.0 + sp) ** 4 - 1.0, linewidth=2, label=y_label(yi))
    ax.set_ylim(0.0, 0.1)
    ax.set_xlabel("Debt Next Period (B')")
    ax.set_ylabel("Spread")
    ax.legend()
    fig.tight_layout()
    figures.append(fig)

    # 5. b' choice probabilities at fixed current debt
    fig, ax = plt.subplots()
    for yi in test_ys:
        ax.plot(
            b_grid,
            b_pol[yi, b_fix_ix, :],
            "o-",
            linewidth=2,
            label=y_label(yi),
        )
    ax.axvline(b_grid[b_fix_ix], linewidth=1)
    if not DEBUG_SMALL:
        ax.set_xlim(0.17, 0.24)
    ax.set_xlabel("Debt")
    ax.set_ylabel("Choice Probability")
    ax.legend()
    fig.tight_layout()
    figures.append(fig)

    # 6. Expected next debt
    fig, ax = plt.subplots()
    for yi in test_ys:
        tmp = eb_pr[yi, :].copy()
        tmp[d_pol[yi, :] > 0.75] = np.nan
        ax.plot(b_grid, tmp, linewidth=2, label=y_label(yi))
    ax.plot(b_grid, b_grid, "--k", linewidth=1)
    ax.set_xlabel("Current Debt")
    ax.set_ylabel("Next Period Debt")
    if not DEBUG_SMALL:
        ax.set_xlim(0.0, 0.5)
    ax.legend()
    fig.tight_layout()
    figures.append(fig)

    with PdfPages(pdf_path) as pdf:
        for fig in figures:
            pdf.savefig(fig)
            plt.close(fig)
    print(f"Saved plots to {pdf_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    eq = solve_model()

    print("Simulate...")
    (
        y_sim,
        b_sim,
        b_pr_sim,
        d_sim,
        sp_sim,
        c_sim,
        gdp_sim,
        tb_sim,
    ) = simulate(
        eq["y_pi"],
        eq["y_grid"],
        eq["b_grid"],
        eq["dPol"],
        eq["bPol"],
        eq["q"],
        GAMM,
        KAPPA,
        DELTA,
        VERY_NEGATIVE,
        LBD0,
        LBD1,
        SIM_SZ,
        RNG_SEED,
    )
    # Materialize sim arrays
    y_sim, b_sim, b_pr_sim, d_sim, sp_sim, c_sim, gdp_sim, tb_sim = (
        np.asarray(y_sim),
        np.asarray(b_sim),
        np.asarray(b_pr_sim),
        np.asarray(d_sim),
        np.asarray(sp_sim),
        np.asarray(c_sim),
        np.asarray(gdp_sim),
        np.asarray(tb_sim),
    )

    # Match Fortran: write/use periods from t=300 (1-based) onward
    burn = 299
    y_sim = y_sim[burn:]
    b_sim = b_sim[burn:]
    b_pr_sim = b_pr_sim[burn:]
    d_sim = d_sim[burn:]
    sp_sim = sp_sim[burn:]
    c_sim = c_sim[burn:]
    gdp_sim = gdp_sim[burn:]
    tb_sim = tb_sim[burn:]

    print("\nMoments (valid sample, no recent default):")
    valid = print_moments(
        eq["b_grid"], b_sim, b_pr_sim, d_sim, sp_sim, c_sim, gdp_sim, tb_sim
    )

    sim = {
        "y_sim": y_sim,
        "b_sim": b_sim,
        "b_pr_sim": b_pr_sim,
        "d_sim": d_sim,
        "sp_sim": sp_sim,
        "c_sim": c_sim,
        "gdp_sim": gdp_sim,
        "tb_sim": tb_sim,
    }
    print("\nPlotting...")
    make_plots(eq, sim, valid)
    print("The end! ^_^")


if __name__ == "__main__":
    main()
