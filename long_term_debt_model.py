"""
Sovereign default model with long-term debt (taste shocks / discrete choice).

Python port of the Fortran implementation in _fortran/src/, with moments and
plots from _fortran/results/*.m. Uses NumPy, matplotlib, and Numba.
"""

from __future__ import annotations

import math

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from numba import njit, prange
from scipy.stats import norm

PLOT_PDF_PATH = "long_term_debt_model_plots.pdf"

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
def u_fun(cons: np.ndarray | float) -> np.ndarray | float:
    cons = np.asarray(cons, dtype=np.float64)
    if np.any(cons <= 0.0):
        raise ValueError("Negative consumption passed to u(.)")
    return (cons ** (1.0 - CRRA) - 1.0) / (1.0 - CRRA)


def h_fun(y_val: np.ndarray | float) -> np.ndarray | float:
    y_val = np.asarray(y_val, dtype=np.float64)
    return y_val - np.maximum(0.0, LBD0 * y_val + LBD1 * y_val**2)


def discretize_ar1(
    mean: float,
    rho: float,
    sigma: float,
    n: int,
    n_sds: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Tauchen-style AR(1) discretization (mid-point bins), as in NL.f90."""
    if n == 1:
        grid = np.array([mean], dtype=np.float64)
        tran = np.array([[1.0]], dtype=np.float64)
        stationary = np.array([1.0], dtype=np.float64)
        return grid, tran, stationary

    half_width = n_sds * sigma / math.sqrt(1.0 - rho * rho)
    grid = np.linspace(mean - half_width, mean + half_width, n, dtype=np.float64)
    tran = np.zeros((n, n), dtype=np.float64)

    for i in range(n):
        cond_mean = (1.0 - rho) * mean + rho * grid[i]
        for j in range(n):
            if j == 0:
                z = ((grid[1] + grid[0]) * 0.5 - cond_mean) / sigma
                tran[i, j] = norm.cdf(z)
            elif j == n - 1:
                z = ((grid[n - 1] + grid[n - 2]) * 0.5 - cond_mean) / sigma
                tran[i, j] = 1.0 - norm.cdf(z)
            else:
                z_lo = ((grid[j] + grid[j - 1]) * 0.5 - cond_mean) / sigma
                z_hi = ((grid[j + 1] + grid[j]) * 0.5 - cond_mean) / sigma
                tran[i, j] = norm.cdf(z_hi) - norm.cdf(z_lo)

    # Stationary distribution by power iteration
    stationary = np.zeros(n, dtype=np.float64)
    stationary[max(0, n // 2)] = 1.0
    err = 1.0
    while err > 1.0e-12:
        stt0 = stationary @ tran
        err = float(np.max(np.abs(stt0 - stationary)))
        stationary = stt0

    return grid, tran, stationary


def prepare_shocks_and_grids() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    b_grid = np.linspace(B_MIN, B_MAX, B_SZ, dtype=np.float64)
    y_grid, y_pi, y_stat = discretize_ar1(0.0, RHO_Y, SIGMA_Y, Y_SZ, 3.0)
    # Level adjustment as in defMod.f90 prepareShocksAndGrids
    y_grid = np.exp(y_grid - 0.5 * SIGMA_Y**2 / (1.0 - RHO_Y**2))
    return y_grid, y_pi, y_stat, b_grid


# ---------------------------------------------------------------------------
# Numba-accelerated VFI kernels
# ---------------------------------------------------------------------------
@njit(cache=True)
def _u_scalar(cons: float, crra: float) -> float:
    return (cons ** (1.0 - crra) - 1.0) / (1.0 - crra)


@njit(parallel=True, cache=True)
def update_default_value(
    uhy: np.ndarray,
    y_pi: np.ndarray,
    v0: np.ndarray,
    vd0: np.ndarray,
    beta: float,
    gamm: float,
    vd1: np.ndarray,
) -> None:
    y_sz = uhy.shape[0]
    for y_ix in prange(y_sz):
        cont = 0.0
        for yp in range(y_sz):
            cont += y_pi[y_ix, yp] * (gamm * v0[yp, 0] + (1.0 - gamm) * vd0[yp])
        vd1[y_ix] = uhy[y_ix] + beta * cont


@njit(parallel=True, cache=True)
def update_values(
    y_grid: np.ndarray,
    b_grid: np.ndarray,
    y_pi: np.ndarray,
    v0: np.ndarray,
    vd1: np.ndarray,
    q0: np.ndarray,
    beta: float,
    delta: float,
    kappa: float,
    rho_b: float,
    rho_d: float,
    crra: float,
    very_negative: float,
    vr: np.ndarray,
    d_pol: np.ndarray,
    v1: np.ndarray,
    b_pol: np.ndarray,
) -> None:
    y_sz = y_grid.shape[0]
    b_sz = b_grid.shape[0]

    for y_ix in prange(y_sz):
        # E[V0(y', b') | y] for each b'
        ev = np.empty(b_sz, dtype=np.float64)
        for bp in range(b_sz):
            s = 0.0
            for yp in range(y_sz):
                s += y_pi[y_ix, yp] * v0[yp, bp]
            ev[bp] = s

        y = y_grid[y_ix]
        for b_ix in range(b_sz):
            b = b_grid[b_ix]
            w_vec = np.empty(b_sz, dtype=np.float64)
            wbar = very_negative
            for bp in range(b_sz):
                cc = y - kappa * b + q0[y_ix, bp] * (b_grid[bp] - (1.0 - delta) * b)
                if cc <= 0.0:
                    w_vec[bp] = very_negative
                else:
                    w_vec[bp] = _u_scalar(cc, crra) + beta * ev[bp]
                if w_vec[bp] > wbar:
                    wbar = w_vec[bp]

            the_sum = 0.0
            for bp in range(b_sz):
                ex = math.exp((w_vec[bp] - wbar) / rho_b)
                b_pol[y_ix, b_ix, bp] = ex
                the_sum += ex
            inv_sum = 1.0 / the_sum
            for bp in range(b_sz):
                b_pol[y_ix, b_ix, bp] *= inv_sum
            vr[y_ix, b_ix] = wbar + rho_b * math.log(the_sum)

            # Default vs repay
            vd = vd1[y_ix]
            vrep = vr[y_ix, b_ix]
            wbar2 = vd if vd > vrep else vrep
            e_d = math.exp((vd - wbar2) / rho_d)
            e_r = math.exp((vrep - wbar2) / rho_d)
            the_sum2 = e_d + e_r
            d_pol[y_ix, b_ix] = e_d / the_sum2
            v1[y_ix, b_ix] = wbar2 + rho_d * math.log(the_sum2)


@njit(parallel=True, cache=True)
def update_prices(
    y_pi: np.ndarray,
    d_pol: np.ndarray,
    b_pol: np.ndarray,
    q0: np.ndarray,
    kappa: float,
    delta: float,
    rf: float,
    q1: np.ndarray,
) -> None:
    y_sz = y_pi.shape[0]
    b_sz = d_pol.shape[1]
    for y_ix in prange(y_sz):
        for bp in range(b_sz):
            s = 0.0
            for yp in range(y_sz):
                # expected continuation price: bPol(y', B') · q0(y', :)
                qbar = 0.0
                for bpp in range(b_sz):
                    qbar += b_pol[yp, bp, bpp] * q0[yp, bpp]
                s += y_pi[y_ix, yp] * (1.0 - d_pol[yp, bp]) * (
                    kappa + (1.0 - delta) * qbar
                )
            q1[y_ix, bp] = s / (1.0 + rf)


@njit(parallel=True, cache=True)
def expected_b_prime(b_grid: np.ndarray, b_pol: np.ndarray, eb_pr: np.ndarray) -> None:
    y_sz = b_pol.shape[0]
    b_sz = b_pol.shape[1]
    for y_ix in prange(y_sz):
        for b_ix in range(b_sz):
            s = 0.0
            for bp in range(b_sz):
                s += b_pol[y_ix, b_ix, bp] * b_grid[bp]
            eb_pr[y_ix, b_ix] = s


# ---------------------------------------------------------------------------
# Simulation (Numba)
# ---------------------------------------------------------------------------
@njit(cache=True)
def _sim_discrete(pmf: np.ndarray) -> int:
    draw = np.random.random()
    the_sum = 0.0
    n = pmf.shape[0]
    for i in range(n):
        the_sum += pmf[i]
        if draw <= the_sum:
            return i
    return n - 1


@njit(cache=True)
def simulate(
    y_pi: np.ndarray,
    y_grid: np.ndarray,
    b_grid: np.ndarray,
    d_pol: np.ndarray,
    b_pol: np.ndarray,
    q1: np.ndarray,
    gamm: float,
    kappa: float,
    delta: float,
    very_negative: float,
    lbd0: float,
    lbd1: float,
    sim_sz: int,
    seed: int,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    np.random.seed(seed)
    y_sz = y_grid.shape[0]
    b_sz = b_grid.shape[0]

    y_sim = np.empty(sim_sz, dtype=np.int64)
    b_sim = np.empty(sim_sz, dtype=np.int64)
    b_pr_sim = np.empty(sim_sz, dtype=np.int64)
    d_sim = np.empty(sim_sz, dtype=np.int64)
    sp_sim = np.empty(sim_sz, dtype=np.float64)
    c_sim = np.empty(sim_sz, dtype=np.float64)
    gdp_sim = np.empty(sim_sz, dtype=np.float64)
    tb_sim = np.empty(sim_sz, dtype=np.float64)

    # CEILING(ySz/2) in Fortran (1-based) → 0-based mid index
    y_sim[0] = (y_sz + 1) // 2 - 1
    b_sim[0] = 0
    b_pr_sim[0] = 0
    d_sim[0] = 0
    sp_sim[0] = very_negative
    c_sim[0] = 0.0
    gdp_sim[0] = y_grid[y_sim[0]]
    tb_sim[0] = 0.0

    for t in range(1, sim_sz):
        if d_sim[t - 1] == 1:
            if np.random.random() <= gamm:
                d_sim[t] = 0
                b_sim[t] = 0
            else:
                d_sim[t] = 1
                b_sim[t] = b_sim[t - 1]
        else:
            b_sim[t] = b_pr_sim[t - 1]
            d_sim[t] = 0

        y_sim[t] = _sim_discrete(y_pi[y_sim[t - 1], :])

        if d_sim[t] == 0:
            if np.random.random() <= d_pol[y_sim[t], b_sim[t]]:
                d_sim[t] = 1

        if d_sim[t] == 0:
            b_pr_sim[t] = _sim_discrete(b_pol[y_sim[t], b_sim[t], :])
            q = q1[y_sim[t], b_pr_sim[t]]
            sp_sim[t] = kappa * (1.0 / q - 1.0)
            gdp_sim[t] = y_grid[y_sim[t]]
            c_sim[t] = (
                gdp_sim[t]
                - kappa * b_grid[b_sim[t]]
                + q * (b_grid[b_pr_sim[t]] - (1.0 - delta) * b_grid[b_sim[t]])
            )
            tb_sim[t] = gdp_sim[t] - c_sim[t]
        else:
            b_pr_sim[t] = b_sim[t]
            sp_sim[t] = very_negative
            yv = y_grid[y_sim[t]]
            pen = lbd0 * yv + lbd1 * yv * yv
            if pen < 0.0:
                pen = 0.0
            gdp_sim[t] = yv - pen
            c_sim[t] = gdp_sim[t]
            tb_sim[t] = 0.0

    return y_sim, b_sim, b_pr_sim, d_sim, sp_sim, c_sim, gdp_sim, tb_sim


# ---------------------------------------------------------------------------
# Solve
# ---------------------------------------------------------------------------
def solve_model():
    print("Start!")
    y_grid, y_pi, y_stat, b_grid = prepare_shocks_and_grids()
    hy = h_fun(y_grid)
    uhy = np.asarray(u_fun(hy), dtype=np.float64)

    v0 = np.empty((Y_SZ, B_SZ), dtype=np.float64)
    v1 = np.empty_like(v0)
    vd0 = uhy.copy()
    vd1 = np.empty(Y_SZ, dtype=np.float64)
    vr = np.empty((Y_SZ, B_SZ), dtype=np.float64)
    d_pol = np.empty((Y_SZ, B_SZ), dtype=np.float64)
    q0 = np.ones((Y_SZ, B_SZ), dtype=np.float64)
    q1 = np.empty_like(q0)
    b_pol = np.empty((Y_SZ, B_SZ, B_SZ), dtype=np.float64)
    eb_pr = np.empty((Y_SZ, B_SZ), dtype=np.float64)

    for y_ix in range(Y_SZ):
        for b_ix in range(B_SZ):
            cons = max(y_grid[y_ix] - KAPPA * b_grid[b_ix], 0.01)
            v0[y_ix, b_ix] = float(u_fun(cons))

    # Warm-up JIT (first call compiles)
    print("Compiling Numba kernels (first iteration may be slow)...")

    it = 1
    err_v = 1.0
    err_q = 1.0
    while it <= MAX_ITER and (err_v > TOL_ERR_V or err_q > TOL_ERR_Q):
        update_default_value(uhy, y_pi, v0, vd0, BETA, GAMM, vd1)
        update_values(
            y_grid,
            b_grid,
            y_pi,
            v0,
            vd1,
            q0,
            BETA,
            DELTA,
            KAPPA,
            RHO_B,
            RHO_D,
            CRRA,
            VERY_NEGATIVE,
            vr,
            d_pol,
            v1,
            b_pol,
        )
        update_prices(y_pi, d_pol, b_pol, q0, KAPPA, DELTA, RF, q1)

        err_v = max(float(np.max(np.abs(v1 - v0))), float(np.max(np.abs(vd1 - vd0))))
        err_q = float(np.max(np.abs(q1 - q0)))
        if it % 10 == 0:
            print(f"{it:5d}{err_v:20.5e}{err_q:20.5e}")
        it += 1
        v0[:] = v1
        vd0[:] = vd1
        q0[:] = q1

    print(f"Converged (or max iter) after {it - 1} iterations: errV={err_v:.5e}, errQ={err_q:.5e}")
    expected_b_prime(b_grid, b_pol, eb_pr)

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
    sz = d_sim.shape[0]
    valid = np.zeros(sz, dtype=bool)
    for ix in range(k + n, sz):  # MATLAB: for ix = K+N+1:sz with 1-based → 0-based k+n
        # MATLAB sum(dSimIx(ix-N:ix)) == 0 → indices ix-N .. ix inclusive
        # 0-based: from ix-n to ix inclusive
        if np.sum(d_sim[ix - n : ix + 1]) == 0:
            valid[ix] = True
    return valid


def print_moments(
    b_grid: np.ndarray,
    b_sim: np.ndarray,
    b_pr_sim: np.ndarray,
    d_sim: np.ndarray,
    sp_sim: np.ndarray,
    c_sim: np.ndarray,
    gdp_sim: np.ndarray,
    tb_sim: np.ndarray,
) -> np.ndarray:
    # Drop burn-in as Fortran sim.tab (t from 300), then annualize spreads
    # Fortran writes from t=300 (1-based) → indices 299..end; we already trim in main
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
    y_sz = eq["y_grid"].shape[0]
    b_grid = eq["b_grid"]
    q = eq["q"]
    d_pol = eq["dPol"]
    b_pol = eq["bPol"]
    eb_pr = eq["EbPr"]

    sp_sim = sim["sp_sim"]
    b_sim = sim["b_sim"]
    # annualized for histogram (already done in moments path; recompute safely)
    sp_ann = (1.0 + sp_sim) ** 4 - 1.0

    ry_ix = int(round(y_sz / 2.0)) - 1  # MATLAB round(ySz/2) then 0-based
    if DEBUG_SMALL:
        # keep indices in range for small grids
        step = max(1, y_sz // 3)
        test_ys = [max(0, ry_ix - step), ry_ix, min(y_sz - 1, ry_ix + step)]
        b_fix_ix = min(149, B_SZ - 1)
    else:
        test_ys = [ry_ix - 10, ry_ix, ry_ix + 10]
        b_fix_ix = 149  # Fortran/MATLAB index 150

    y_grid = eq["y_grid"]
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
