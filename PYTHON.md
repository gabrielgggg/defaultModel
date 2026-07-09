# Python and JAX ports of the long-term debt model

This document describes the economic model and the work done to port the original **Fortran** implementation to **Python (NumPy/Numba)** and then to **JAX**. The ports are intended to reproduce the same quantitative object: a sovereign default model with long-term debt, solved with taste shocks / discrete-choice methods.

---

## 1. Which model is this?

### 1.1 Summary

This is a **sovereign default model with long-term debt**, solved by value-function iteration (VFI) with **extreme-value taste shocks** on both the default decision and the choice of next-period debt \(B'\). The continuous (or fine-grid discrete) policy is replaced by closed-form **logit / soft-max choice probabilities**, which is the discrete-choice approach common in quantitative macro (cf. taste-shock / “trembling hand” regularizations of discrete choice).

The Fortran reference tree lives in ``. Supporting write-ups are:

| File | Role |
|------|------|
| `README.md` | Model equations (primary narrative) |
| `longTermDebtModel.tex` | Same equations as LaTeX (Gabriel Mihalache, September 2022) |
| `src/` | Fortran implementation (OpenMP) |
| `results/*.m` | MATLAB post-processing (moments, plots) |

### 1.2 Economic environment

- **Agent:** Small open economy / sovereign that issues non-contingent long-term bonds to risk-neutral foreign creditors.
- **State:** Endowment (income) \(y\) and outstanding debt stock \(B\).
- **Income:** Markov process for \(\log y\), discretized by a Tauchen-style method (mid-point bins). After discretization, levels are adjusted so that the AR(1) on logs is consistent with the usual variance correction used in the Fortran code.
- **Long-term debt:** One-period price \(q(y, B')\) of a bond that pays a constant coupon stream and amortizes. Parameters \(\delta\) (decay / principal amortization) and \(\kappa = \delta + r_f\) are set so that, in a risk-free world, the bond has a given **Macaulay duration** (here 5 years at quarterly frequency → 20 quarters).
- **Default:** Binary default decision \(d \in \{0,1\}\). Default triggers a temporary **output penalty** and exclusion from markets until re-entry.
- **Exclusion / re-entry:** While excluded, the sovereign receives \(h(y)\) (penalized income). Each period it re-enters with probability \(\gamma\) and zero debt; otherwise it stays excluded.
- **Taste shocks:** i.i.d. extreme-value shocks on (i) default vs repay and (ii) each discrete \(B'\) alternative, with scales \(\rho_D\) and \(\rho_B\). This yields smooth choice probabilities and a closed-form “log-sum-exp” value of the choice set.

### 1.3 Recursive formulation

**Ex ante value** (default vs repay, soft-max with scale \(\rho_D\)):

\[
V(y,B)=\rho_D\log\left\{\exp\!\left[\frac{V^d(y)}{\rho_D}\right]+\exp\!\left[\frac{V^r(y,B)}{\rho_D}\right]\right\}
\]

**Default probability:**

\[
\Pr(d=1\mid y,B)=\frac{\exp\!\left[V^d(y)/\rho_D\right]}{\exp\!\left[V^d(y)/\rho_D\right]+\exp\!\left[V^r(y,B)/\rho_D\right]}
\]

**Value in default / exclusion:**

\[
V^d(y)=u\bigl[h(y)\bigr]+\beta\,\mathbb{E}_{y'|y}\Bigl\{\gamma\,V(y',0)+(1-\gamma)\,V^d(y')\Bigr\}
\]

**Continuation payoff of choosing next debt \(B'\) when repaying:**

\[
W(y,B,B')=u\Bigl[y-\kappa B+q(y,B')\bigl(B'-(1-\delta)B\bigr)\Bigr]+\beta\,\mathbb{E}_{y'|y}\,V(y',B')
\]

**Repayment value** (soft-max over \(B'\) with scale \(\rho_B\)):

\[
V^r(y,B)=\rho_B\log\sum_{B'}\exp\!\left[\frac{W(y,B,B')}{\rho_B}\right]
\]

**Debt choice probabilities:**

\[
\Pr(B'=x\mid y,B)=\frac{\exp\!\left[W(y,B,x)/\rho_B\right]}{\sum_i\exp\!\left[W(y,B,i)/\rho_B\right]}
\]

**Bond price** (risk-neutral creditors, recovery via non-default and remaining coupon stream):

\[
q(y,B')=\frac{1}{1+r_f}\,\mathbb{E}_{y'|y}\,\Pr(d=0\mid y',B')\Bigl[\kappa+(1-\delta)\sum_{B''}\Pr(B''\mid y',B')\,q(y',B'')\Bigr]
\]

### 1.4 Functional forms

| Object | Form |
|--------|------|
| Utility | \(u(c)=(c^{1-\sigma}-1)/(1-\sigma)\) (CRRA) |
| Default output | \(h(y)=y-\max\{0,\lambda_0 y+\lambda_1 y^2\}\) |
| Income process (logs) | \(\log y' = -(1-\rho)\frac{\sigma_y^2}{2(1-\rho^2)}+\rho\log y+\sigma_y\varepsilon\), \(\varepsilon\sim\mathcal{N}(0,1)\) |

(The level grid in code is obtained by discretizing an AR(1) on logs then applying the exponential / variance correction used in `prepareShocksAndGrids`.)

### 1.5 Calibration (code parameters)

These match `src/defMod.f90` and both Python ports:

| Symbol | Code name | Value | Notes |
|--------|-----------|------:|-------|
| \(\sigma\) | `CRRA` | 2.0 | Risk aversion |
| \(\lambda_0,\lambda_1\) | `LBD0`, `LBD1` | −0.48, 0.525 | Output cost of default |
| \(\rho_y,\sigma_y\) | `RHO_Y`, `SIGMA_Y` | 0.95, 0.005 | Income AR(1) (quarterly) |
| \(r_f\) | `RF` | \(1.04^{1/4}-1\) | Risk-free rate |
| \(\beta\) | `BETA` | 0.9775 | Discount factor |
| \(\gamma\) | `GAMM` | \(1/(2\cdot 4)\) | Re-entry hazard (mean exclusion ≈ 2 years) |
| Macaulay duration | `MACAULAY` | \(5\cdot 4\) quarters | Sets \(\delta,\kappa\) |
| \(\delta\) | `DELTA` | \((1+r_f)/\mathrm{Mac}-r_f\) | Bond amortization |
| \(\kappa\) | `KAPPA` | \(\delta+r_f\) | Coupon / payment rate |
| \(\rho_D,\rho_B\) | `RHO_D`, `RHO_B` | \(5\cdot 10^{-4}\), \(10^{-5}\) | Taste-shock scales |
| Grid \(y\) | `Y_SZ` | 31 | Income points |
| Grid \(B\) | `B_SZ` | 600 | Debt points on \([0, 0.75]\) |
| Simulation length | `SIM_SZ` | 100 000 | Periods |
| VFI tolerances | `TOL_ERR_V`, `TOL_ERR_Q` | \(10^{-6}\) | Max abs change in \(V\) and \(q\) |
| Max iterations | `MAX_ITER` | 1000 | Safety cap |
| RNG seed | `RNG_SEED` | 1989 | Simulation (stream differs by backend) |

**Frequency:** quarterly. Spreads in moments/plots are often **annualized** as \((1+\mathrm{sp})^4-1\). Debt/GDP in moments is scaled by 4 (quarterly debt stock vs annualized GDP convention matching the MATLAB scripts).

### 1.6 Solution algorithm (all three codes)

1. Build grids: \(B\) linspace; Tauchen-style AR(1) for income; apply level adjustment.
2. Initialize \(V_0(y,B)\approx u(\max\{y-\kappa B,0.01\})\), \(V^d_0=u(h(y))\), \(q_0\equiv 1\).
3. Iterate until \(\|V_1-V_0\|_\infty\) and \(\|q_1-q_0\|_\infty\) are below tolerance (or max iter):
   - Update \(V^d\) from \(V_0,V^d_0\).
   - Update \(W\), \(V^r\), \(b\mathrm{Pol}\), \(d\mathrm{Pol}\), \(V_1\) from \(V_0,q_0,V^d_1\).
   - Update \(q_1\) from \(d\mathrm{Pol}\), \(b\mathrm{Pol}\), \(q_0\).
4. Compute expected next debt \(\mathbb{E}[B'|y,B]=\sum_{B'} b\mathrm{Pol}\cdot B'\).
5. Simulate a long path with the equilibrium policies and prices.
6. Drop burn-in (Fortran writes from period \(t=300\) onward → index 299 in 0-based Python).
7. Compute **moments** on a “valid” sample: no default in a rolling window (MATLAB `valid` / `moments.m`).
8. Make **diagnostic plots** (histograms, \(q\), spreads, choice probs, policy).

---

## 2. Repository layout

```
default_grok/
├──                       # Original reference implementation
│   ├── README.md                  # Model equations
│   ├── longTermDebtModel.tex      # Same equations (PDF companion)
│   ├── Makefile
│   ├── src/
│   │   ├── defaultModel.f90       # Main VFI + driver
│   │   ├── defMod.f90             # Parameters, grids, simulate, I/O
│   │   ├── NL.f90                 # Numerics (e.g. AR(1) discretization)
│   │   └── sim.f90                # Random draws helpers
│   └── results/
│       ├── moments.m              # Business-cycle / debt moments
│       ├── simplePlots.m          # Diagnostic figures
│       └── …                      # loaders, binary I/O helpers
├── long_term_debt_model.py        # NumPy + Numba port
├── long_term_debt_model_jax.py    # JAX port (CPU/GPU-ready)
├── long_term_debt_model_plots.pdf     # Plots from Numba run (if executed)
├── long_term_debt_model_jax_plots.pdf # Plots from JAX run
└── PYTHON.md                      # This document
```

---

## 3. Port 1: Fortran → Python (`long_term_debt_model.py`)

### 3.1 Goals

- Reproduce the Fortran **economic logic**, **parameters**, **grids**, **VFI updates**, **simulation timing**, **moments**, and **plots** in a single Python module.
- Make the hot loops fast enough for the full grid (`Y=31`, `B=600`) on a laptop/desktop CPU without requiring a Fortran toolchain.
- Keep the code readable as a research reference (one file, clear sections).

### 3.2 Stack

| Layer | Choice | Role |
|-------|--------|------|
| Arrays | NumPy | Grids, host data, moments |
| JIT | Numba (`@njit`, `prange`) | VFI kernels + simulation loop |
| Distributions | SciPy `norm.cdf` | Tauchen transition matrix (setup only) |
| Plots | Matplotlib + `PdfPages` | Multi-page PDF matching MATLAB figures |

### 3.3 What was ported, piece by piece

| Fortran / MATLAB | Python | Notes |
|------------------|--------|--------|
| `defMod.f90` parameters | Module-level constants | Identical numerical values |
| `prepareShocksAndGrids` | `prepare_shocks_and_grids` | Linspace \(B\); AR(1); `exp` level fix |
| `discretizeAR1` (`NL.f90`) | `discretize_ar1` | Mid-point bin Tauchen; stationary power iteration |
| `uFun`, `hFun` | `u_fun`, `h_fun` | CRRA and quadratic default cost |
| VFI: \(V^d\) update | `update_default_value` | Numba parallel over \(y\) |
| VFI: \(W,V^r,bPol,dPol,V\) | `update_values` | Numba parallel over \(y\); inner loops over \(B,B'\) |
| VFI: \(q\) update | `update_prices` | Nested expectations over \(y',B''\) |
| `EbPr` | `expected_b_prime` | \(\sum bPol\cdot bGrid\) |
| `simulate` in `defMod` | `simulate` | Same branching: re-entry, default draw, \(B'\) draw |
| Burn-in \(t\ge 300\) | `burn = 299` in `main` | 0-based indexing |
| `moments.m` | `print_moments` + `compute_valid_mask` | Valid sample: no default in window |
| `simplePlots.m` | `make_plots` | Six figures → one PDF |

### 3.4 Implementation details

**In-place Numba kernels.**  
Fortran updates arrays in place with OpenMP. The Numba kernels take preallocated output arrays (`vd1`, `v1`, `b_pol`, …) and write into them. This avoids large temporary allocations in tight loops.

**Log-sum-exp stability.**  
For both \(\rho_B\) and \(\rho_D\) soft-maxes, the code subtracts `Wbar = max(·)` before `exp`, matching the Fortran pattern (`Wbar + rho * log(sum(exp((W-Wbar)/rho)))`).

**Invalid consumption.**  
If \(c\le 0\) for a candidate \(B'\), the continuation payoff is set to `VERY_NEGATIVE` (\(-10^6\)) rather than calling \(u(c)\), as in Fortran’s `WHERE (cc <= 0)`.

**Simulation discrete draws.**  
`_sim_discrete(pmf)` uses inverse CDF with `np.random.random()` inside Numba (seeded with `np.random.seed(RNG_SEED)`).

**DEBUG_SMALL.**  
Optional coarse grids (`Y=7`, `B=50`, shorter sim) for smoke tests without changing the full calibration defaults.

**Output.**  
Plots go to `long_term_debt_model_plots.pdf`. Moments print to stdout in the same labels/units as `moments.m`.

### 3.5 Indexing conventions

- Fortran / MATLAB: **1-based** arrays; first debt point is zero debt.
- Python: **0-based**; `b=0` is zero debt.
- Mid income index: Fortran `CEILING(ySz/2)` → Python `(y_sz + 1) // 2 - 1`.
- Plot slice for choice probabilities: MATLAB index 150 → Python `b_fix_ix = 149`.

---

## 4. Port 2: Python/Numba → JAX (`long_term_debt_model_jax.py`)

### 4.1 Goals

- Keep **the same model and parameters** as the Numba port (and thus Fortran).
- Rewrite the computational core in **pure functional array code** that XLA can compile.
- Enable future **GPU/TPU** runs without a separate code path: install CUDA-enabled JAX and the same script uses the visible accelerator.
- Leave `long_term_debt_model.py` and `` unchanged.

### 4.2 Stack

| Layer | Choice | Role |
|-------|--------|------|
| Arrays / AD | `jax.numpy` | Device arrays (float64 enabled) |
| Compile | `jax.jit` | Fuse VFI step and simulation |
| Control flow | `jax.lax.while_loop`, `jax.lax.scan` | Stationary dist; time path |
| Random | `jax.random` | Simulation uniforms |
| SciPy-on-JAX | `jax.scipy.stats.norm`, `logsumexp` | Transitions; soft-max |
| Host edge | NumPy + Matplotlib | Moments (`ddof=1`), PDF plots |

Install (CPU example):

```bash
pip install -U "jax[cpu]"
```

GPU: use the platform-specific wheels from the [JAX installation docs](https://docs.jax.dev/en/latest/installation.html).

### 4.3 What changed relative to Numba

| Concern | Numba | JAX |
|---------|--------|-----|
| Parallel VFI | `prange` triple loops | Broadcast + matmul |
| Mutation | In-place array writes | Pure functions returning new arrays |
| Soft-max | Manual `exp` / `sum` | `logsumexp` (same math) |
| Simulation RNG | `np.random` inside `@njit` | `jax.random` + `lax.scan` |
| Discrete choice draw | Sequential CDF loop | `cumsum` + `searchsorted` |
| Convergence loop | Python `while` | Python `while` over jitted `vfi_step` (host prints) |
| Float width | float64 | float64 (`jax_enable_x64`) |
| Plot PDF name | `long_term_debt_model_plots.pdf` | `long_term_debt_model_jax_plots.pdf` |

### 4.4 Vectorized kernels (algebraic form)

**Default value update**

```text
vd1 = uhy + beta * (y_pi @ (gamm * v0[:, 0] + (1 - gamm) * vd0))
```

**Repayment / policies** (shapes \(Y\times B\times B\) for \(W\))

```text
ev = y_pi @ v0
cc[y,b,bp] = y[y] - kappa*b[b] + q0[y,bp]*(b'[bp] - (1-delta)*b[b])
w = where(cc > 0, u(cc) + beta * ev[:,None,:], very_negative)
b_pol, vr  = softmax / logsumexp over bp with temperature rho_b
d_pol, v1  = binary soft-max between vd1 and vr with rho_d
```

Peak temporary storage for \(W\) is about \(31\times 600\times 600\times 8\) bytes ≈ **89 MB** in float64—acceptable on modern CPUs and GPUs.

**Prices**

```text
qbar = sum_bpp b_pol[yp, bp, bpp] * q0[yp, bpp]   # (Y, B)
q1   = (y_pi @ ((1 - d_pol) * (kappa + (1 - delta) * qbar))) / (1 + rf)
```

**Expected \(B'\)**

```text
eb_pr = b_pol @ b_grid
```

**Simulation**  
`lax.scan` over \(t=1,\ldots,T-1\) with carry `(y, b, b', d, rng_key)`. Branching (re-entry, default, market vs exclusion) is written with `jnp.where` so the step is pure and jittable. Path dependence prevents parallelizing over time; GPU still helps the VFI tensors heavily.

### 4.5 Intentional non-parity

- **RNG streams differ.** Moments from simulation are statistically comparable, not bit-identical across Fortran / Numba / JAX.
- **Reduction order.** Matmul / reduction floating-point order can differ slightly from serial Fortran loops; with float64 and the same tolerances, VFI still converges to the same criteria.
- **No Fortran binary I/O.** Python/JAX keep results in memory and write plots/PDF only (no `results.mat` unless added later).

### 4.6 Full-grid run (reference)

A full-grid JAX solve on CPU (JAX 0.10.x, `Y=31`, `B=600`) was executed successfully in this project:

| Item | Result |
|------|--------|
| Backend | CPU (`jax` default) |
| Convergence | **428** iterations |
| Final `errV` | \(\approx 9.91\times 10^{-7}\) |
| Final `errQ` | \(\approx 4.35\times 10^{-12}\) |
| Wall time (order of magnitude) | ~1–2 minutes on the development machine (includes compile) |

Example moments from that run (valid sample; simulation RNG = JAX):

| Moment | Value (×100 where MATLAB prints percent-style) |
|--------|-----------------------------------------------:|
| Mean Debt to GDP | 7.88 |
| Mean Spread | 2.10 |
| Std Spread | 0.82 |
| Std log C | 1.72 |
| Std log GDP | 1.50 |
| Corr Sp, GDP | −43.93 |
| Corr TB/GDP, GDP | −29.31 |

Treat these as a **regression snapshot**, not a published calibration table: re-runs with different hardware or JAX versions may shift sim moments slightly.

---

## 5. How to run

### 5.1 Numba Python

```bash
python long_term_debt_model.py
```

Dependencies: `numpy`, `scipy`, `numba`, `matplotlib`.

### 5.2 JAX Python

```bash
pip install -U "jax[cpu]"   # or CUDA wheels for GPU
python long_term_debt_model_jax.py
```

Dependencies: `jax`, `jaxlib`, `numpy`, `matplotlib` (and SciPy transitively via JAX).

### 5.3 Smoke test

In either file, set:

```python
DEBUG_SMALL = True
```

This shrinks grids and iteration/sim caps for a fast end-to-end check. Restore `False` for the research-scale grid.

### 5.4 Fortran (reference)

See `Makefile` and `osc_slrm.sh` for the original build/run workflow (OpenMP).

---

## 6. Mapping of symbols to code arrays

| Theory | Fortran | Python / JAX dict or array |
|--------|---------|----------------------------|
| \(y\) grid | `yGrid` | `y_grid` / `eq["y_grid"]` |
| \(\Pi(y'|y)\) | `yPi` | `y_pi` / `eq["y_pi"]` |
| \(B\) grid | `bGrid` | `b_grid` / `eq["b_grid"]` |
| \(V\) | `V0`/`V1` | `V` / `eq["V"]` |
| \(V^r\) | `Vr` | `eq["Vr"]` |
| \(V^d\) | `Vd0`/`Vd1` | `eq["Vd"]` |
| \(q\) | `q0`/`q1` | `eq["q"]` |
| \(\Pr(d=1)\) | `dPol` | `eq["dPol"]` |
| \(\Pr(B')\) | `bPol` | `eq["bPol"]` |
| \(\mathbb{E}[B']\) | `EbPr` | `eq["EbPr"]` |

---

## 7. Design choices and future work

**Why taste shocks?**  
They smooth discrete choices, deliver closed-form policies, and make the VFI contraction better behaved than pure \(\arg\max\) with ties / bang-bang debt.

**Why long-term debt?**  
Duration (\(\delta,\kappa\)) makes spreads and incentives richer than one-period bonds: default risk is priced into a decaying coupon claim, and issuance is \(B'-(1-\delta)B\).

**Why JAX after Numba?**  
Numba is excellent for CPU loops that mirror Fortran. JAX rewrites the same math as **array programs** that:

- fuse into XLA kernels,
- run on GPU without OpenMP/CUDA C,
- compose with automatic differentiation if future estimation / GMM / gradient-based calibration is desired.

**Possible extensions (not implemented):**

- float32 / mixed precision for larger grids on GPU  
- `lax.while_loop` entirely on device (no host print each 10 iters)  
- Bit-reproducible cross-backend RNG  
- Export of equilibrium arrays to `.npz` / `.mat` for comparison with Fortran binaries  
- Automatic differentiation through VFI for sensitivity or estimation  

---

## 8. Authorship and provenance

- **Model write-up:** Gabriel Mihalache, *Long-Term Debt Model* (September 2022), `longTermDebtModel.tex` / `README.md`.
- **Original solver:** Fortran 90 + OpenMP in `src/`, with MATLAB post-processing in `results/`.
- **Python (NumPy/Numba) port:** `long_term_debt_model.py` — line-by-line economic parity with Fortran VFI and simulation, moments/plots aligned to the MATLAB scripts.
- **JAX port:** `long_term_debt_model_jax.py` — same model, vectorized/jitted for accelerator readiness; documented in this file.

---

## 9. Quick checklist: “are the ports the same model?”

| Check | Fortran | Numba | JAX |
|-------|:-------:|:-----:|:---:|
| Same parameters (`defMod`) | ✓ | ✓ | ✓ |
| Same \(V,V^d,V^r,q\) recursion | ✓ | ✓ | ✓ |
| Taste shocks \(\rho_D,\rho_B\) | ✓ | ✓ | ✓ |
| Long-term debt \(q,\delta,\kappa\) | ✓ | ✓ | ✓ |
| Same grid sizes (default) | 31×600 | 31×600 | 31×600 |
| Same valid-sample moments definition | MATLAB | ✓ | ✓ |
| Same diagnostic plot set | MATLAB | ✓ | ✓ |
| Identical sim RNG path | — | ≠ Fortran | ≠ Numba/Fortran |
| OpenMP / GPU acceleration | OpenMP | Numba `prange` | XLA (CPU/GPU) |

If VFI errors print below \(10^{-6}\) and policy/price plots look like the MATLAB `simplePlots` figures, the port is operating as the intended long-term debt default model.
