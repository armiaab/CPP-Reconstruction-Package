# MATLAB Interface

## Available Functions

The project root contains the full original MATLAB reconstruction package alongside the C++ backend.  Add the project root (and `Burg/`) to the MATLAB path once:

```matlab
addpath('/path/to/CPP_Reconstruction');
addpath('/path/to/CPP_Reconstruction/Burg');
```

| Function | Source | Description |
|----------|--------|-------------|
| `euler_reconstruction` | C++ wrapper | Spline-mode Euler reconstruction — this repo's main contribution |
| `hermite_reconstruction` | Original MATLAB | Hermite method for low-resolution / sparse data |
| `legitimate_points` | Original MATLAB | Find valid Hermite starting points from an Euler result |
| `plot_results` | Original MATLAB | Visualise fitted drift and diffusion |
| `simulate` | Original MATLAB | Forward Euler–Maruyama simulation |
| `eulergrad` / `eulergrad_spline` | Original MATLAB | Analytic gradients for parametric Euler fits |
| `eulerhess` / `eulerhess_spline` | Original MATLAB | Analytic Hessians |
| `ME_TimeScale` | Original MATLAB | Markov-Einstein time-scale (requires `Burg/` on path) |
| `RelaxationTime` | Original MATLAB | Data resolution diagnostic |
| `sample_stratified` | Original MATLAB | Stratified mini-batch sampling for large datasets |
| `AllFigures` / `test_reconst` | Original MATLAB | Reproduce paper figures / comprehensive test script |

The `Examples/` folder contains ten runnable scripts that demonstrate the full workflow.

---

## `euler_reconstruction` — C++ backed wrapper

| File | Role |
|------|------|
| `euler_reconstruction.m` | **Drop-in replacement** for the original package — same function name, same options, same result struct |
| `langevin_mex` (compiled) | C++ MEX backend called by the wrapper — can also be called directly |

> **Parametric mode is not supported.**  Passing `'mu'` or `'sigma'` function handles raises an error.  Use `'nKnots'` / `'spline'` (spline mode) instead.  `'gradient_fun'` is silently ignored.

To use this C++ backend in place of the original MATLAB `euler_reconstruction`, add the project root to the MATLAB path — no other code changes are needed.

---

## `euler_reconstruction` (MATLAB wrapper — use this)

This is the primary entry point.  Its calling convention is **identical** to the original `euler_reconstruction.m` from the MATLAB reconstruction package.

---

## Compilation

```matlab
mex -O CXXFLAGS='$CXXFLAGS -std=c++17' ...
    -I/path/to/eigen3 ...
    langevin_mex.cpp langevin_core.cpp gsl_impl.cpp ...
    -lnlopt -o langevin_mex
```

Refer to the project `README.md` for platform-specific build commands.

---

## `euler_reconstruction` usage

```matlab
result = euler_reconstruction(data, dt, 'Name', Value, ...)
```

### Required inputs

| Argument | Type | Description |
|----------|------|-------------|
| `data` | `double` vector | Observed time series (may contain `NaN`) |
| `dt` | `double` scalar | Time step between observations |

### Name–value options (all keys case-insensitive)

| Name | Type | Default | Description |
|------|------|---------|-------------|
| `'nKnots'` | scalar or `[n_mu, n_sig]` | `[8 8]` | Number of knots for drift and diffusion. |
| `'knots'` | double vector | — | Explicit flat knot vector: first `n_mu` → mu knots, next `n_sig` → sigma knots. Overrides `'nKnots'`. |
| `'spline'` | 2-char string | `'CC'` | Spline type per function. `'L'`=linear, `'C'`=cubic, `'Q'`=order-3 (mapped to cubic), `'P'`=PCHIP. |
| `'L'` | double | `[]` | Left boundary: observations below `L` are discarded. |
| `'R'` | double | `[]` | Right boundary: observations above `R` are discarded. |
| `'lb'` | double vector | auto | Lower bounds for all parameters, length `n_mu + n_sig`. |
| `'ub'` | double vector | auto | Upper bounds for all parameters, length `n_mu + n_sig`. |
| `'solver'` | string | `'fmincon'` | `'fmincon'`/`'PS'` → SLSQP (gradient, local); `'GWO'` → ISRES (global, slow); `'lbfgs'`; `'isres'`. |
| `'search_agents'` | integer | `1` | Number of random starts. Increase for robustness. |
| `'maxiter'` | integer | `10000` | Max function evaluations per start. Pass `realmax` for unlimited. |
| `'ftol'` | double | `1e-9` | Function-value convergence tolerance. |
| `'useparallel'` | — | — | Silently ignored (C++ backend is single-threaded). |
| `'reconst_fraction'` | — | — | Silently ignored. |
| `'gradient_fun'`, `'hessian_fun'` | — | — | Silently ignored. |
| `'mu'`, `'sigma'` | — | — | **Raises an error** — parametric mode is not supported by the C++ backend. |

### Result struct

All fields match the original MATLAB package exactly.

| Field | Type | Description |
|-------|------|-------------|
| `estimated_par` | `1 × (n_mu+n_sig)` | Fitted spline coefficients |
| `negloglik` | scalar | Negative log-likelihood at optimum |
| `datarange` | `[1 × 2]` | `[min max]` of observations after NaN/boundary removal |
| `dt` | scalar | Time step (echoed) |
| `mufun` | `@(x, par)` | Drift evaluation handle |
| `sigmafun` | `@(x, par)` | Diffusion evaluation handle |
| `computation_time` | scalar | Wall-clock seconds |
| `knots` | `{1×2}` cell | `{knots_mu, knots_sig}` |
| `options` | struct | Parsed options (`solver`, `spline`, `l`, `r`, `lb`, `ub`, `nknots`, `search_agents`, `maxiter`) |

### Examples

Runnable scripts are in the `Examples/` folder.  Each script adds the project root and `Burg/` to the MATLAB path automatically via `mfilename`.

**Spline fit (example_5 pattern — QQ, 5 starts, Hermite follow-up):**
```matlab
% Low-resolution data — shrink state space first
L = -2.5;  R = 2.5;
data(data < L | data > R) = nan;
nmu = 7;  nsigma = 7;

result_euler = euler_reconstruction(data, dt, ...
    'nKnots', [nmu nsigma], ...
    'spline', 'QQ', ...
    'L', L, 'R', R, ...
    'lb', [zeros(1,nmu)-10, zeros(1,nsigma)+eps], ...
    'ub',  zeros(1,nmu+nsigma)+10, ...
    'solver', 'fmincon', 'search_agents', 5);

legpts = legitimate_points(data, dt, 'prev', result_euler, ...
    'prev_range', 0.5, 'j', 3, 'k', 9);
result_hermite = hermite_reconstruction(data, dt, 'prev', legpts, ...
    'solver', 'fmincon');

xplot = linspace(L, R, 2000);
plot_results(result_euler,  xplot);
plot_results(result_hermite, xplot);
```

**High-resolution spline fit (example_1 pattern — CC, boundary-trimmed):**
```matlab
L = -2;  R = 1.8;
data(data < L | data > R) = nan;
nmu = 8;  nsigma = 8;

result = euler_reconstruction(data, dt, ...
    'nKnots', [nmu nsigma], ...
    'spline', 'CC', ...
    'L', L, 'R', R, ...
    'lb', [zeros(1,nmu)-10, zeros(1,nsigma)+eps], ...
    'ub',  zeros(1,nmu+nsigma)+10, ...
    'solver', 'fmincon');

plot_results(result, linspace(L, R, 500));
```

**Additive (constant) diffusion:**
```matlab
nmu = 8;  nsigma = 1;   % nsigma == 1 forces constant diffusion

result = euler_reconstruction(data, dt, ...
    'nKnots', [nmu nsigma], ...
    'spline', 'QQ', ...
    'L', L, 'R', R, ...
    'lb', [zeros(1,nmu)-10, eps], ...
    'ub',  [ones(1,nmu)*10, 10], ...
    'solver', 'fmincon', 'search_agents', 5);
```

---

## `langevin_mex` (direct MEX — low-level)

`langevin_mex` can be called directly when you need precise control over knot positions or want raw numeric outputs without a result struct.

### Reconstruction mode

```matlab
[par, nll, time_sec, n_evals] = langevin_mex(data, dt, Name, Value, ...)
```

Key options: `'knots_mu'`, `'knots_sig'`, `'spline'` (`'cc'`/`'ll'`/`'pp'`/`'ca'` etc.), `'solver'` (`'slsqp'`/`'lbfgs'`/`'isres'`/`'gwo'`), `'nstarts'` or `'search_agents'`, `'maxiter'`, `'ftol'`, `'l'`, `'r'`, `'lb'`, `'ub'`.

### Simulation mode

```matlab
traj = langevin_mex('simulate', par, knots_mu, knots_sig, dt, x0, T, Name, Value, ...)
```

Options: `'spline'`, `'l'`, `'r'`, `'seed'`.

### Example

```matlab
n_mu = 6;  n_sig = 4;
knots_mu  = linspace(min(data)-0.02*range(data), max(data)+0.02*range(data), n_mu);
knots_sig = linspace(min(data)-0.02*range(data), max(data)+0.02*range(data), n_sig);

[par, nll] = langevin_mex(data, dt, ...
    'knots_mu', knots_mu, 'knots_sig', knots_sig, ...
    'spline', 'cc', 'nstarts', 10);

traj = langevin_mex('simulate', par, knots_mu, knots_sig, dt, data(1), 5000, ...
    'spline', 'cc', 'seed', 42);
```
