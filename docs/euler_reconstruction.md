# Euler Reconstruction — Algorithm Deep Dive

> This document describes the algorithm implemented in the C++ backend (`langevin_core.cpp`).  It covers spline mode only.  For the full MATLAB interface — including parametric mode, Hermite reconstruction, and the other package functions — see [matlab_interface.md](matlab_interface.md).

## Background

The library solves the **inverse problem** for a continuous-time diffusion process.  We assume the true process satisfies the Itô SDE

```
dX_t = μ(x) dt  +  σ(x) dW_t
```

and that we observe it at equally-spaced times with step `dt`, yielding the sequence `{x_0, x_1, …, x_{n-1}}`.  The goal is to recover μ and σ from those observations.

---

## 1. Euler–Maruyama Discretisation

The Euler–Maruyama scheme approximates one time step as:

```
X_{k+1} = X_k + μ(X_k) · dt  +  σ(X_k) · √dt · Z_k,    Z_k ~ N(0,1)
```

Conditioned on `X_k = x_k` this means `X_{k+1}` is Gaussian:

```
X_{k+1} | X_k  ~  N( x_k + μ(x_k)·dt,  σ(x_k)²·dt )
```

The transition log-density for a single step is therefore

```
log p(x_{k+1}|x_k) = −½ · z_k²  −  log(σ_k · √dt)  −  ½ · log(2π)
```

where the **standardised residual** is

```
z_k = (x_{k+1} − x_k − μ_k · dt) / (σ_k · √dt)
```

and `μ_k = μ(x_k)`, `σ_k = σ(x_k)`.

---

## 2. Negative Log-Likelihood

Summing over all `n-1` consecutive pairs and dropping the constant `½ log(2π)`:

```
NLL(β_μ, β_σ) = Σ_k [ ½ · z_k²  +  log(σ_k)  +  ½ · log(dt) ]
```

Minimising NLL over the spline coefficients `(β_μ, β_σ)` gives the **maximum-likelihood estimate**.

---

## 3. Spline Parametrisation

Rather than estimating μ and σ pointwise, each function is represented by a spline on a fixed grid of **knots**:

```
μ(x; β_μ) = Σ_j  β_μ[j] · φ_j(x)          (n_mu basis functions)
σ(x; β_σ) = Σ_j  β_σ[j] · ψ_j(x)          (n_sig basis functions)
```

The supported basis types are linear (tent functions), cubic spline, PCHIP/Steffen, and additive (constant).  See [splines.md](splines.md) for details.

This reduces the optimisation to a finite-dimensional problem of dimension `n_mu + n_sig`.

---

## 4. Analytic Gradients (Linear & Cubic Splines)

For linear and cubic splines the basis functions are known analytically.  The library pre-computes two **basis matrices** before optimisation begins:

```
B_μ  ∈ R^{(n-1) × n_mu},    B_μ[k, j] = φ_j(x_k)
B_σ  ∈ R^{(n-1) × n_sig},   B_σ[k, j] = ψ_j(x_k)
```

Evaluating the splines at all data points then reduces to matrix–vector products:

```
μ_vec = B_μ · β_μ
σ_vec = B_σ · β_σ
```

The gradient of NLL with respect to `β_μ` and `β_σ` follows by the chain rule:

```
∂NLL/∂β_μ[j] = Σ_k  (−z_k / (σ_k√dt)) · B_μ[k, j]
             = B_μᵀ · (−z / (σ·√dt))

∂NLL/∂β_σ[j] = Σ_k  (−z_k² + 1) / σ_k · B_σ[k, j]
             = B_σᵀ · ((1 − z²) / σ)
```

Both require only matrix–vector multiplications — O(n · k) work instead of O(n · k²).

For **PCHIP splines** the basis functions are not separable, so gradients are computed by **central differences** with `h = 1×10⁻⁶`.

---

## 5. Optimisation Strategy

### 5.1 Multi-Start

A single gradient-based run is likely to converge to a local minimum.  The library therefore:

1. Generates `n_starts` starting points using a **Sobol quasi-random sequence** (implemented as Halton sequences over the first 64 primes).
2. Each candidate is scaled to the parameter bounds `[lb, ub]`.
3. Any candidate where σ(x) ≤ 0 somewhere on the data range is rejected and resampled (up to 100 tries per start).
4. The full NLopt optimiser is run from every accepted starting point.
5. The solution with the lowest NLL is kept.

### 5.2 Solvers

Three NLopt algorithms are available:

| ID | NLopt Algorithm | Type | When to Use |
|----|----------------|------|-------------|
| `LC_SLSQP` | `LD_SLSQP` | Local, gradient-based, constrained | Default; fast and robust |
| `LC_LBFGS` | `LD_LBFGS` | Local, quasi-Newton, unconstrained | When bounds are not active |
| `LC_ISRES` | `GN_ISRES` | Global, evolutionary, derivative-free | Exploratory; very slow |

### 5.3 Stopping Criteria

Optimisation stops when either:
- The number of objective evaluations reaches `max_iter` (default 10 000), or
- The change in objective value falls below `ftol` (default 1×10⁻⁹).

---

## 6. Positivity Constraint on σ

The diffusion must remain strictly positive for the likelihood to be well-defined.  After fitting, and also before each optimisation start, the library evaluates σ on a fine uniform grid across the data range and rejects any parameter set that gives σ ≤ 0 anywhere.

---

## 7. Data Preprocessing

Before optimisation:

1. **NaN removal** — Any `x_k` equal to `NaN` is dropped.
2. **Range filtering** — If `L` and/or `R` are provided, observations outside `[L, R]` are discarded.
3. **Pair extraction** — Consecutive pairs `(x_k, x_{k+1})` are formed from the cleaned sequence.  The `n-1` pairs form the input to the likelihood.

---

## 8. Automatic Knot and Bound Selection

When `knots_mu` / `knots_sig` are `NULL`, the library places `n_mu` / `n_sig` knots uniformly across the data range with a 2 % margin:

```
margin = 0.02 * (max(data) - min(data))
knots  = linspace(min(data) - margin,  max(data) + margin,  n_knots)
```

When `lb` / `ub` are `NULL`, default bounds are used:

| Parameter | Lower | Upper |
|-----------|-------|-------|
| Drift β_μ | −10 | +10 |
| Diffusion β_σ | 0.01 | 10 |
| Diffusion β_σ (additive) | 0.05 | 10 |

---

## 9. End-to-End Data Flow

```
Input:  {x_0 … x_{n-1}},  dt,  n_mu,  n_sig,  sp_mu,  sp_sig,  solver,  n_starts
         │
         ▼
_preprocess()          Filter NaN & out-of-range; form (x_k, x_{k+1}) pairs
         │
         ▼
_auto_knots()          Place knots uniformly across cleaned data range
_auto_bounds()         Set [lb, ub] for (β_μ, β_σ)
         │
         ▼
Prob constructor       Build B_μ, B_σ basis matrices (once, reused every eval)
         │
         ▼
Sobol sequence ──► for each start 1…n_starts:
         │              scale to [lb, ub]
         │              check σ > 0  (retry up to 100×)
         │              run NLopt ──► eval_with_grad() called repeatedly
         │                              μ = B_μ · β_μ
         │                              σ = B_σ · β_σ
         │                              z = (Δx − μ·dt) / (σ·√dt)
         │                              NLL = Σ(½z² + log σ)
         │                              grad via Bᵀ multiplication
         │              keep if NLL < best so far
         ▼
Return best (β_μ, β_σ) and NLL
```

---

## 10. Simulation (Forward Direction)

Given estimated parameters, `euler_simulate()` draws a trajectory by iterating:

```
x[t] = x[t-1]  +  μ(x[t-1]) · dt  +  σ(x[t-1]) · √dt · Z_t
```

where `Z_t ~ N(0,1)` is generated via the **Box–Muller transform** applied to `std::mt19937` uniform samples.  If `L` and/or `R` are finite, the state is reflected at the boundaries.
