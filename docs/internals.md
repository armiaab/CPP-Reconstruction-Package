# Implementation Internals

> This document covers the C++ backend only (`langevin_core.cpp` and `gsl_impl.cpp`).  The MATLAB layer — `euler_reconstruction.m`, `hermite_reconstruction.m`, and the rest of the original package — is separate; see [matlab_interface.md](matlab_interface.md).

---

## 1. The `PPoly` Struct

`PPoly` is the runtime representation of a fitted spline.  It stores the piecewise-polynomial coefficients so that evaluation reduces to a single cubic polynomial per query — no GSL calls needed at inference time.

```cpp
struct PPoly {
    std::vector<double> x;           // Knot positions (m+1 values)
    std::vector<double> a, b, c, d;  // Polynomial coefficients (m segments)
    size_t              m;           // Number of polynomial segments
    LC_SplineType       sp;          // Which interpolation type
    double              scalar;      // For LC_ADDITIVE: the single constant value
};
```

On the i-th segment `[x[i], x[i+1]]`, the polynomial is evaluated as:

```
h = x_query − x[i]
p = a[i]  +  b[i]·h  +  c[i]·h²  +  d[i]·h³
```

For `LC_LINEAR`, `c` and `d` are always zero.  For `LC_ADDITIVE`, `b`, `c`, `d` are zero and `a[0]` holds the single parameter value.

---

## 2. Building `PPoly` — `_build_pp()`

Converts a `gsl_spline` object (which stores GSL-native derivatives) into a `PPoly`.

Steps:
1. Allocate vectors of length `m` (number of segments = number of knots − 1).
2. For each segment `i`:
   - `a[i]` = spline value at left knot `x[i]`
   - Retrieve the first and second derivative from GSL at `x[i]`
   - Back-compute `b`, `c`, `d` from the Hermite conditions at the two endpoints
3. For `LC_ADDITIVE`, skip all of this and store the single scalar.

---

## 3. Evaluating `PPoly` — `_eval_pp()`

```cpp
void _eval_pp(const PPoly &pp, const double *xs, int n, double *out);
```

For each query point `xs[k]`:
1. **Binary search** the knot array for the enclosing interval.
2. Compute `h = xs[k] − x[left_knot]`.
3. Evaluate the cubic using Horner's rule: `a + h*(b + h*(c + h*d))`.
4. If the query is outside `[x[0], x[m]]`, write `NaN` — signals an invalid evaluation.

The binary search is O(log m), which is negligible compared to the likelihood computation for typical knot counts (≤ 20).

---

## 4. Basis Matrices — `_build_basis_matrix()`

For analytic gradients the library builds one **basis matrix** per spline (drift and diffusion) once, before any optimisation iterations.

```
B ∈ R^{n_pairs × n_knots}
B[k, j] = value of the j-th basis function at x_k
```

Each column j corresponds to setting `β[j] = 1` and all other coefficients to zero, and evaluating the spline at all data points.  Because spline evaluation is linear in the coefficients:

```
spline_values = B · β
```

This makes every evaluation inside the optimiser a simple matrix–vector multiply, avoiding repeated spline reconstruction.

---

## 5. The `Prob` Struct

`Prob` owns everything the optimiser needs during its iterations:

```cpp
struct Prob {
    // Data
    const double *x0, *x1;   // Consecutive pairs (x_k, x_{k+1})
    int           n_pairs;

    double dt, sqrt_dt;

    // Knots
    std::vector<double> knots_mu, knots_sig;
    LC_SplineType        sp_mu, sp_sig;

    // Basis matrices (Eigen)
    Eigen::MatrixXd B_mu, B_sig;

    // Working buffers
    Eigen::VectorXd mu_vec, sig_vec, z_vec;

    // Best result tracking
    double          best_nll;
    std::vector<double> best_par;

    // Evaluation counter
    int n_evals;
};
```

### `Prob::eval_with_grad()`

Called by NLopt on every objective evaluation.

```
1.  Split parameter vector θ into (β_μ, β_σ)
2.  mu_vec  = B_mu  * β_μ
    sig_vec = B_sig * β_σ
3.  z_vec[k] = (x1[k] − x0[k] − mu_vec[k]*dt) / (sig_vec[k]*sqrt_dt)
4.  NLL = Σ_k ( 0.5*z[k]² + log(sig[k]) )  +  const
5.  If gradient requested:
      ∂NLL/∂β_μ  = B_muᵀ  · (−z / (sig * sqrt_dt))    [analytic]
      ∂NLL/∂β_σ  = B_sigᵀ · ((1 − z²) / sig)           [analytic]
      (or central differences for PCHIP)
6.  Update best_nll / best_par if this is a new minimum
7.  Return NLL
```

---

## 6. Optimisation Loop

Implemented in `euler_reconstruct()` after `Prob` is constructed:

```
sobol = Halton quasi-random generator (dim = n_par)
best_nll = +∞

for start = 1 … n_starts:
    tries = 0
    do:
        draw next Sobol point u ∈ [0,1]^{n_par}
        θ = lb + u * (ub − lb)          // scale to bounds
        tries++
    while σ(θ) ≤ 0 anywhere  AND  tries < 100

    if tries == 100: skip this start

    run NLopt with initial point θ
    if result.nll < best_nll:
        best_nll = result.nll
        best_par = result.par
```

---

## 7. NLopt Callback

NLopt calls a static `_obj()` function which immediately forwards to `Prob::eval_with_grad()`.  The `void*` user-data pointer carries the `Prob*`.

Bounds are passed to NLopt via `nlopt_set_lower_bounds` / `nlopt_set_upper_bounds` for the `LC_SLSQP` solver.  For `LC_LBFGS` (unconstrained), bounds are encoded as a penalty term instead.

---

## 8. GSL Mini-Library (`gsl_impl.cpp`)

A self-contained implementation of the GSL spline and quasi-random APIs to avoid a system GSL dependency.

### Spline Memory Layout

```c
struct gsl_spline {
    double *x, *y;      // Input knots and values
    double *dy, *d2y;   // First and second derivatives at knots (computed on init)
    size_t  size;       // Number of knots
    const gsl_interp_type *type;  // Pointer to method vtable
};
```

### Derivative Computation at Initialisation

| Method | Function | Algorithm |
|--------|----------|-----------|
| Linear | `linear_derivs()` | Forward differences; constant at boundaries |
| Cubic | `cspline_derivs()` | Solve tridiagonal system for M_i (Thomas algorithm) |
| PCHIP / Steffen | `steffen_derivs()` | Monotone slope limiting to prevent overshoots |

#### Thomas Algorithm (cubic spline)

Solves the `(n-2) × (n-2)` tridiagonal system arising from natural cubic spline conditions in O(n) time:

```
Forward sweep:   c'[i] = c[i] / (b[i] − a[i]·c'[i−1])
                 d'[i] = (d[i] − a[i]·d'[i−1]) / (b[i] − a[i]·c'[i−1])
Back substitution: M[n-2] = d'[n-2]
                   M[i]   = d'[i] − c'[i]·M[i+1]
```

where `a, b, c, d` are the standard tridiagonal entries from the second-derivative conditions.

### Evaluation

`gsl_spline_eval()` calls the appropriate method for the spline type and uses a `gsl_interp_accel` cache to remember the last interval — O(1) amortised for sequential queries.

### Quasi-Random Numbers — Halton Sequence

```c
gsl_qrng_alloc(gsl_qrng_halton, dim);
```

Returns a quasi-random generator using the **Halton sequence**: for dimension `d`, use the `d`-th prime `p_d` as the base:

```
H_n(p) = Σ_{k≥0}  floor(n / p^k) mod p  ·  p^{-(k+1)}
```

The first 64 primes `{2, 3, 5, 7, …, 311}` are hard-coded in the implementation.  Halton sequences have lower discrepancy than pseudo-random uniform samples, which means that `n_starts` optimisation starts cover the parameter space more evenly.

---

## 9. Key Numerical Details

| Detail | Value / Approach |
|--------|-----------------|
| Gradient finite-difference step | `h = 1×10⁻⁶` (PCHIP only) |
| σ positivity check grid | 200 uniform points across data range |
| σ validation threshold | `σ > 0` strictly |
| Out-of-range extrapolation | Returns `NaN` (caught by NLL computation) |
| Normal RNG for simulation | `std::mt19937` + Box–Muller |
| Eigen build flag | `EIGEN_NO_DEBUG` (no bounds checks in release) |
| SIMD target | AVX2 (MSVC `/arch:AVX2`, GCC `-march=native`) |
