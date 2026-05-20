# C/C++ API Reference

All public symbols are declared in `langevin_core.h` and compiled with `extern "C"` linkage so the library is callable from plain C as well as C++.

---

## Enumerations

### `LC_SplineType`

Controls how drift and diffusion functions are interpolated between knots.

```c
typedef enum {
    LC_LINEAR   = 0,   // Piecewise-linear (tent functions)
    LC_CUBIC    = 1,   // Natural cubic spline (C² continuous)
    LC_PCHIP    = 2,   // Steffen/PCHIP — shape-preserving, no overshoots
    LC_ADDITIVE = 3    // Single constant value (1 parameter regardless of n_knots)
} LC_SplineType;
```

See [splines.md](splines.md) for detailed guidance on choosing a type.

---

### `LC_Solver`

Selects the NLopt optimisation algorithm.

```c
typedef enum {
    LC_SLSQP = 0,   // Sequential Least Squares (local, gradient, constrained) — default
    LC_LBFGS = 1,   // Limited-memory BFGS (local, gradient, unconstrained)
    LC_ISRES = 2    // Improved Stochastic Ranking ES (global, derivative-free)
} LC_Solver;
```

---

## Structs

### `LC_Result`

Returned by `euler_reconstruct()`. Heap-allocated; free with `lc_result_free()`.

```c
typedef struct {
    double *par;        // Fitted parameters: [β_μ₀…β_μ_{n_mu-1}, β_σ₀…β_σ_{n_sig-1}]
    int     n_par;      // Total number of parameters (n_mu + n_sig)
    double  negloglik;  // Negative log-likelihood at the optimum
    double  time_sec;   // Wall-clock time in seconds
    int     n_evals;    // Number of objective function evaluations (total, all starts)
    int     n_starts;   // Number of optimisation starts actually run
    int     status;     // NLopt return code (positive = success, see nlopt.h)
} LC_Result;
```

---

## Functions

### `euler_reconstruct`

Fit drift and diffusion splines to an observed time series.

```c
LC_Result *euler_reconstruct(
    const double *data,        // [n]     Input time series (may contain NaN)
    int           n,           //         Length of data
    double        dt,          //         Time step between observations

    const double *knots_mu,    // [n_mu]  Knot positions for drift spline  (NULL = auto)
    int           n_mu,        //         Number of drift knots / coefficients
    const double *knots_sig,   // [n_sig] Knot positions for diffusion spline (NULL = auto)
    int           n_sig,       //         Number of diffusion knots / coefficients

    const double *lb,          // [n_mu+n_sig] Lower bounds (NULL = auto)
    const double *ub,          // [n_mu+n_sig] Upper bounds (NULL = auto)

    double        L,           //         Left data boundary (NaN = no constraint)
    double        R,           //         Right data boundary (NaN = no constraint)

    LC_SplineType sp_mu,       //         Interpolation type for drift
    LC_SplineType sp_sig,      //         Interpolation type for diffusion

    LC_Solver     solver,      //         Optimisation algorithm
    int           n_starts,    //         Number of random starts (≥ 1)
    int           max_iter,    //         Max objective evaluations per start
    double        ftol         //         Function-value convergence tolerance
);
```

**Returns** a heap-allocated `LC_Result*` on success, or `NULL` on allocation failure.  The caller must call `lc_result_free()` when done.

**Parameter notes:**
- When `knots_mu` / `knots_sig` are `NULL`, knots are placed uniformly with a 2 % margin around the observed data range.
- When `lb` / `ub` are `NULL`, default bounds `[-10, 10]` (drift) and `[0.01, 10]` (diffusion) are used.
- `L` / `R` of `NaN` mean no boundary clipping.
- `n_starts = 1` with `LC_SLSQP` is fast but may find a local minimum; increase for robustness.

---

### `euler_reconstruct_flat`

Flat C interface — writes results into caller-allocated arrays instead of returning a struct.  Used internally by the MATLAB MEX gateway.

```c
int euler_reconstruct_flat(
    const double *data,  int n,  double dt,
    const double *knots_mu,  int n_mu,
    const double *knots_sig, int n_sig,
    const double *lb,    const double *ub,
    double L,  double R,
    LC_SplineType sp_mu,  LC_SplineType sp_sig,
    LC_Solver solver,
    int n_starts,  int max_iter,  double ftol,

    /* outputs — must point to pre-allocated buffers */
    double *par_out,        // [n_mu + n_sig] fitted parameters
    double *negloglik_out,  // scalar
    double *time_out,       // scalar (seconds)
    int    *nevals_out      // scalar
);
```

**Returns** the NLopt status code (positive = success, negative = error).

---

### `euler_simulate`

Simulate a trajectory forward in time given fitted (or user-supplied) parameters.

```c
void euler_simulate(
    const double *par,         // [n_mu + n_sig]  Spline coefficients
    const double *knots_mu,    // [n_mu]          Drift knots
    int           n_mu,        //                 Number of drift knots
    const double *knots_sig,   // [n_sig]         Diffusion knots
    int           n_sig,       //                 Number of diffusion knots

    LC_SplineType sp_mu,       //                 Drift interpolation type
    LC_SplineType sp_sig,      //                 Diffusion interpolation type

    double        dt,          //                 Time step
    double        x0,          //                 Initial state x[0]
    double        L,           //                 Left boundary (NaN = none)
    double        R,           //                 Right boundary (NaN = none)
    int           T,           //                 Number of steps to simulate
    long          seed,        //                 RNG seed (use ≥ 0)

    double       *out          // [T]             Output trajectory (pre-allocated)
);
```

Generates `T` states `{x[0], x[1], …, x[T-1]}` using the Euler–Maruyama step with `std::mt19937` and Box–Muller normal sampling.

---

### `lc_result_free`

Releases memory allocated by `euler_reconstruct()`.

```c
void lc_result_free(LC_Result *r);
```

Always call this after you are done with an `LC_Result`.

---

## Typical Usage (C++)

```cpp
#include "langevin_core.h"
#include <cmath>

int main() {
    // --- Fit ---
    double data[1000] = { /* your observations */ };
    int    n           = 1000;
    double dt          = 0.01;

    LC_Result *res = euler_reconstruct(
        data, n, dt,
        NULL, 5,   // 5 drift knots, placed automatically
        NULL, 5,   // 5 diffusion knots
        NULL, NULL,
        NAN, NAN,  // no boundary constraints
        LC_CUBIC, LC_CUBIC,
        LC_SLSQP,
        10,        // 10 random starts
        10000, 1e-9
    );

    // --- Inspect ---
    printf("NLL = %f  (%d evals)\n", res->negloglik, res->n_evals);
    for (int i = 0; i < res->n_par; i++)
        printf("par[%d] = %f\n", i, res->par[i]);

    // --- Simulate ---
    double traj[2000];
    euler_simulate(res->par,
                   NULL, 5,   // knots must match reconstruction
                   NULL, 5,
                   LC_CUBIC, LC_CUBIC,
                   dt, data[0], NAN, NAN,
                   2000, 42,
                   traj);

    lc_result_free(res);
    return 0;
}
```

> **Note:** When using `NULL` knots in `euler_simulate`, you must pass the same knot count used during `euler_reconstruct`, and ensure the knots are consistent.  In practice, save the knots returned by reconstruction and pass them explicitly to simulation.
