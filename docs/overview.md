# CPP_Reconstruction — Project Overview

A hybrid project combining a C++ performance backend with the full original [MATLAB reconstruction package](https://github.com/mshoja/MATLAB-reconstruction-package) for reconstructing Langevin dynamics parameters from observed time-series data.

---

## What It Does

Given a discrete-time trajectory `{x_0, x_1, …, x_n}` sampled from a stochastic process, this library estimates the **drift** function μ(x) and **diffusion** function σ(x) of the underlying Itô SDE:

```
dX_t = μ(X_t) dt  +  σ(X_t) dW_t
```

Both functions are represented as smooth splines whose coefficients are fitted by maximising the Euler–Maruyama log-likelihood.

---

## Architecture

The project has two layers:

| Layer | What it covers | Implementation |
|-------|---------------|----------------|
| **C++ backend** | `euler_reconstruction` (spline mode) | `langevin_core.cpp` + MEX wrapper |
| **MATLAB layer** | All other functions — `hermite_reconstruction`, `legitimate_points`, `plot_results`, `eulergrad`, diagnostic tools, etc. | Original MATLAB package files, copied verbatim |

`euler_reconstruction.m` is the junction point: it is a drop-in replacement for the original with an identical calling convention, but internally calls the C++ MEX binary instead of pure MATLAB.

> **Limitation:** The C++ backend supports **spline mode only** (`nKnots` / `spline` options). Passing `'mu'` or `'sigma'` function handles raises an error. See [matlab_interface.md](matlab_interface.md) for details.

---

## Repository Layout

```
CPP_Reconstruction/
│
├── ── C++ backend ──────────────────────────────────────────────────────────
├── CMakeLists.txt          Build configuration (C++17, Eigen3, NLopt)
├── build.bat               MSVC convenience build script
├── README.md               Quick-start build instructions
│
├── langevin_core.h         Public C API — types, structs, function declarations
├── langevin_core.cpp       Full C++ implementation
├── langevin_mex.cpp        MATLAB MEX gateway
├── gsl_impl.cpp            Bundled mini-GSL (splines + quasi-random numbers)
├── langevin_mex.mexw64     Pre-compiled MEX binary (Windows x64)
│
├── gsl/
│   ├── gsl_errno.h
│   ├── gsl_interp.h
│   ├── gsl_spline.h
│   └── gsl_qrng.h
│
├── ── MATLAB entry point ───────────────────────────────────────────────────
├── euler_reconstruction.m  C++ MEX wrapper (spline mode; drop-in replacement)
│
├── ── Original MATLAB functions ────────────────────────────────────────────
├── hermite_reconstruction.m    Hermite-based reconstruction for low-resolution data
├── legitimate_points.m         Find valid initial points for Hermite reconstruction
├── plot_results.m              Visualise drift/diffusion estimates
├── simulate.m                  Forward simulation of fitted SDE
│
├── eulergrad.m / eulergrad_spline.m     Analytic gradients for parametric fits
├── eulerhess.m / eulerhess_spline.m     Analytic Hessians
├── gradest.m / hessian.m / hessdiag.m  Numerical differentiation utilities
├── derivest.m                  General derivative estimator
├── fixpars.m                   Parameter constraint helper
├── lagtable.m                  Lag-table construction
│
├── ME_TimeScale.m              Markov-Einstein time-scale estimation
├── RelaxationTime.m            Data resolution / relaxation time
├── armasel_s.m / armodel.m     AR model selection (standalone wrappers)
│
├── Coefficients_Spline.m       Pre-tabulated spline coefficient lookup
├── DataStandardizationForSplines.m
├── OrganizeData.m
├── prepare_replicateData.m
├── sample_stratified.m         Stratified mini-batch sampling
├── AllFigures.m                Reproduce all paper figures
├── test_reconst.m              Comprehensive test script
├── MyMatlabFunction.m / g_LL.m / Uncertainty.m
│
├── +coefficients/              MATLAB package — pre-tabulated spline coefficients
│   ├── +ETA/                   ETA coefficient functions
│   └── +EZ/                    EZ coefficient functions
│
├── Burg/                       ARMASA toolbox (required by ME_TimeScale)
│   ├── armasel.m               Main AR model selection function
│   ├── armasel_s.m             Multi-series variant
│   └── ... (50+ support files)
│
├── ── Examples ─────────────────────────────────────────────────────────────
└── Examples/
    ├── example_1.m   OU model — parametric & spline fits        [spline parts only]
    ├── example_2.m   Nonlinear May model — parametric & spline  [spline parts only]
    ├── example_3.m   Cyanobacteria data — spline reconstruction
    ├── example_4.m   Low-resolution OU — Hermite reconstruction [spline parts only]
    ├── example_5.m   Low-resolution OU spline + Hermite
    ├── example_6.m   Additive-noise spline model + Hermite
    ├── example_7.m   May model standardised + Hermite           [spline parts only]
    ├── example_8.m   May model spline + Hermite
    ├── example_9.m   NGRIP climate record — spline + Hermite
    └── example_10.m  Large-data mini-batch strategy
│
├── ── Data files ───────────────────────────────────────────────────────────
├── OUdata1D.mat / MayData1D.mat / MayData1D_Replicate.mat / MayData_LeftTail.mat
├── legpoints.mat
├── BGA_stdlevel_2011.csv / NGRIP20.csv
└── Mendota2015_Phyco.csv / Mendota2018_Phyco.csv / New_Cyanobacteria.csv
```

> Examples marked `[spline parts only]` contain both parametric-mode calls (which will error with the C++ backend) and spline-mode calls (which work). Run each section independently or skip the parametric blocks.

---

## External Dependencies

| Dependency | Version | Role |
|------------|---------|------|
| **Eigen3** | 3.4.0 | Dense linear algebra (basis matrices, vector operations) |
| **NLopt** | 2.7.1 | Nonlinear optimisation (SLSQP, L-BFGS, ISRES) |
| **ARMASA / Burg** | bundled | AR model selection — used by `ME_TimeScale` |

Eigen3 and NLopt are downloaded automatically by CMake via `FetchContent`. ARMASA is bundled in the `Burg/` directory.

---

## Key Concepts at a Glance

| Concept | Description |
|---------|-------------|
| **Euler–Maruyama likelihood** | Transitions modelled as Gaussians; NLL minimised over spline coefficients |
| **Spline parametrisation** | Drift and diffusion each represented by a spline with `n_mu` / `n_sig` knots |
| **Multi-start optimisation** | Sobol quasi-random starts explore the parameter space before local refinement |
| **Analytic gradients** | Pre-computed basis matrices allow O(n · k) gradient evaluation for linear/cubic splines |
| **Hermite reconstruction** | Higher-order method for low-resolution data; requires Euler result as warm start |
| **Boundary constraints** | Hard `[L, R]` limits on the state space; out-of-range observations are excluded |

---

## Document Index

| File | Contents |
|------|----------|
| [euler_reconstruction.md](euler_reconstruction.md) | Algorithm derivation and step-by-step walkthrough |
| [api_reference.md](api_reference.md) | Full C/C++ public API (types, structs, functions) |
| [internals.md](internals.md) | Implementation details: PPoly, basis matrices, optimisation loop |
| [splines.md](splines.md) | Supported interpolation methods and when to use each |
| [matlab_interface.md](matlab_interface.md) | MEX function usage and all options |
| [building.md](building.md) | CMake build instructions for Windows and Linux |
