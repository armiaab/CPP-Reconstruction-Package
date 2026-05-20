# Spline Types

The drift μ(x) and diffusion σ(x) are each represented as a piecewise polynomial (spline) whose coefficients are optimised.  The library supports four interpolation types, selectable independently for each function.

---

## Overview

| `LC_SplineType` | Name | Smoothness | Parameters | Gradients | Notes |
|----------------|------|-----------|------------|-----------|-------|
| `LC_LINEAR` (0) | Piecewise linear | C⁰ (continuous) | `n_knots` | Analytic | Fastest; coarse |
| `LC_CUBIC` (1) | Natural cubic spline | C² (twice differentiable) | `n_knots` | Analytic | Best default choice |
| `LC_PCHIP` (2) | Steffen / PCHIP | C¹ (once differentiable) | `n_knots` | Numeric | Shape-preserving |
| `LC_ADDITIVE` (3) | Constant | C∞ | 1 (always) | Analytic | Special case |

---

## `LC_LINEAR` — Piecewise Linear

Each interval `[x_i, x_{i+1}]` is a straight line:

```
μ(x) = μ_i  +  (μ_{i+1} − μ_i) · (x − x_i) / (x_{i+1} − x_i)
```

The `n_knots` coefficients are simply the function values at the knots.

**Pros:** Simple; very fast; analytic gradients.
**Cons:** Kinks at every knot — not physically realistic for smooth drift/diffusion.  Needs more knots to achieve the same accuracy as cubic.

**Use when:** You want the simplest possible model, or are doing a quick sanity check.

---

## `LC_CUBIC` — Natural Cubic Spline

A classic C² spline: on each interval, a cubic polynomial, with the constraint that the first and second derivatives match at interior knots.  The "natural" boundary conditions set the second derivative to zero at the endpoints.

The second derivatives at the knots are found by solving the tridiagonal system (see [internals.md](internals.md) §8).  The `n_knots` free parameters are the function values at the knots; the spline is fully determined from these.

**Pros:** Smooth (C²); analytic gradients; accurate with few knots; default choice.
**Cons:** Can oscillate (Runge phenomenon) with too many knots or unevenly-spaced data.

**Use when:** In doubt.  This is the right choice for most physical systems.

---

## `LC_PCHIP` — Steffen / PCHIP

Piecewise Cubic Hermite Interpolation (PCHIP) using the Steffen monotone method.  The derivatives at each knot are chosen to preserve the shape (monotonicity, local extrema) of the data.

Unlike a natural cubic spline, PCHIP does **not** minimise global curvature, but it also does **not** overshoot between data points.

**Gradient note:** The basis functions for PCHIP are not easily separable, so gradients are computed by central finite differences (h = 10⁻⁶) rather than analytically.  This makes each evaluation ~2×n_par times slower.

**Pros:** No spurious oscillations; good for monotone or near-monotone functions.
**Cons:** Numeric gradients — slower optimisation; not C².

**Use when:** The drift or diffusion is expected to be monotone and you want a smooth fit that doesn't oscillate.

---

## `LC_ADDITIVE` — Constant

A degenerate case where the spline reduces to a single constant value regardless of the number of knots specified:

```
μ(x) = β₀   (same everywhere)
```

Only one parameter is optimised, even if `n_mu > 1`.

**Use when:** You want to test the hypothesis that drift or diffusion is constant (e.g., additive noise σ(x) = const).

---

## Choosing Knot Count

The number of knots controls how flexibly the spline can adapt to the data.  More knots → more flexibility, but also more parameters to optimise and greater risk of overfitting.

Typical starting values:

| System type | Suggested `n_mu` | Suggested `n_sig` |
|-------------|-----------------|------------------|
| Linear drift, constant diffusion | 2–3 | 1 (additive) |
| Nonlinear drift, smooth diffusion | 4–8 | 3–5 |
| Complex multi-well potential | 8–15 | 4–8 |

When knots are not supplied explicitly, they are placed **uniformly** across the data range with a 2 % margin (see [euler_reconstruction.md § 8](euler_reconstruction.md)).

---

## Mixed-Type Examples

```c
// Cubic drift, constant diffusion (very common for Ornstein-Uhlenbeck-like systems)
euler_reconstruct(data, n, dt, NULL, 5, NULL, 1,
                  NULL, NULL, NAN, NAN,
                  LC_CUBIC, LC_ADDITIVE, ...);

// PCHIP drift (expected to be monotone), cubic diffusion
euler_reconstruct(data, n, dt, NULL, 6, NULL, 4,
                  NULL, NULL, NAN, NAN,
                  LC_PCHIP, LC_CUBIC, ...);
```
