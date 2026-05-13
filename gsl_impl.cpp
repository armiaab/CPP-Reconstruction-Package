#ifdef _MSC_VER
#  define _USE_MATH_DEFINES
#endif
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <vector>

#include "gsl/gsl_errno.h"
#include "gsl/gsl_interp.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_qrng.h"

extern "C" gsl_error_handler_t gsl_set_error_handler(gsl_error_handler_t h)
{ (void)h; return nullptr; }

static const gsl_interp_type _linear  = {0, "linear"};
static const gsl_interp_type _cspline = {1, "cspline"};
static const gsl_interp_type _steffen = {2, "steffen"};

const gsl_interp_type *gsl_interp_linear  = &_linear;
const gsl_interp_type *gsl_interp_cspline = &_cspline;
const gsl_interp_type *gsl_interp_steffen = &_steffen;

extern "C" gsl_interp_accel *gsl_interp_accel_alloc(void) {
    gsl_interp_accel *a = (gsl_interp_accel*)std::malloc(sizeof(gsl_interp_accel));
    if (a) { a->cache = 0; a->x_prev = 0.0; }
    return a;
}
extern "C" void gsl_interp_accel_free(gsl_interp_accel *a) { std::free(a); }
extern "C" void gsl_interp_accel_reset(gsl_interp_accel *a) {
    if (a) { a->cache = 0; a->x_prev = 0.0; }
}

static void cspline_derivs(const double *x, const double *y,
                            size_t n, double *dy)
{
    if (n == 0) return;
    if (n == 1) { dy[0] = 0.0; return; }
    if (n == 2) {
        dy[0] = dy[1] = (y[1]-y[0]) / (x[1]-x[0]);
        return;
    }

    size_t m = n - 1;
    std::vector<double> h(m), s(m);
    for (size_t i = 0; i < m; ++i) {
        h[i] = x[i+1] - x[i];
        s[i] = (y[i+1] - y[i]) / h[i];
    }

    size_t ns = n - 2;
    std::vector<double> M(n, 0.0);

    if (ns > 0) {
        std::vector<double> diag(ns), sup(ns, 0.0), sub(ns, 0.0), rhs(ns);
        for (size_t i = 0; i < ns; ++i) {
            diag[i] = 2.0 * (h[i] + h[i+1]);
            if (i < ns-1) sup[i] = h[i+1];
            if (i > 0)    sub[i] = h[i];
            rhs[i]  = 6.0 * (s[i+1] - s[i]);
        }

        std::vector<double> c2(ns, 0.0), d2(ns);
        c2[0] = sup[0] / diag[0];
        d2[0] = rhs[0] / diag[0];
        for (size_t i = 1; i < ns; ++i) {
            double w = diag[i] - sub[i] * c2[i-1];
            c2[i] = (i < ns-1) ? sup[i] / w : 0.0;
            d2[i] = (rhs[i] - sub[i] * d2[i-1]) / w;
        }

        M[ns] = d2[ns-1];
        for (int i = (int)ns-2; i >= 0; --i)
            M[i+1] = d2[i] - c2[i] * M[i+2];
    }

    for (size_t i = 0; i < m; ++i)
        dy[i] = s[i] - h[i] * (2.0*M[i] + M[i+1]) / 6.0;
    dy[n-1] = s[m-1] + h[m-1] * (M[m-1] + 2.0*M[n-1]) / 6.0;
}

static void steffen_derivs(const double *x, const double *y,
                            size_t n, double *dy)
{
    if (n == 0) return;
    if (n == 1) { dy[0] = 0.0; return; }
    if (n == 2) {
        dy[0] = dy[1] = (y[1]-y[0])/(x[1]-x[0]);
        return;
    }

    size_t m = n - 1;
    std::vector<double> h(m), s(m);
    for (size_t i = 0; i < m; ++i) {
        h[i] = x[i+1] - x[i];
        s[i] = (y[i+1] - y[i]) / h[i];
    }

    for (size_t i = 1; i < n-1; ++i) {
        double pi = (s[i-1]*h[i] + s[i]*h[i-1]) / (h[i-1] + h[i]);
        if (s[i-1] * s[i] <= 0.0) {
            dy[i] = 0.0;
        } else {
            double lim = 2.0 * std::min(std::abs(s[i-1]), std::abs(s[i]));
            dy[i] = std::copysign(std::min(std::abs(pi), lim), pi);
        }
    }

    auto boundary = [&](size_t k, size_t k2, double *out) {
        double p = ((2.0*h[k] + h[k2])*s[k] - h[k]*s[k2]) / (h[k] + h[k2]);
        if (std::copysign(1.0, p) != std::copysign(1.0, s[k])) p = 0.0;
        else if (std::copysign(1.0, s[k]) != std::copysign(1.0, s[k2])
                 && std::abs(p) > 3.0*std::abs(s[k]))
            p = 3.0 * s[k];
        *out = p;
    };
    boundary(0, 1, &dy[0]);
    boundary(m-1, m-2, &dy[n-1]);
}

static void linear_derivs(const double *x, const double *y,
                           size_t n, double *dy)
{
    if (n == 0) return;
    if (n == 1) { dy[0] = 0.0; return; }
    for (size_t i = 0; i < n-1; ++i)
        dy[i] = (y[i+1]-y[i]) / (x[i+1]-x[i]);
    dy[n-1] = dy[n-2];
}

extern "C" {

gsl_spline *gsl_spline_alloc(const gsl_interp_type *type, size_t n) {
    gsl_spline *spl = (gsl_spline*)std::malloc(sizeof(gsl_spline));
    if (!spl) return nullptr;
    spl->type = type;
    spl->size = n;
    spl->xa = (double*)std::malloc(n * sizeof(double));
    spl->ya = (double*)std::malloc(n * sizeof(double));
    spl->dy = (double*)std::malloc(n * sizeof(double));
    if (!spl->xa || !spl->ya || !spl->dy) {
        std::free(spl->xa); std::free(spl->ya); std::free(spl->dy);
        std::free(spl);
        return nullptr;
    }
    return spl;
}

int gsl_spline_init(gsl_spline *spl,
                    const double *xa, const double *ya, size_t n)
{
    spl->size = n;
    std::memcpy(spl->xa, xa, n * sizeof(double));
    std::memcpy(spl->ya, ya, n * sizeof(double));
    switch (spl->type->id) {
        case 1:  cspline_derivs(xa, ya, n, spl->dy); break;
        case 2:  steffen_derivs(xa, ya, n, spl->dy); break;
        default: linear_derivs (xa, ya, n, spl->dy); break;
    }
    return GSL_SUCCESS;
}

double gsl_spline_eval(const gsl_spline *spl, double x, gsl_interp_accel *)
{
    size_t n = spl->size;
    if (n == 0) return 0.0;
    if (n == 1) return spl->ya[0];

    if (x <= spl->xa[0])   return spl->ya[0];
    if (x >= spl->xa[n-1]) return spl->ya[n-1];

    size_t lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        size_t mid = (lo + hi) / 2;
        if (spl->xa[mid] <= x) lo = mid; else hi = mid;
    }
    size_t k = lo;
    double h  = spl->xa[k+1] - spl->xa[k];
    double t  = x - spl->xa[k];

    if (spl->type->id == 0)
        return spl->ya[k] + t * spl->dy[k];

    double y0 = spl->ya[k],  y1 = spl->ya[k+1];
    double d0 = spl->dy[k],  d1 = spl->dy[k+1];
    double c  = (3.0*(y1-y0)/h - 2.0*d0 - d1) / h;
    double d  = (2.0*(y0-y1)/h +     d0 + d1) / (h*h);
    return y0 + t*(d0 + t*(c + t*d));
}

double gsl_spline_eval_deriv(const gsl_spline *spl, double x, gsl_interp_accel *)
{
    size_t n = spl->size;
    if (n == 0) return 0.0;
    if (n == 1) return spl->dy[0];

    size_t lo = 0, hi = n - 1;
    while (hi - lo > 1) {
        size_t mid = (lo + hi) / 2;
        if (spl->xa[mid] <= x) lo = mid; else hi = mid;
    }
    if (hi < n &&
        std::abs(x - spl->xa[hi]) < std::abs(x - spl->xa[lo]))
        return spl->dy[hi];
    return spl->dy[lo];
}

void gsl_spline_free(gsl_spline *spl) {
    if (!spl) return;
    std::free(spl->xa); std::free(spl->ya); std::free(spl->dy);
    std::free(spl);
}

}

static const gsl_qrng_type _sobol_type = {0, "halton"};
const gsl_qrng_type *gsl_qrng_sobol = &_sobol_type;

static const int kPrimes[64] = {
     2,  3,  5,  7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
    53, 59, 61, 67, 71, 73, 79, 83, 89, 97,101,103,107,109,113,
   127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,
   199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,
   283,293,307,311
};

static double halton_val(unsigned n, int base) {
    double result = 0.0, f = 1.0;
    while (n > 0) {
        f      /= base;
        result += f * (int)(n % (unsigned)base);
        n      /= (unsigned)base;
    }
    return result;
}

extern "C" {

gsl_qrng *gsl_qrng_alloc(const gsl_qrng_type *, unsigned dim) {
    gsl_qrng *q = (gsl_qrng*)std::malloc(sizeof(gsl_qrng));
    if (q) { q->dim = dim; q->seq = 1; }
    return q;
}

int gsl_qrng_get(gsl_qrng *q, double *x) {
    unsigned n = q->seq++;
    unsigned maxd = (unsigned)(sizeof(kPrimes)/sizeof(kPrimes[0]));
    for (unsigned d = 0; d < q->dim; ++d)
        x[d] = halton_val(n, kPrimes[d < maxd ? d : maxd-1]);
    return GSL_SUCCESS;
}

void gsl_qrng_free(gsl_qrng *q) { std::free(q); }

}
