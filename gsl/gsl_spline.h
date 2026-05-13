#pragma once
#include "gsl_interp.h"
#include "gsl_errno.h"

typedef struct {
    const gsl_interp_type *type;
    size_t size;
    double *xa;
    double *ya;
    double *dy;
} gsl_spline;

#ifdef __cplusplus
extern "C" {
#endif
gsl_spline *gsl_spline_alloc(const gsl_interp_type *type, size_t n);
int         gsl_spline_init(gsl_spline *spl, const double *xa, const double *ya, size_t n);
double      gsl_spline_eval(const gsl_spline *spl, double x, gsl_interp_accel *acc);
double      gsl_spline_eval_deriv(const gsl_spline *spl, double x, gsl_interp_accel *acc);
void        gsl_spline_free(gsl_spline *spl);
#ifdef __cplusplus
}
#endif
