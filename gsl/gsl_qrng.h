#pragma once
#include <stddef.h>
#include "gsl_errno.h"

typedef struct { int id; const char *name; } gsl_qrng_type;

typedef struct {
    unsigned dim;
    unsigned seq;
} gsl_qrng;

#ifdef __cplusplus
extern "C" {
#endif
extern const gsl_qrng_type *gsl_qrng_sobol;

gsl_qrng *gsl_qrng_alloc(const gsl_qrng_type *type, unsigned dim);
int       gsl_qrng_get  (gsl_qrng *q, double *x);
void      gsl_qrng_free (gsl_qrng *q);
#ifdef __cplusplus
}
#endif
