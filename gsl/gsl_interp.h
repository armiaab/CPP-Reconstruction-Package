#pragma once
#include <stddef.h>

typedef struct { int id; const char *name; } gsl_interp_type;
typedef struct { size_t cache; double x_prev; } gsl_interp_accel;

#ifdef __cplusplus
extern "C" {
#endif
extern const gsl_interp_type *gsl_interp_linear;
extern const gsl_interp_type *gsl_interp_cspline;
extern const gsl_interp_type *gsl_interp_steffen;

gsl_interp_accel *gsl_interp_accel_alloc(void);
void              gsl_interp_accel_free(gsl_interp_accel *a);
void              gsl_interp_accel_reset(gsl_interp_accel *a);
#ifdef __cplusplus
}
#endif
