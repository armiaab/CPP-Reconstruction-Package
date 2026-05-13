#pragma once
#define GSL_SUCCESS 0

#ifdef __cplusplus
extern "C" {
#endif
typedef void (*gsl_error_handler_t)(const char*, const char*, int, int);
gsl_error_handler_t gsl_set_error_handler(gsl_error_handler_t h);
#ifdef __cplusplus
}
#endif
