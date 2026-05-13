#ifndef LANGEVIN_CORE_H
#define LANGEVIN_CORE_H

#if defined(_WIN32) || defined(__CYGWIN__)
#  ifdef LANGEVIN_BUILD_DLL
#    define LC_API __declspec(dllexport)
#  elif defined(LANGEVIN_STATIC)
#    define LC_API
#  else
#    define LC_API __declspec(dllimport)
#  endif
#else
#  define LC_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    LC_LINEAR   = 0,
    LC_CUBIC    = 1,
    LC_PCHIP    = 2,
    LC_ADDITIVE = 3
} LC_SplineType;

typedef enum {
    LC_SLSQP = 0,
    LC_LBFGS = 1,
    LC_ISRES = 2
} LC_Solver;

typedef struct {
    double *par;
    int     n_par;
    double  negloglik;
    double  time_sec;
    int     n_evals;
    int     n_starts;
    int     status;
} LC_Result;

void lc_result_free(LC_Result *r);

LC_Result *euler_reconstruct(
    const double *data,       int n,
    double        dt,
    const double *knots_mu,   int n_mu,
    const double *knots_sig,  int n_sig,
    const double *lb,
    const double *ub,
    double        L,          double R,
    LC_SplineType sp_mu,
    LC_SplineType sp_sig,
    LC_Solver     solver,
    int           n_starts,
    int           max_iter,
    double        ftol);

void euler_simulate(
    const double *par,
    const double *knots_mu,   int n_mu,
    const double *knots_sig,  int n_sig,
    LC_SplineType sp_mu,
    LC_SplineType sp_sig,
    double dt, double x0, double L, double R,
    int T, long seed,
    double *out);

int euler_reconstruct_flat(
    const double *data,      int n,       double dt,
    const double *knots_mu,  int n_mu,
    const double *knots_sig, int n_sig,
    const double *lb,
    const double *ub,
    double        L,         double R,
    int           sp_mu,     int sp_sig,
    int           solver,
    int           n_starts,  int max_iter, double ftol,
    double       *par_out,
    double       *negloglik_out,
    double       *time_out,
    int          *n_evals_out);

#ifdef __cplusplus
}
#endif
#endif
