#include "mex.h"
#include "langevin_core.h"
#include <cstring>
#include <cmath>
#include <string>

static std::string _str(const mxArray *a) {
    if (!mxIsChar(a)) return "";
    char buf[128]; mxGetString(a, buf, sizeof(buf));
    for (char &c : buf[0] == '\0' ? buf : buf) if (!c) break; else c = (char)tolower(c);
    return std::string(buf);
}

static const mxArray *_opt(int nrhs, const mxArray **prhs, int i0, const char *key)
{
    std::string k(key);
    for (int i = i0; i < nrhs - 1; i += 2)
        if (mxIsChar(prhs[i]) && _str(prhs[i]) == k) return prhs[i + 1];
    return nullptr;
}

static double _optd(int nr, const mxArray **pr, int i0, const char *k, double def) {
    auto *v = _opt(nr, pr, i0, k); return v ? mxGetScalar(v) : def;
}
static int _opti(int nr, const mxArray **pr, int i0, const char *k, int def) {
    auto *v = _opt(nr, pr, i0, k); return v ? (int)mxGetScalar(v) : def;
}
static std::string _opts(int nr, const mxArray **pr, int i0, const char *k, const char *def) {
    auto *v = _opt(nr, pr, i0, k);
    if (!v) return def;
    char buf[64]; mxGetString(v, buf, sizeof(buf)); return buf;
}

static LC_SplineType _sp(char c) {
    switch (tolower(c)) {
        case 'l': return LC_LINEAR;
        case 'c': return LC_CUBIC;
        case 'q': return LC_QUADRATIC;   /* piecewise-quadratic, matches MATLAB 'Q' */
        case 'p': return LC_PCHIP;
        default:  return LC_CUBIC;
    }
}
static LC_Solver _slv(const std::string &s) {
    if (s == "lbfgs")                         return LC_LBFGS;
    if (s == "isres" || s == "gwo")           return LC_ISRES;
    // fmincon / slsqp / ps / anything else → SLSQP
    return LC_SLSQP;
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    bool sim = nrhs >= 1 && mxIsChar(prhs[0]) && _str(prhs[0]) == "simulate";

    if (!sim) {
        if (nrhs < 2)
            mexErrMsgIdAndTxt("lc:usage",
                "[par,nll,t,nevals] = langevin_mex(data, dt, [opts])");
        if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]))
            mexErrMsgIdAndTxt("lc:type", "data must be a real double vector.");

        const double *data = mxGetPr(prhs[0]);
        int    n    = (int)mxGetNumberOfElements(prhs[0]);
        double dt   = mxGetScalar(prhs[1]);
        int    opt0 = 2;

        int n_mu = 8, n_sig = 8;
        auto *nk = _opt(nrhs, prhs, opt0, "nknots");
        if (nk) {
            if (mxGetNumberOfElements(nk) == 1) n_mu = n_sig = (int)mxGetScalar(nk);
            else { n_mu = (int)mxGetPr(nk)[0]; n_sig = (int)mxGetPr(nk)[1]; }
        }

        std::string sp_s = _opts(nrhs, prhs, opt0, "spline", "cc");
        LC_SplineType sp_mu  = _sp(sp_s[0]);
        LC_SplineType sp_sig = sp_s.size() > 1 ? _sp(sp_s[1]) : sp_mu;
        if (n_sig == 1) sp_sig = LC_ADDITIVE;

        LC_Solver solver = _slv(_opts(nrhs, prhs, opt0, "solver", "slsqp"));

        // 'search_agents' is the original package name; 'nstarts' is our alias
        int nstarts = _opti(nrhs, prhs, opt0, "nstarts", -1);
        if (nstarts < 0) nstarts = _opti(nrhs, prhs, opt0, "search_agents", 1);
        int    maxiter = _opti(nrhs, prhs, opt0, "maxiter", 10000);
        double ftol    = _optd(nrhs, prhs, opt0, "ftol",    1e-9);
        double L       = _optd(nrhs, prhs, opt0, "l",       (double)NAN);
        double R       = _optd(nrhs, prhs, opt0, "r",       (double)NAN);

        auto *mk = _opt(nrhs, prhs, opt0, "knots_mu");
        auto *sk = _opt(nrhs, prhs, opt0, "knots_sig");
        auto *ml = _opt(nrhs, prhs, opt0, "lb");
        auto *mu = _opt(nrhs, prhs, opt0, "ub");
        const double *km = mk ? mxGetPr(mk) : nullptr;
        const double *ks = sk ? mxGetPr(sk) : nullptr;
        const double *lb = ml ? mxGetPr(ml) : nullptr;
        const double *ub = mu ? mxGetPr(mu) : nullptr;
        if (mk) n_mu  = (int)mxGetNumberOfElements(mk);
        if (sk) n_sig = (int)mxGetNumberOfElements(sk);

        int    dim = n_mu + n_sig;
        double *par_buf = new double[dim]();
        double  nll = 0, tsec = 0; int evals = 0;

        euler_reconstruct_flat(data, n, dt,
            km, n_mu, ks, n_sig, lb, ub, L, R,
            (int)sp_mu, (int)sp_sig, (int)solver,
            nstarts, maxiter, ftol,
            par_buf, &nll, &tsec, &evals);

        if (nlhs >= 1) {
            plhs[0] = mxCreateDoubleMatrix(1, dim, mxREAL);
            memcpy(mxGetPr(plhs[0]), par_buf, dim * sizeof(double));
        }
        if (nlhs >= 2) plhs[1] = mxCreateDoubleScalar(nll);
        if (nlhs >= 3) plhs[2] = mxCreateDoubleScalar(tsec);
        if (nlhs >= 4) plhs[3] = mxCreateDoubleScalar((double)evals);
        delete[] par_buf;
        return;
    }

    if (nrhs < 7)
        mexErrMsgIdAndTxt("lc:usage",
            "traj = langevin_mex('simulate',par,km,ks,dt,x0,T,[opts])");

    const double *par = mxGetPr(prhs[1]);
    const double *km  = mxGetPr(prhs[2]); int n_mu  = (int)mxGetNumberOfElements(prhs[2]);
    const double *ks  = mxGetPr(prhs[3]); int n_sig = (int)mxGetNumberOfElements(prhs[3]);
    double dt  = mxGetScalar(prhs[4]);
    double x0  = mxGetScalar(prhs[5]);
    int    T   = (int)mxGetScalar(prhs[6]);
    int    opt0 = 7;

    std::string sp_s = _opts(nrhs, prhs, opt0, "spline", "cc");
    LC_SplineType sp_mu  = _sp(sp_s[0]);
    LC_SplineType sp_sig = sp_s.size() > 1 ? _sp(sp_s[1]) : sp_mu;
    if (n_sig == 1) sp_sig = LC_ADDITIVE;

    double L    = _optd(nrhs, prhs, opt0, "l", (double)NAN);
    double R    = _optd(nrhs, prhs, opt0, "r", (double)NAN);
    long   seed = (long)_optd(nrhs, prhs, opt0, "seed", 42.0);

    mxArray *traj = mxCreateDoubleMatrix(T, 1, mxREAL);
    euler_simulate(par, km, n_mu, ks, n_sig,
                   sp_mu, sp_sig, dt, x0, L, R, T, seed,
                   mxGetPr(traj));
    if (nlhs >= 1) plhs[0] = traj;
}
