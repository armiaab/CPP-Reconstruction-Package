#include "langevin_core.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <ctime>
#include <limits>
#include <memory>
#include <vector>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_spline.h>

#include <Eigen/Dense>

#include <nlopt.h>

using Eigen::ArrayXd;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Map;

static const double LOG_SQRT_2PI = 0.9189385332046728;
static const double DBL_INF      = std::numeric_limits<double>::max();
static const double NaN_         = std::numeric_limits<double>::quiet_NaN();

static void _gsl_noop(const char*, const char*, int, int) {}
static void _gsl_init() {
    static bool done = false;
    if (!done) { gsl_set_error_handler(_gsl_noop); done = true; }
}

static const gsl_interp_type *_gsl_itype(LC_SplineType sp) {
    switch (sp) {
        case LC_LINEAR: return gsl_interp_linear;
        case LC_CUBIC:  return gsl_interp_cspline;
        case LC_PCHIP:  return gsl_interp_steffen;
        default:        return gsl_interp_linear;
    }
}

struct PPoly {
    std::vector<double> x;
    std::vector<double> a, b, c, d;
    size_t        m  = 0;
    LC_SplineType sp = LC_LINEAR;
    double scalar    = 0.0;
};

static void _build_pp(gsl_spline *spl,
                       const double *kn, const double *vals,
                       size_t m, LC_SplineType sp, PPoly &pp)
{
    pp.m  = m;
    pp.sp = sp;
    pp.x.assign(kn, kn + m);
    pp.a.resize(m-1); pp.b.resize(m-1);
    pp.c.resize(m-1); pp.d.resize(m-1);

    if (sp == LC_LINEAR) {
        for (size_t i = 0; i < m-1; ++i) {
            double h = kn[i+1] - kn[i];
            pp.a[i] = vals[i];
            pp.b[i] = (vals[i+1] - vals[i]) / h;
            pp.c[i] = pp.d[i] = 0.0;
        }
        return;
    }

    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    std::vector<double> dy(m);
    for (size_t i = 0; i < m; ++i)
        dy[i] = gsl_spline_eval_deriv(spl, kn[i], acc);
    gsl_interp_accel_free(acc);

    for (size_t i = 0; i < m-1; ++i) {
        double h  = kn[i+1] - kn[i];
        double y0 = vals[i],  y1 = vals[i+1];
        double d0 = dy[i],    d1 = dy[i+1];
        pp.a[i] = y0;
        pp.b[i] = d0;
        pp.c[i] = (3.0*(y1-y0)/h - 2.0*d0 - d1) / h;
        pp.d[i] = (2.0*(y0-y1)/h +   d0 + d1) / (h*h);
    }
}

static void _eval_pp(const PPoly &pp, const double *xq, size_t n, double *out)
{
    if (pp.sp == LC_ADDITIVE) {
        std::fill(out, out + n, pp.scalar);
        return;
    }

    const double lo = pp.x.front(), hi = pp.x.back();
    const int    M  = (int)pp.m - 2;

    for (size_t i = 0; i < n; ++i) {
        double xi = xq[i];
        if (!std::isfinite(xi) || xi < lo || xi > hi) { out[i] = NaN_; continue; }
        auto it = std::upper_bound(pp.x.begin(), pp.x.end(), xi);
        int k = std::min((int)(it - pp.x.begin()) - 1, M);
        if (k < 0) k = 0;
        double h = xi - pp.x[k];
        out[i] = pp.a[k] + h * (pp.b[k] + h * (pp.c[k] + h * pp.d[k]));
    }
}

static MatrixXd _build_basis_matrix(
    const double *kn, int nk, LC_SplineType sp,
    const double *xq, size_t n)
{
    if (sp == LC_ADDITIVE)
        return MatrixXd::Ones((int)n, 1);

    MatrixXd B((int)n, nk);

    if (sp == LC_LINEAR) {
        B.setZero();
        for (size_t i = 0; i < n; ++i) {
            double xi = xq[i];
            if (!std::isfinite(xi) || xi < kn[0] || xi > kn[nk-1]) continue;
            int k = 0;
            while (k < nk - 2 && xi > kn[k+1]) ++k;
            double t = (xi - kn[k]) / (kn[k+1] - kn[k]);
            B((int)i, k)   = 1.0 - t;
            B((int)i, k+1) = t;
        }
        return B;
    }

    gsl_spline       *spl = gsl_spline_alloc(_gsl_itype(sp), nk);
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    std::vector<double> v(nk, 0.0);

    for (int k = 0; k < nk; ++k) {
        std::fill(v.begin(), v.end(), 0.0);
        v[k] = 1.0;
        gsl_spline_init(spl, kn, v.data(), nk);
        gsl_interp_accel_reset(acc);
        for (size_t i = 0; i < n; ++i) {
            double xi = xq[i];
            if (!std::isfinite(xi) || xi < kn[0] || xi > kn[nk-1])
                B((int)i, k) = 0.0;
            else
                B((int)i, k) = gsl_spline_eval(spl, xi, acc);
        }
    }

    gsl_interp_accel_free(acc);
    gsl_spline_free(spl);
    return B;
}

static double _euler_nll(const double *X,  const double *X0,
                           const double *mu, const double *sig,
                           double dt, size_t n)
{
    Map<const ArrayXd> Xm(X,n), X0m(X0,n), muv(mu,n), sigv(sig,n);

    ArrayXd sd = sigv * std::sqrt(dt);
    if ((sd <= 0.0).any()) return DBL_INF;

    ArrayXd z = (Xm - X0m - muv * dt) / sd;
    ArrayXd t = 0.5 * z.square() + sd.log() + LOG_SQRT_2PI;
    double  s = (t.isFinite().cast<double>() * t).sum();
    return std::isfinite(s) ? s : DBL_INF;
}

struct Prob {
    const double *X0, *X;
    size_t        nd;
    const double *km, *ks;
    size_t        nm, ns;
    double        dt;
    LC_SplineType sm, ss;

    MatrixXd B_mu, B_sig;
    bool     analytic_grad;
    double   nll_const;

    gsl_spline *spl_m = nullptr, *spl_s = nullptr;
    std::vector<double> mbuf, sbuf;
    PPoly pp_m, pp_s;

    ArrayXd mu_buf, sig_buf, res_buf, z_buf;

    int    n_evals = 0;
    double best_f  = DBL_INF;
    std::vector<double> best_par;

    Prob(const double *X0_, const double *X_, size_t nd_,
         const double *km_, size_t nm_,
         const double *ks_, size_t ns_,
         double dt_, LC_SplineType sm_, LC_SplineType ss_)
        : X0(X0_), X(X_), nd(nd_),
          km(km_), ks(ks_), nm(nm_), ns(ns_),
          dt(dt_), sm(sm_), ss(ss_),
          mbuf(nd_), sbuf(nd_), best_par(nm_ + ns_, 0.0),
          mu_buf((int)nd_), sig_buf((int)nd_),
          res_buf((int)nd_), z_buf((int)nd_)
    {
        nll_const = (double)nd_ * (0.5 * std::log(dt_) + LOG_SQRT_2PI);

        analytic_grad =
            (sm_ == LC_LINEAR || sm_ == LC_CUBIC || sm_ == LC_ADDITIVE) &&
            (ss_ == LC_LINEAR || ss_ == LC_CUBIC || ss_ == LC_ADDITIVE);

        if (analytic_grad) {
            B_mu  = _build_basis_matrix(km_, (int)nm_, sm_, X0_, nd_);
            B_sig = _build_basis_matrix(ks_, (int)ns_, ss_, X0_, nd_);
        } else {
            if (sm_ != LC_ADDITIVE) spl_m = gsl_spline_alloc(_gsl_itype(sm_), nm_);
            if (ss_ != LC_ADDITIVE) spl_s = gsl_spline_alloc(_gsl_itype(ss_), ns_);
        }
    }

    ~Prob() {
        if (spl_m) gsl_spline_free(spl_m);
        if (spl_s) gsl_spline_free(spl_s);
    }

    double eval_with_grad(const double *par, double *grad) {
        ++n_evals;
        const size_t dim = nm + ns;

        if (analytic_grad) {
            Map<const VectorXd> pm(par,      (int)nm);
            Map<const VectorXd> ps(par + nm, (int)ns);

            mu_buf  = (B_mu  * pm).array();
            sig_buf = (B_sig * ps).array();

            if (sig_buf.minCoeff() <= 0.0) {
                if (grad) std::fill(grad, grad + dim, 0.0);
                return DBL_INF;
            }

            Map<const ArrayXd> X0m(X0, (int)nd), X1m(X, (int)nd);
            const double sqdt = std::sqrt(dt);

            res_buf = X1m - X0m - mu_buf * dt;
            z_buf   = res_buf / (sig_buf * sqdt);

            double f = sig_buf.log().sum()
                     + 0.5 * z_buf.square().sum()
                     + nll_const;

            if (!std::isfinite(f)) {
                if (grad) std::fill(grad, grad + dim, 0.0);
                return DBL_INF;
            }

            if (grad) {
                ArrayXd dL_dmu  = -res_buf / sig_buf.square();
                ArrayXd dL_dsig = (1.0 - z_buf.square()) / sig_buf;

                Map<VectorXd>(grad,      (int)nm) =
                    B_mu.transpose()  * dL_dmu.matrix();
                Map<VectorXd>(grad + nm, (int)ns) =
                    B_sig.transpose() * dL_dsig.matrix();
            }

            if (f < best_f) {
                best_f = f;
                std::copy(par, par + dim, best_par.begin());
            }
            return f;

        } else {
            double f = _eval_core(par);
            if (f < best_f) {
                best_f = f;
                std::copy(par, par + dim, best_par.begin());
            }
            if (grad) {
                const double h = 1e-6;
                std::vector<double> ph(par, par + dim);
                for (size_t j = 0; j < dim; ++j) {
                    double orig = ph[j];
                    ph[j] = orig + h; double fp = _eval_core(ph.data()); n_evals++;
                    ph[j] = orig - h; double fm = _eval_core(ph.data()); n_evals++;
                    ph[j] = orig;
                    grad[j] = (std::isfinite(fp) && std::isfinite(fm))
                                ? (fp - fm) / (2.0 * h) : 0.0;
                }
            }
            return f;
        }
    }

private:
    double _eval_core(const double *par) {
        const double *pm = par, *ps = par + nm;

        if (sm == LC_ADDITIVE) { pp_m.sp = LC_ADDITIVE; pp_m.scalar = pm[0]; }
        else {
            if (gsl_spline_init(spl_m, km, pm, nm) != GSL_SUCCESS) return DBL_INF;
            _build_pp(spl_m, km, pm, nm, sm, pp_m);
        }
        if (ss == LC_ADDITIVE) { pp_s.sp = LC_ADDITIVE; pp_s.scalar = ps[0]; }
        else {
            if (gsl_spline_init(spl_s, ks, ps, ns) != GSL_SUCCESS) return DBL_INF;
            _build_pp(spl_s, ks, ps, ns, ss, pp_s);
        }

        _eval_pp(pp_m, X0, nd, mbuf.data());
        _eval_pp(pp_s, X0, nd, sbuf.data());
        return _euler_nll(X, X0, mbuf.data(), sbuf.data(), dt, nd);
    }
};

static double _obj(unsigned, const double *par, double *grad, void *ud)
{
    return static_cast<Prob*>(ud)->eval_with_grad(par, grad);
}

static int _run(Prob &p, LC_Solver solver,
                 const double *lb, const double *ub,
                 int dim, int maxeval, double ftol,
                 std::vector<double> &par)
{
    nlopt_algorithm alg;
    switch (solver) {
        case LC_LBFGS: alg = NLOPT_LD_LBFGS;  break;
        case LC_ISRES: alg = NLOPT_GN_ISRES;   break;
        default:       alg = NLOPT_LD_SLSQP;   break;
    }

    nlopt_opt opt = nlopt_create(alg, (unsigned)dim);
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);
    nlopt_set_min_objective(opt, _obj, &p);
    nlopt_set_maxeval(opt, maxeval);
    nlopt_set_ftol_abs(opt, ftol);
    nlopt_set_ftol_rel(opt, ftol);

    double fv = DBL_INF;
    int rc = nlopt_optimize(opt, par.data(), &fv);
    nlopt_destroy(opt);

    par = p.best_par;
    return rc;
}

static void _preprocess(const double *data, int n, double L, double R,
                          std::vector<double> &X0, std::vector<double> &X)
{
    X0.clear(); X.clear();
    X0.reserve(n); X.reserve(n);
    for (int i = 0; i < n-1; ++i) {
        double a = data[i], b = data[i+1];
        if (!std::isfinite(a) || !std::isfinite(b)) continue;
        if (std::isfinite(L) && (a < L || b < L))   continue;
        if (std::isfinite(R) && (a > R || b > R))   continue;
        X0.push_back(a);
        X.push_back(b);
    }
}

static std::vector<double> _auto_knots(const double *data, int n,
                                        double L, double R, int nk)
{
    double lo, hi;
    if (std::isfinite(L) && std::isfinite(R)) {
        lo = L; hi = R;
    } else {
        lo = 1e300; hi = -1e300;
        for (int i = 0; i < n; ++i) {
            if (!std::isfinite(data[i])) continue;
            if (data[i] < lo) lo = data[i];
            if (data[i] > hi) hi = data[i];
        }
    }
    double m = (hi - lo) * 0.02;
    lo -= m; hi += m;
    std::vector<double> k(nk);
    for (int i = 0; i < nk; ++i)
        k[i] = lo + (hi - lo) * i / (nk - 1.0);
    return k;
}

static void _auto_bounds(int nm, int ns, LC_SplineType ss,
                           std::vector<double> &lb, std::vector<double> &ub)
{
    lb.assign(nm + ns, 0.0);
    ub.assign(nm + ns, 0.0);
    std::fill(lb.begin(),      lb.begin() + nm, -10.0);
    std::fill(ub.begin(),      ub.begin() + nm,  10.0);
    double sig_lo = (ss == LC_ADDITIVE) ? 0.05 : 0.01;
    std::fill(lb.begin() + nm, lb.end(),  sig_lo);
    std::fill(ub.begin() + nm, ub.end(),  10.0);
}

static bool _sigma_ok(const double *par, size_t nm, size_t ns,
                       LC_SplineType ss, const double *ks,
                       const double *fine, size_t nf)
{
    if (ss == LC_ADDITIVE) return par[nm] > 0.0;

    gsl_spline *spl = gsl_spline_alloc(_gsl_itype(ss), ns);
    gsl_spline_init(spl, ks, par + nm, ns);
    PPoly pp;
    _build_pp(spl, ks, par + nm, ns, ss, pp);
    gsl_spline_free(spl);

    std::vector<double> out(nf);
    _eval_pp(pp, fine, nf, out.data());
    for (double v : out)
        if (!std::isfinite(v) || v <= 0.0) return false;
    return true;
}

extern "C"
LC_Result *euler_reconstruct(
    const double *data, int n, double dt,
    const double *knots_mu,  int n_mu,
    const double *knots_sig, int n_sig,
    const double *lb_in, const double *ub_in,
    double L, double R,
    LC_SplineType sp_mu, LC_SplineType sp_sig,
    LC_Solver solver,
    int n_starts, int max_iter, double ftol)
{
    _gsl_init();
    auto t0 = std::chrono::steady_clock::now();

    std::vector<double> X0, X1;
    _preprocess(data, n, L, R, X0, X1);
    if (X0.empty()) return nullptr;

    std::vector<double> km_own, ks_own;
    const double *km = knots_mu, *ks = knots_sig;
    if (!km) { km_own = _auto_knots(data, n, L, R, n_mu); km = km_own.data(); }
    if (n_sig == 1) sp_sig = LC_ADDITIVE;
    if (!ks) {
        ks_own = (n_sig == 1) ? std::vector<double>(1, km[0])
                               : _auto_knots(data, n, L, R, n_sig);
        ks = ks_own.data();
    }

    std::vector<double> lb_own, ub_own;
    const double *lb = lb_in, *ub = ub_in;
    if (!lb || !ub) {
        _auto_bounds(n_mu, n_sig, sp_sig, lb_own, ub_own);
        lb = lb_own.data(); ub = ub_own.data();
    }
    int dim = n_mu + n_sig;

    const int NF = 500;
    double data_lo = *std::min_element(X0.begin(), X0.end());
    double data_hi = *std::max_element(X0.begin(), X0.end());
    std::vector<double> fine(NF);
    for (int i = 0; i < NF; ++i)
        fine[i] = data_lo + (data_hi - data_lo) * i / (NF - 1.0);

    Prob prob(X0.data(), X1.data(), X0.size(),
              km, n_mu, ks, n_sig, dt, sp_mu, sp_sig);

    gsl_qrng *qrng = gsl_qrng_alloc(gsl_qrng_sobol, (unsigned)dim);
    std::vector<double> best_par(dim, 0.0);
    double best_f = DBL_INF;
    int    best_rc = 0;
    std::vector<double> sv(dim), par0(dim);

    for (int s = 0; s < n_starts; ++s) {
        gsl_qrng_get(qrng, sv.data());
        for (int j = 0; j < dim; ++j)
            par0[j] = lb[j] + (ub[j] - lb[j]) * sv[j];

        int tries = 0;
        while (!_sigma_ok(par0.data(), n_mu, n_sig, sp_sig, ks,
                           fine.data(), NF) && tries < 100) {
            gsl_qrng_get(qrng, sv.data());
            for (int j = 0; j < dim; ++j)
                par0[j] = lb[j] + (ub[j] - lb[j]) * sv[j];
            ++tries;
        }

        int rc = _run(prob, solver, lb, ub, dim, max_iter, ftol, par0);

        if (prob.best_f < best_f) {
            best_f   = prob.best_f;
            best_par = prob.best_par;
            best_rc  = rc;
        }
        prob.best_f = DBL_INF;
    }
    gsl_qrng_free(qrng);

    double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();

    LC_Result *r = new (std::nothrow) LC_Result{};
    if (!r) return nullptr;
    r->par = new double[dim];
    std::copy(best_par.begin(), best_par.end(), r->par);
    r->n_par      = dim;
    r->negloglik  = best_f;
    r->time_sec   = elapsed;
    r->n_evals    = prob.n_evals;
    r->n_starts   = n_starts;
    r->status     = best_rc;
    return r;
}

extern "C"
void lc_result_free(LC_Result *r) {
    if (!r) return;
    delete[] r->par;
    delete r;
}

extern "C"
void euler_simulate(
    const double *par,
    const double *km, int nm,
    const double *ks, int ns,
    LC_SplineType sm, LC_SplineType ss,
    double dt, double x0, double L, double R,
    int T, long seed, double *out)
{
    _gsl_init();

    auto init_pp = [&](LC_SplineType sp, const double *kn, int nk,
                        const double *vals, PPoly &pp) {
        if (sp == LC_ADDITIVE) {
            pp.sp = LC_ADDITIVE; pp.scalar = vals[0]; return;
        }
        gsl_spline *spl = gsl_spline_alloc(_gsl_itype(sp), nk);
        gsl_spline_init(spl, kn, vals, nk);
        _build_pp(spl, kn, vals, nk, sp, pp);
        gsl_spline_free(spl);
    };

    PPoly pp_m, pp_s;
    init_pp(sm, km, nm, par,      pp_m);
    init_pp(ss, ks, ns, par + nm, pp_s);

    srand((unsigned)(seed < 0 ? (long)time(nullptr) : seed));
    const double sqdt = std::sqrt(dt);
    out[0] = x0;

    for (int t = 1; t < T; ++t) {
        double u1, u2;
        do { u1 = (rand() + 1.0) / (RAND_MAX + 2.0); } while (u1 <= 0.0);
        u2 = (rand() + 1.0) / (RAND_MAX + 2.0);
        double z = std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);

        double xi = out[t-1];
        double mv, sv;
        _eval_pp(pp_m, &xi, 1, &mv);
        _eval_pp(pp_s, &xi, 1, &sv);
        if (!std::isfinite(mv)) mv = 0.0;
        if (!std::isfinite(sv) || sv < 0.0) sv = 0.0;

        double a = xi + mv*dt + sqdt*sv*z;
        if (std::isfinite(L)) a = std::max(a, L);
        if (std::isfinite(R)) a = std::min(a, R);
        out[t] = a;
    }
}

extern "C"
int euler_reconstruct_flat(
    const double *data, int n,       double dt,
    const double *km,   int n_mu,
    const double *ks,   int n_sig,
    const double *lb,   const double *ub,
    double L, double R,
    int sp_mu, int sp_sig, int solver,
    int n_starts, int max_iter, double ftol,
    double *par_out, double *nll_out, double *time_out, int *evals_out)
{
    LC_Result *r = euler_reconstruct(
        data, n, dt,
        km, n_mu, ks, n_sig, lb, ub, L, R,
        (LC_SplineType)sp_mu, (LC_SplineType)sp_sig,
        (LC_Solver)solver,
        n_starts, max_iter, ftol > 0 ? ftol : 1e-9);

    if (!r) return -99;
    std::copy(r->par, r->par + r->n_par, par_out);
    if (nll_out)   *nll_out   = r->negloglik;
    if (time_out)  *time_out  = r->time_sec;
    if (evals_out) *evals_out = r->n_evals;
    int st = r->status;
    lc_result_free(r);
    return st;
}
