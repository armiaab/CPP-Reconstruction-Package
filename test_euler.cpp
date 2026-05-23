/*
 * test_euler.cpp — Verify euler_reconstruct against Tutorial.pdf examples.
 *
 * Compiled as part of the langevin_bench CMake project.
 * Links against langevin_core (which includes gsl_mini + nlopt + Eigen).
 *
 * Tutorial expected outputs (Tutorial.pdf):
 *   Example 3  (CC spline, OUdata1D 20000pts, L=-2, R=1.8, 8k each):
 *     NLL ≈ -17730.17
 *     mu  ≈ [4.89 1.20 1.10 0.25 -0.10 -0.70 -1.32 -6.49]
 *     sig ≈ [0.84 1.02 1.00 0.99  1.00  0.98  1.02  0.80]
 *
 *   Example 5  (CC spline, MayData1D 100000pts, L=min, R=max, 8k each):
 *     NLL ≈ -180071.29
 *     mu  ≈ [0.43  0.11 -0.09  0.002  0.055  0.051 -0.18 -0.58]
 *     sig ≈ [0.38  0.40  0.40  0.41   0.40   0.40   0.40  0.39]
 *
 *   Example 8  (QQ spline, OUdata1D every 100th, L=-2.5,R=2.5, 7k each):
 *     (example_5.m in the repository)
 *     NLL shown in Tutorial pages 20+; reference: params should converge to
 *     approximately mu≈[-x shape] and sigma≈constant near 1.
 */

#include "langevin_core.h"

#include <cstdio>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include <algorithm>
#include <numeric>
#include <cstring>

/* ------------------------------------------------------------------ */
/* Helper: load a binary double file written by extract_data.py        */
/* Format: uint64_t count, then count * double (little-endian)         */
/* ------------------------------------------------------------------ */
static std::vector<double> load_bin(const char *path)
{
    std::ifstream f(path, std::ios::binary);
    if (!f) { fprintf(stderr, "Cannot open %s\n", path); return {}; }
    uint64_t n = 0;
    f.read(reinterpret_cast<char*>(&n), 8);
    std::vector<double> v(n);
    f.read(reinterpret_cast<char*>(v.data()), n * 8);
    return v;
}

/* ------------------------------------------------------------------ */
/* Print result                                                         */
/* ------------------------------------------------------------------ */
static void print_result(const char *tag, LC_Result *r)
{
    if (!r) { printf("[%s]  ERROR: null result\n\n", tag); return; }

    printf("[%s]  NLL = %.4f   (search_agents=%d, evals=%d, t=%.2fs)\n",
           tag, r->negloglik, r->n_starts, r->n_evals, r->time_sec);
    printf("  mu  params: ");
    int half = r->n_par / 2;
    for (int i = 0; i < half; ++i) printf("%.5g ", r->par[i]);
    printf("\n  sig params: ");
    for (int i = half; i < r->n_par; ++i) printf("%.5g ", r->par[i]);
    printf("\n\n");
    lc_result_free(r);
}

/* ------------------------------------------------------------------ */
/* Build equallyspaced knots over [lo, hi], matching MATLAB linspace   */
/* ------------------------------------------------------------------ */
static std::vector<double> linspace(double lo, double hi, int n)
{
    std::vector<double> k(n);
    for (int i = 0; i < n; ++i)
        k[i] = lo + (hi - lo) * i / (n - 1.0);
    return k;
}

/* ------------------------------------------------------------------ */
/* main                                                                 */
/* ------------------------------------------------------------------ */
int main(int argc, char **argv)
{
    /* Resolve data directory (default: same directory as this executable,
       or pass as first argument). */
    std::string dir = ".";
    if (argc >= 2) dir = argv[1];

    auto ou_path  = dir + "/OUdata1D.bin";
    auto may_path = dir + "/MayData1D.bin";

    printf("=== euler_reconstruct verification vs Tutorial.pdf ===\n\n");

    /* ============================================================
     * EXAMPLE 3  (Tutorial p.12)
     * CC spline, OUdata1D first 20000 pts, dt=0.01, L=-2, R=1.8, 8k
     * Expected NLL ≈ -17730.17
     * ============================================================ */
    {
        auto raw = load_bin(ou_path.c_str());
        if (raw.empty()) { fprintf(stderr, "Failed to load OUdata1D\n"); return 1; }

        const int N = 20000;
        std::vector<double> data(raw.begin(), raw.begin() + N);

        /* Set data outside [L,R] to NaN — matches MATLAB: data(data<L|data>R)=nan */
        double L = -2.0, R = 1.8;
        for (double &v : data) if (v < L || v > R) v = std::numeric_limits<double>::quiet_NaN();

        int nk = 8;
        auto km = linspace(L, R, nk);
        auto ks = linspace(L, R, nk);

        /* lb: [-10]*8 ++ [eps]*8,   ub: [10]*16 */
        std::vector<double> lb(2*nk), ub(2*nk, 10.0);
        for (int i = 0;  i < nk; ++i) lb[i]     = -10.0;
        for (int i = nk; i < 2*nk; ++i) lb[i]  =  1e-14;

        printf("--- Example 3: CC spline, OUdata1D 20 000 pts, L=-2 R=1.8, 8 knots ---\n");
        printf("    Expected NLL ≈ -17730.17\n");

        LC_Result *r = euler_reconstruct(
            data.data(), N, /*dt=*/0.01,
            km.data(), nk, ks.data(), nk,
            lb.data(), ub.data(), L, R,
            LC_CUBIC, LC_CUBIC,
            LC_SLSQP, /*n_starts=*/1, /*max_iter=*/10000, /*ftol=*/1e-9);
        print_result("Ex3", r);
    }

    /* ============================================================
     * EXAMPLE 5  (Tutorial p.15)
     * CC spline, MayData1D first 100000 pts, dt=0.01, L=min R=max, 8k
     * Expected NLL ≈ -180071.29
     * ============================================================ */
    {
        auto raw = load_bin(may_path.c_str());
        if (raw.empty()) { fprintf(stderr, "Failed to load MayData1D\n"); return 1; }

        const int N = 100000;
        std::vector<double> data(raw.begin(), raw.begin() + N);

        double L = *std::min_element(data.begin(), data.end());
        double R = *std::max_element(data.begin(), data.end());
        printf("--- Example 5: CC spline, MayData1D 100 000 pts, L=%.4f R=%.4f, 8 knots ---\n", L, R);
        printf("    Expected NLL ≈ -180071.29\n");

        int nk = 8;
        auto km = linspace(L, R, nk);
        auto ks = linspace(L, R, nk);

        std::vector<double> lb(2*nk), ub(2*nk, 10.0);
        for (int i = 0;  i < nk; ++i) lb[i]    = -10.0;
        for (int i = nk; i < 2*nk; ++i) lb[i] =  1e-14;

        LC_Result *r = euler_reconstruct(
            data.data(), N, /*dt=*/0.01,
            km.data(), nk, ks.data(), nk,
            lb.data(), ub.data(), L, R,
            LC_CUBIC, LC_CUBIC,
            LC_SLSQP, /*n_starts=*/1, /*max_iter=*/10000, /*ftol=*/1e-9);
        print_result("Ex5", r);
    }

    /* ============================================================
     * EXAMPLE 8 / example_5.m in repo  (Tutorial p.20)
     * QQ spline, OUdata1D every 100th pt, dt=1, L=-2.5, R=2.5, 7k
     * search_agents=5 (5 multi-starts)
     * ============================================================ */
    {
        auto raw = load_bin(ou_path.c_str());
        if (raw.empty()) { fprintf(stderr, "Failed to load OUdata1D\n"); return 1; }

        /* Every 100th data point */
        std::vector<double> data;
        for (size_t i = 0; i < raw.size(); i += 100) data.push_back(raw[i]);

        double L = -2.5, R = 2.5;
        for (double &v : data) if (v < L || v > R) v = std::numeric_limits<double>::quiet_NaN();

        int nk = 7;
        auto km = linspace(L, R, nk);
        auto ks = linspace(L, R, nk);

        std::vector<double> lb(2*nk), ub(2*nk, 10.0);
        for (int i = 0;  i < nk; ++i) lb[i]    = -10.0;
        for (int i = nk; i < 2*nk; ++i) lb[i] =  1e-14;

        printf("--- Example 8 (example_5.m): QQ spline, OUdata1D every-100th, L=-2.5 R=2.5, 7 knots ---\n");
        printf("    5 multi-starts; sigma should be near 1 (OU noise level)\n");

        LC_Result *r = euler_reconstruct(
            data.data(), (int)data.size(), /*dt=*/1.0,
            km.data(), nk, ks.data(), nk,
            lb.data(), ub.data(), L, R,
            LC_QUADRATIC, LC_QUADRATIC,
            LC_SLSQP, /*n_starts=*/5, /*max_iter=*/10000, /*ftol=*/1e-9);
        print_result("Ex8/ex5", r);
    }

    return 0;
}
