// @library_deps: nmath
// @calls_nmath: dbinom_raw
// @depends_nmath: dbinom
// @calls_opencl_builtin: (none)
// @all_depends_nmath_count: 18
// @all_depends_nmath: dpq, refactored, Rmath, nmath, stirlerr_cycle_free, chebyshev, cospi, fmax2, gammalims, lgammacor, log1p, gamma, lgamma, pgamma_utils, stirlerr_cycle_dependent, bd0, stirlerr, dbinom

#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#pragma OPENCL EXTENSION cl_khr_printf : enable   // for printf

#define MAX_L2 64   // upper bound on l2; tune as needed

// Matches CPU dbinom_glmb + glm convention: y is success proportion, wt is trial count.
static inline double nll_binomial_glmb_ocl(double y_prop, double wt, double mean_p_raw) {
    // int trials  = (int)round(wt);
    // int success = (int)round(y_prop * wt);
    double trials  = wt;
    double success = y_prop * wt;
    double p = fmin(1.0, fmax(0.0, mean_p_raw));
    double q = 1.0 - p;
    // double logpmf = dbinom_raw((double)success, (double)trials, p, q, 1);
    double logpmf = dbinom_raw(success, trials, p, q, 1);
    return -logpmf;
}


__kernel void f2_f3_binomial_logit(
    __global const double* X,      // design matrix, size = l1×l2, col-major
    __global const double* B,      // grid,        size = m1×l2, row-major per‐grid
    __global const double* mu,     // prior mean,  length = l2
    __global const double* P,      // prior prec., size = l2×l2, row-major
    __global const double* alpha,  // offset,      length = l1
    __global const double* y,      // proportion of successes (glm / dbinom_glmb)
    __global const double* wt,     // trial counts / weights
    __global double*       qf,     // out: quadratics,    length = m1
    __global double*       grad,   // out: dfdB,        size = m1×l2
    const int l1,
    const int l2,
    const int m1
) {
    int j = get_global_id(0);
    if (j >= m1) return;

    double tmp[MAX_L2];
    for (int k = 0; k < l2; ++k) {
        double acc = 0.0;
        for (int ℓ = 0; ℓ < l2; ++ℓ) {
            acc += P[k*l2 + ℓ] * (B[j*l2 + ℓ] - mu[ℓ]);
        }
        tmp[k] = acc;
    }

    double qsum = 0.0;
    for (int k = 0; k < l2; ++k) {
        double d_k = B[j*l2 + k] - mu[k];
        qsum += d_k * tmp[k];
    }
    double res_acc = 0.5 * qsum;

    double g_loc[MAX_L2];
    for (int k = 0; k < l2; ++k) {
        g_loc[k] = tmp[k];
    }

    double p, q, e;

    for (int i = 0; i < l1; ++i) {
        double dot = -alpha[i];
        for (int k = 0; k < l2; ++k) {
            dot -= X[k*l1 + i] * B[j*l2 + k];
        }

        if (dot <= 0) {
            e = exp(dot);
            p = 1.0 / (1.0 + e);
            q = e / (1.0 + e);
        } else {
            e = exp(-dot);
            p = e / (1.0 + e);
            q = 1.0 / (1.0 + e);
        }

        res_acc += nll_binomial_glmb_ocl(y[i], wt[i], p);

        double resid = (p - y[i]) * wt[i];
        for (int k = 0; k < l2; ++k) {
            g_loc[k] += X[k*l1 + i] * resid;
        }
    }

    qf[j] = res_acc;

    for (int k = 0; k < l2; ++k) {
        grad[k * m1 + j] = g_loc[k];
    }
}
