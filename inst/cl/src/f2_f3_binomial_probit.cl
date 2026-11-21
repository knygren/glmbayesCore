// @library_deps: nmath
// @calls_nmath: dbinom_raw, pnorm5, dnorm4
// @depends_nmath: dbinom, pnorm, dnorm
// @calls_opencl_builtin: (none)
// @all_depends_nmath_count: 20
// @all_depends_nmath: dpq, refactored, Rmath, nmath, stirlerr_cycle_free, chebyshev, cospi, dnorm, fmax2, gammalims, lgammacor, log1p, pnorm, gamma, lgamma, pgamma_utils, stirlerr_cycle_dependent, bd0, stirlerr, dbinom

#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#pragma OPENCL EXTENSION cl_khr_printf : enable   // for printf

#define MAX_L2 64   // upper bound on l2; tune as needed

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


__kernel void f2_f3_binomial_probit(
    __global const double* X,
    __global const double* B,
    __global const double* mu,
    __global const double* P,
    __global const double* alpha,
    __global const double* y,
    __global const double* wt,
    __global double*       qf,
    __global double*       grad,
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

    for (int i = 0; i < l1; ++i) {
        double eta = alpha[i];
        for (int k = 0; k < l2; ++k) {
            eta += X[k*l1 + i] * B[j*l2 + k];
        }

        double p1 = pnorm5(eta, 0.0, 1.0, 1, 0);
        double p2 = pnorm5(-eta, 0.0, 1.0, 1, 0);
        double d  = dnorm4(eta, 0.0, 1.0, 0);

        res_acc += nll_binomial_glmb_ocl(y[i], wt[i], p1);

        double resid = ((y[i] * d / p1) - ((1.0 - y[i]) * d / p2)) * wt[i];

        for (int k = 0; k < l2; ++k) {
            g_loc[k] -= X[k*l1 + i] * resid;
        }
    }

    qf[j] = res_acc;
    for (int k = 0; k < l2; ++k) {
        grad[k * m1 + j] = g_loc[k];
    }
}
