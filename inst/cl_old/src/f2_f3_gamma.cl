// f2_binomial_logit_prep_parallel.cl

#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#pragma OPENCL EXTENSION cl_khr_printf : enable   // for printf

#define MAX_L2 64   // upper bound on l2; tune as needed


__kernel void f2_f3_gamma(
    __global const double* X,      // design matrix:       l1 × l2, column-major
    __global const double* B,      // grid of β’s:         m1 × l2, row-major
    __global const double* mu,     // prior means:         length = l2
    __global const double* P,      // prior precision:     l2 × l2, row-major
    __global const double* alpha,  // offsets α:           length = l1
    __global const double* y,      // responses y:         length = l1
    __global const double* wt,     // weights/shape:       length = l1
    __global double*       qf,     // out: log-posterior:  length = m1
//    __global double*       xb,     // out: xb = μ/ wt:     size = m1 × l1
    __global double*       grad,   // out: ∂ℓ/∂β:           size = m1 × l2
    const int l1,
    const int l2,
    const int m1
) {
    int j = get_global_id(0);
    if (j >= m1) return;

    // 1) Prior: tmp = P * (B_j - mu)
    double tmp[MAX_L2];
    for (int k = 0; k < l2; ++k) {
        double acc = 0.0;
        for (int ell = 0; ell < l2; ++ell) {
            acc += P[k*l2 + ell] * (B[j*l2 + ell] - mu[ell]);
        }
        tmp[k] = acc;
    }

    // 2) Prior quadratic form: 0.5 * (B_j - mu)' * tmp
    double qsum = 0.0;
    for (int k = 0; k < l2; ++k) {
        double d = B[j*l2 + k] - mu[k];
        qsum += d * tmp[k];
    }
    double res_acc = 0.5 * qsum;

    // 3) Initialize gradient with prior part
    double g_loc[MAX_L2];
    for (int k = 0; k < l2; ++k) {
        g_loc[k] = tmp[k];
    }

    // 4) Data term: μ = exp(η), scaled xb, log-lik via dgamma, gradient
    int base = j * l1;
    for (int i = 0; i < l1; ++i) {
        // linear predictor η_i = α[i] + X_i · B_j
        double eta = alpha[i];
        for (int k = 0; k < l2; ++k) {
            eta += X[k*l1 + i] * B[j*l2 + k];
        }

        // μ_i and scaled xb_i
        double mui    = exp(eta);
        double xb_i   = mui / wt[i];
//        xb[base + i]  = xb_i;

        // exact log-likelihood term from your C: ll = dgamma(y, wt, xb, TRUE)
        double ll     = dgamma(y[i], wt[i], xb_i, 1);
        res_acc      += -ll;

        // gradient contribution: (1 - y/μ) * wt
        double resid = (1.0 - (y[i] / mui)) * wt[i];
        for (int k = 0; k < l2; ++k) {
            g_loc[k] += X[k*l1 + i] * resid;
        }
    }

    // 5) Write back total log-posterior
    qf[j] = res_acc;

    // 6) Write back gradient ∂ℓ/∂β for grid-point j
    for (int k = 0; k < l2; ++k) {
        grad[k * m1 + j] = g_loc[k];
    }
}