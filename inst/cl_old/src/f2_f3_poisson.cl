// f2_binomial_logit_prep_parallel.cl

#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#pragma OPENCL EXTENSION cl_khr_printf : enable   // for printf

#define MAX_L2 64   // upper bound on l2; tune as needed


__kernel void f2_f3_poisson(
    __global const double* X,      // design matrix,   l1 × l2, col-major
    __global const double* B,      // grid points,     m1 × l2, row-major per grid
    __global const double* mu,     // prior mean,      length = l2
    __global const double* P,      // prior precision, l2 × l2, row-major
    __global const double* alpha,  // offsets,         length = l1
    __global const double* y,      // counts,          length = l1
    __global const double* wt,     // weights,         length = l1
    __global double*       qf,     // out: neg-log-post, length = m1
//    __global double*       xb,     // out: μ = exp(η),    size = m1 × l1
    __global double*       grad,   // out: ∂(neg-log-post)/∂B, size = l2 × m1 (col-major)
    const int l1,
    const int l2,
    const int m1
) {
    int j = get_global_id(0);
    if (j >= m1) return;

    // 1) Prior term: tmp = P * (B_j - mu)
    double tmp[MAX_L2];
    for (int k = 0; k < l2; ++k) {
        double acc = 0.0;
        for (int ell = 0; ell < l2; ++ell) {
            acc += P[k*l2 + ell] * (B[j*l2 + ell] - mu[ell]);
        }
        tmp[k] = acc;
    }

    // 2) quadratic form qf[j] = 0.5 * (B_j - mu)' * tmp
    double qsum = 0.0;
    for (int k = 0; k < l2; ++k) {
        qsum += (B[j*l2 + k] - mu[k]) * tmp[k];
    }
    qf[j] = 0.5 * qsum;
 

    // 3) Initialize gradient with prior part
    double g_loc[MAX_L2];
    for (int k = 0; k < l2; ++k) {
        g_loc[k] = tmp[k];
    }

    // 4) Data term: Poisson log-link
    int base = j * l1;
    for (int i = 0; i < l1; ++i) {
        // linear predictor η = α[i] + X_i·B_j
        double dot = alpha[i];
        for (int k = 0; k < l2; ++k) {
            dot += X[k*l1 + i] * B[j*l2 + k];
        }

        // fitted mean
        double mui = exp(dot);
//        xb[base + i] = mui;

        // negate log-likelihood contribution

        // Non-integer values requires replacement of dpois with version using lgamma function

//        double logp = dpois(y[i], mui, 1);
          double logp = -mui + y[i] * log(mui) - lgamma(y[i] + 1.0);
        
        qf[j] -= wt[i] * logp;

        // gradient contribution: -(y - μ) * X * wt
        double resid = (y[i] - mui) * wt[i];
        for (int k = 0; k < l2; ++k) {
            g_loc[k] -= X[k*l1 + i] * resid;
        }
    }



    // 5) Write back
  
    for (int k = 0; k < l2; ++k) {
        grad[k*m1 + j] = g_loc[k];
    }
}