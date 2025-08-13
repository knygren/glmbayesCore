// f2_binomial_logit_prep_parallel.cl

#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#pragma OPENCL EXTENSION cl_khr_printf : enable   // for printf

#define MAX_L2 64   // upper bound on l2; tune as needed

#pragma OPENCL EXTENSION cl_khr_fp64   : enable   // for double



__kernel void f2_gamma_prep_grad(
    __global const double* X,      // design matrix, size = l1×l2, column-major
    __global const double* B,      // grid,        size = m1×l2, row-major
    __global const double* mu,     // prior mean,  length = l2
    __global const double* P,      // prior prec., size = l2×l2, row-major
    __global const double* alpha,  // offset,      length = l1
    __global const double* y,      // response,    length = l1
    __global const double* wt,     // weights,     length = l1
    __global double*       qf,     // out: prior quadratic terms, length = m1
    __global double*       xb,     // out: μ = exp(η),            size = m1×l1
    __global double*       grad,   // out: gradient dfdB,         size = m1×l2 (column-major by coef)
    const int l1,
    const int l2,
    const int m1
) {
    int j = get_global_id(0);
    if (j >= m1) return;

    // 1) tmp = P * (B_j - mu)   (prior precision times deviation)
    double tmp[MAX_L2];
    for (int k = 0; k < l2; ++k) {
        double acc = 0.0;
        for (int ell = 0; ell < l2; ++ell) {
            // P is row-major: row k, col ell
            acc += P[k*l2 + ell] * (B[j*l2 + ell] - mu[ell]);
        }
        tmp[k] = acc;
    }

    // 2) prior quadratic form: qf[j] = 0.5 * (B_j - mu)' * tmp
    double qsum = 0.0;
    for (int k = 0; k < l2; ++k) {
        qsum += (B[j*l2 + k] - mu[k]) * tmp[k];
    }
    qf[j] = 0.5 * qsum;

    // 3) initialize gradient accumulator with prior part
    double g_loc[MAX_L2];
    for (int k = 0; k < l2; ++k) {
        g_loc[k] = tmp[k];
    }

    // 4) Gamma log-link data-term
    //    η_i = α[i] + x_i · B_j
    //    μ_i = exp(η_i)
    //    gradient contribution per obs: X^T * ((1 - y/μ) * wt)
    int base = j * l1; // base index for xb[j, *]
    for (int i = 0; i < l1; ++i) {
        double eta = alpha[i];
        // X is column-major: column k, row i at X[k*l1 + i]
        for (int k = 0; k < l2; ++k) {
            eta += X[k*l1 + i] * B[j*l2 + k];
        }

        // μ_i
        double mui = exp(eta);
        xb[base + i] = mui;  // store μ_i (unscaled)

        // residual for gamma with log link: (1 - y/μ) * wt
        // Note: mui > 0 by construction. If needed, add small epsilon for safety.
        double resid = (1.0 - (y[i] / mui)) * wt[i];

        // accumulate gradient: g += X^T * resid
        for (int k = 0; k < l2; ++k) {
            g_loc[k] += X[k*l1 + i] * resid;
        }
    }

    // 5) write back gradient for grid point j
    // grad is stored with coefficient-major leading dimension: grad[k, j] at grad[k*m1 + j]
    for (int k = 0; k < l2; ++k) {
        grad[k * m1 + j] = g_loc[k];
    }
}