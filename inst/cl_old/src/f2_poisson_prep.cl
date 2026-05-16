// f2_binomial_logit_prep_parallel.cl

#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#pragma OPENCL EXTENSION cl_khr_printf : enable   // for printf

#define MAX_L2 64   // upper bound on l2; tune as needed




__kernel void f2_poisson_prep_grad(
    __global const double* X,      // design matrix, size = l1×l2, column-major
    __global const double* B,      // grid,        size = m1×l2, row-major
    __global const double* mu,     // prior mean,  length = l2
    __global const double* P,      // prior prec., size = l2×l2, row-major
    __global const double* alpha,  // offset,      length = l1
    __global const double* y,      // response,    length = l1
    __global const double* wt,     // weights,     length = l1
    __global double*       qf,     // out: quadratics,    length = m1
    __global double*       xb,     // out: μ = exp(η),     size = m1×l1
    __global double*       grad,   // out: dfdB,           size = m1×l2
    const int l1,
    const int l2,
    const int m1
) {
    int j = get_global_id(0);
    if (j >= m1) return;

    // 1) compute tmp = P * (B_j - mu)
    double tmp[MAX_L2];
    for (int k = 0; k < l2; ++k) {
        double acc = 0.0;
        for (int ell = 0; ell < l2; ++ell) {
            acc += P[k*l2 + ell] * (B[j*l2 + ell] - mu[ell]); // P is row-major
        }
        tmp[k] = acc;
    }

    // 2) quadratic form qf[j] = 0.5 * (B_j - mu)' * tmp
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

    // 4) Poisson log link data-term
    int base = j * l1; // for xb[j, *]
    for (int i = 0; i < l1; ++i) {
        // η_i = α[i] + x_i · B_j
        double dot = alpha[i];
        for (int k = 0; k < l2; ++k) {
            dot += X[k*l1 + i] * B[j*l2 + k];
        }

        double mui = exp(dot);           // μ_i
        xb[base + i] = mui;              // store μ_i

        double resid = (y[i] - mui) * wt[i];  // residual
        for (int k = 0; k < l2; ++k) {
            g_loc[k] -= X[k*l1 + i] * resid;  // accumulate gradient
        }
    }

    // 5) write back gradient row for grid point j
    for (int k = 0; k < l2; ++k) {
        grad[k * m1 + j] = g_loc[k];  // column-major layout
    }
}