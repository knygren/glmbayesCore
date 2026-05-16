// f2_f3_gaussian.cl

#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#pragma OPENCL EXTENSION cl_khr_printf : enable   // for printf

#define MAX_L2 64   // upper bound on l2; tune as needed

__kernel void f2_f3_gaussian(
    __global const double* X,    // l1 x l2, column-major: X[k*l1 + i]
    __global const double* B,    // m1 x l2, row-major:    B[j*l2 + k]
    __global const double* mu,   // l2
    __global const double* P,    // l2 x l2, row-major
    __global const double* alpha,// l1
    __global const double* y,    // l1
    __global const double* wt,   // l1 (precision)
    __global double*       qf,   // m1
    __global double*       grad, // m1 x l2, column-major: grad[k*m1 + j]
    const int l1,
    const int l2,
    const int m1
) {
    int j = get_global_id(0);
    if (j >= m1) return;

    // tmp = P * (B_j - mu)
    double tmp[MAX_L2];
    for (int k = 0; k < l2; ++k) {
        double acc = 0.0;
        for (int ell = 0; ell < l2; ++ell) {
            acc += P[k*l2 + ell] * (B[j*l2 + ell] - mu[ell]);
        }
        tmp[k] = acc;
    }

    // objective: 0.5*(B_j - mu)'*tmp
    double qsum = 0.0;
    for (int k = 0; k < l2; ++k) {
        double dk = B[j*l2 + k] - mu[k];
        qsum += dk * tmp[k];
    }
    double res_acc = 0.5 * qsum;

    // gradient starts with prior term
    double g_loc[MAX_L2];
    for (int k = 0; k < l2; ++k) g_loc[k] = tmp[k];

    // data term
    for (int i = 0; i < l1; ++i) {
        // dot = alpha[i] + X[i,*] %* B_j
        double dot = alpha[i];
        for (int k = 0; k < l2; ++k) {
            dot += X[k*l1 + i] * B[j*l2 + k];
        }

        // objective: -log dnorm(y | mean=dot, sd=1/sqrt(wt))
        double wi   = wt[i];
        double sd_i = 1.0 / sqrt(wi);
        double ll   = dnorm4(y[i], dot, sd_i, 1);
        res_acc    += -ll;

        // gradient contribution: Xᵀ * (wt * (dot - y))
        double resid = wi * (dot - y[i]);
        for (int k = 0; k < l2; ++k) {
            g_loc[k] += X[k*l1 + i] * resid;
        }
    }

    // write back
    qf[j] = res_acc;
    for (int k = 0; k < l2; ++k) {
        grad[k*m1 + j] = g_loc[k];
    }
}