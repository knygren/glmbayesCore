// f2_binomial_logit_prep_parallel.cl

#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#pragma OPENCL EXTENSION cl_khr_printf : enable   // for printf

#define MAX_L2 64   // upper bound on l2; tune as needed


#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#define MAX_L2 128

// f2 + f3 for Binomial–cloglog, single‐pass: prior + data‐term + gradient
__kernel void f2_f3_binomial_cloglog(
    __global const double* X,      // design matrix, size = l1×l2, col-major
    __global const double* B,      // grid,        size = m1×l2, row-major per‐grid
    __global const double* mu,     // prior mean,  length = l2
    __global const double* P,      // prior prec., size = l2×l2, row-major
    __global const double* alpha,  // offset,      length = l1
    __global const double* y,      // response,    length = l1
    __global const double* wt,     // weights,     length = l1
    __global double*       qf,     // out: quadratics,    length = m1
//    __global double*       xb,     // out: p = cloglog,   size = m1×l1
    __global double*       grad,   // out: dfdB,          size = m1×l2
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
        for (int ℓ = 0; ℓ < l2; ++ℓ) {
            acc += P[k*l2 + ℓ] * (B[j*l2 + ℓ] - mu[ℓ]);
        }
        tmp[k] = acc;
    }

    // 2) quadratic form qf[j] = 0.5 * (B_j - mu)' * tmp
    double qsum = 0.0;
    for (int k = 0; k < l2; ++k) {
        qsum += (B[j*l2 + k] - mu[k]) * tmp[k];
    }
    
    qf[j] = 0.5 * qsum;
    double res_acc = 0.5 * qsum;


    // 3) initialize gradient accumulator with prior part
    double g_loc[MAX_L2];
    for (int k = 0; k < l2; ++k) {
        g_loc[k] = tmp[k];
    }

    // 4) cloglog prep + data‐term for gradient
    int base = j * l1;
    for (int i = 0; i < l1; ++i) {
        // compute linpred = α[i] + x_i·B_j as dot
        double dot = alpha[i];
        for (int k = 0; k < l2; ++k) {
            dot += X[k*l1 + i] * B[j*l2 + k];
        }

//        double p1 = 1.0 - exp(-exp(dot));
        double p1 = -expm1(-exp(dot));            // 1 - exp(-exp(dot)), cancellation-safe
        double p2 = exp(-exp(dot));
        double atemp = exp(dot - exp(dot));
//        xb[base + i] = p1;


//      int trials  = (int)(rint(wt[i]));
//      int success = (int)(rint(y[i] * wt[i]));

//      double p = fmin(1.0, fmax(0.0, p1));  // clamp to [0, 1]

//      double ll = dbinom((double)success, (double)trials, p, /*give_log=*/1);
//      res_acc -= ll;


        // use dbinom for log-likelihood
  //      double ll = dbinom(y[i], wt[i], p1, /*give_log=*/1);
        double ll = dbinom_raw(y[i], wt[i], p1, p2,/*give_log=*/1);

        res_acc -= ll;


        double resid = ((y[i] * atemp / p1) - ((1.0 - y[i]) * atemp / p2)) * wt[i];
        for (int k = 0; k < l2; ++k) {
            g_loc[k] -= X[k*l1 + i] * resid;
        }
    }

    qf[j] = res_acc;

    // 5) write back gradient row for grid‐point j
    for (int k = 0; k < l2; ++k) {
        grad[k * m1 + j] = g_loc[k];  // column-major layout
    }
}