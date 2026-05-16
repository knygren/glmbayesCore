// f2_binomial_logit_prep_parallel.cl

#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#pragma OPENCL EXTENSION cl_khr_printf : enable   // for printf

#define MAX_L2 64   // upper bound on l2; tune as needed



__kernel void f2_f3_binomial_logit(
    __global const double* X,      // design matrix, size = l1×l2, col-major
    __global const double* B,      // grid,        size = m1×l2, row-major per‐grid
    __global const double* mu,     // prior mean,  length = l2
    __global const double* P,      // prior prec., size = l2×l2, row-major
    __global const double* alpha,  // offset,      length = l1
    __global const double* y,      // response,    length = l1
    __global const double* wt,     // weights,     length = l1
    __global double*       qf,     // out: quadratics,    length = m1
//    __global double*       xb,     // out: p = logistic, size = m1×l1
    __global double*       grad,   // out: dfdB,        size = m1×l2
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
        double d_k = B[j*l2 + k] - mu[k];
        qsum += d_k * tmp[k];
    }
//    qf[j] = 0.5 * qsum;
    // res_acc will hold 0.5 * (B_j−mu)'P(B_j−mu)  plus  sum of –log‐likelihoods
    double res_acc = 0.5 * qsum;


    // 3) initialize gradient accumulator with prior part
    double g_loc[MAX_L2];
    for (int k = 0; k < l2; ++k) {
        g_loc[k] = tmp[k];
    }

    double p,q, e;

    // 4) logistic prep + dbinom call + data‐term for gradient
    int base = j * l1;
    for (int i = 0; i < l1; ++i) {
        // compute linpred = α[i] + X[i,·]·B_j as dot
        double dot = -alpha[i];
        for (int k = 0; k < l2; ++k) {
            dot -= X[k*l1 + i] * B[j*l2 + k];
        }

        if(dot<=0){
        e=exp(dot);
        p = 1.0 / (1.0 + e);
        q = e / (1.0 + e);
        }
        else{
        e=exp(-dot);
        p = e / (1.0 + e);
        q = 1.0 / (1.0 + e);
        }


//        double p = 1.0 / (1.0 + exp(dot));

//        xb[base + i] = p;

        // call dbinom on log‐scale and negate to match yy = -dbinom_glmb(...)
        // give_log=1 returns log-density
        // double ll = dbinom(y[i], wt[i], p, 1);
        double ll = dbinom_raw(y[i], wt[i], p,q, 1);
//        double yy = -ll;
        res_acc += -ll;

        // for now we ignore yy; later hook it to an output buffer

        // accumulate gradient: X[i,·]^T * ((p - y[i]) * wt[i])
        
        double resid;
        if(p<0.5){
        
        resid = p * wt[i] - y[i]* wt[i] ;
        }
        else{
        
        resid = (1- y[i])* wt[i] -q * wt[i] ;
        
        }
        
        for (int k = 0; k < l2; ++k) {
            g_loc[k] += X[k*l1 + i] * resid;
        }
    }

    // 5) write back final objective = prior + data‐term
    qf[j] = res_acc;

    // 6) write back gradient row for grid‐point j
    for (int k = 0; k < l2; ++k) {
        grad[k * m1 + j] = g_loc[k];  // column-major layout
    }
}