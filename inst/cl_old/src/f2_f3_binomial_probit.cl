// f2_binomial_logit_prep_parallel.cl

#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#pragma OPENCL EXTENSION cl_khr_printf : enable   // for printf

#define MAX_L2 64   // upper bound on l2; tune as needed



// f2 + f3 for Binomial–Probit, single‐pass: prior + data‐term + gradient
__kernel void f2_f3_binomial_probit(
    __global const double* X,      // design matrix,    l1 × l2, column-major
    __global const double* B,      // grid points,      m1 × l2, row-major per grid
    __global const double* mu,     // prior mean,       length = l2
    __global const double* P,      // prior precision,  l2 × l2, row-major
    __global const double* alpha,  // offsets,          length = l1
    __global const double* y,      // successes (0/1),  length = l1
    __global const double* wt,     // trials,           length = l1
    __global double*       qf,     // out: neg‐log‐posterior, length = m1
//    __global double*       xb,     // out: Φ(η),           size = m1 × l1
    __global double*       grad,   // out: ∂(neg‐log‐post)/∂B, size = m1 × l2 (col-major)
    const int l1,                  // # observations
    const int l2,                  // # predictors
    const int m1                   // # grid points
) {
    int j = get_global_id(0);
    if (j >= m1) return;

    // 1) Prior term: tmp[k] = [P × (B_j – mu)]_k
    double tmp[MAX_L2];
    for (int k = 0; k < l2; ++k) {
        double acc = 0.0;
        for (int ℓ = 0; ℓ < l2; ++ℓ) {
            acc += P[k*l2 + ℓ] * (B[j*l2 + ℓ] - mu[ℓ]);
        }
        tmp[k] = acc;
    }

    // 2) Quadratic form: 0.5 * (B_j – mu)' P (B_j – mu)
    double qsum = 0.0;
    for (int k = 0; k < l2; ++k) {
        double d_k = B[j*l2 + k] - mu[k];
        qsum += d_k * tmp[k];
    }
    // initialize neg-log-posterior accumulator
    double res_acc = 0.5 * qsum;

    // 3) Gradient accumulator starts with prior part
    double g_loc[MAX_L2];
    for (int k = 0; k < l2; ++k) {
        g_loc[k] = tmp[k];
    }

    // 4) Data term: loop over each observation
    int base = j * l1;
    for (int i = 0; i < l1; ++i) {
        // linear predictor η_i = α[i] + X[i,·]·B_j
        double eta = alpha[i];
        for (int k = 0; k < l2; ++k) {
            eta += X[k*l1 + i] * B[j*l2 + k];
        }

        // rename for clarity
        double dot = eta;

        // Compute Φ(dot), Φ(–dot), and φ(dot) using high-accuracy routines
        double p1 = pnorm5(dot, 0.0, 1.0, 1, 0);   // Φ(dot)
        double p2 = pnorm5(-dot, 0.0, 1.0, 1, 0);  // Φ(–dot)
        double d  = dnorm4(dot, 0.0, 1.0, 0);      // φ(dot)

        // Store probit probability for this grid point and observation
//        xb[base + i] = p1;

       // <— replace manual log‐lik term with dbinom’s log‐binomial
//        double ll = dbinom(y[i], wt[i], p1, /*give_log=*/1);
        double ll = dbinom_raw(y[i], wt[i], p1, p2, /*give_log=*/1);

        res_acc -= ll;


        // Compute gradient residual for observation i
        double resid = ((y[i] * d / p1) - ((1.0 - y[i]) * d / p2)) * wt[i];

        // Accumulate gradient: g_loc += X[i,·]ᵗ × resid
        for (int k = 0; k < l2; ++k) {
            g_loc[k] -= X[k*l1 + i] * resid;
        }


    }

    // 5) Write back results
    qf[j] = res_acc;
    for (int k = 0; k < l2; ++k) {
        grad[k * m1 + j] = g_loc[k];
    }
}