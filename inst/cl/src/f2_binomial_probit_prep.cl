// f2_binomial_logit_prep_parallel.cl

#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#pragma OPENCL EXTENSION cl_khr_printf : enable   // for printf

#define MAX_L2 64   // upper bound on l2; tune as needed



__kernel void f2_binomial_probit_prep_grad(
    __global const double* X,      // Design matrix: l1 × l2, column-major
    __global const double* B,      // Grid points: m1 × l2, row-major per grid
    __global const double* mu,     // Prior mean: length l2
    __global const double* P,      // Prior precision matrix: l2 × l2, row-major
    __global const double* alpha,  // Offset vector: length l1
    __global const double* y,      // Response vector: length l1
    __global const double* wt,     // Weights: length l1
    __global double*       qf,     // Output: quadratic form for each grid point
    __global double*       xb,     // Output: probit probabilities Φ(η), size m1 × l1
    __global double*       grad,   // Output: gradient matrix, size m1 × l2
    const int l1,                  // Number of observations
    const int l2,                  // Number of predictors
    const int m1                   // Number of grid points
) {
    int j = get_global_id(0);      // Grid point index
    if (j >= m1) return;

    // ───────────────────────────────────────────────────────
    // Step 1: Compute tmp = P × (B_j – mu)
    // This is the prior gradient component for grid point j
    // ───────────────────────────────────────────────────────
    double tmp[MAX_L2];
    for (int k = 0; k < l2; ++k) {
        double acc = 0.0;
        for (int ℓ = 0; ℓ < l2; ++ℓ) {
            acc += P[k*l2 + ℓ] * (B[j*l2 + ℓ] - mu[ℓ]);
        }
        tmp[k] = acc;
    }

    // ───────────────────────────────────────────────────────
    // Step 2: Compute quadratic form qf[j] = 0.5 × (B_j – mu)' × P × (B_j – mu)
    // Used in envelope construction
    // ───────────────────────────────────────────────────────
    double qsum = 0.0;
    for (int k = 0; k < l2; ++k) {
        double d_k = B[j*l2 + k] - mu[k];
        qsum += d_k * tmp[k];
    }
    qf[j] = 0.5 * qsum;

    // ───────────────────────────────────────────────────────
    // Step 3: Initialize gradient accumulator with prior part
    // ───────────────────────────────────────────────────────
    double g_loc[MAX_L2];
    for (int k = 0; k < l2; ++k) {
        g_loc[k] = tmp[k];
    }

    // ───────────────────────────────────────────────────────
    // Step 4: Loop over observations to compute probit probabilities
    //         and accumulate gradient contributions
    // ───────────────────────────────────────────────────────
    int base = j * l1;  // Offset for xb[j,·]
    for (int i = 0; i < l1; ++i) {
        // Compute linear predictor η_i = α_i + x_i · B_j
        double dot = alpha[i];
        for (int k = 0; k < l2; ++k) {
            dot += X[k*l1 + i] * B[j*l2 + k];
        }

        // Compute Φ(dot), Φ(–dot), and φ(dot) using high-accuracy routines
        double p1 = pnorm5(dot, 0.0, 1.0, 1, 0);   // Φ(dot)
        double p2 = pnorm5(-dot, 0.0, 1.0, 1, 0);  // Φ(–dot)
        double d  = dnorm4(dot, 0.0, 1.0, 0);      // φ(dot)

        // Store probit probability for this grid point and observation
        xb[base + i] = p1;

        // Compute gradient residual for observation i
        double resid = ((y[i] * d / p1) - ((1.0 - y[i]) * d / p2)) * wt[i];

        // Accumulate gradient: g_loc += X[i,·]ᵗ × resid
        for (int k = 0; k < l2; ++k) {
            g_loc[k] -= X[k*l1 + i] * resid;
        }
    }

    // ───────────────────────────────────────────────────────
    // Step 5: Write back gradient vector for grid point j
    //         Stored column-major: grad[k × m1 + j]
    // ───────────────────────────────────────────────────────
    for (int k = 0; k < l2; ++k) {
        grad[k * m1 + j] = g_loc[k];
    }
}