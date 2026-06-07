## Conditionally independent block Gaussian draw (3 school groups)

set.seed(42)

## Simulate a simple two-level data set: 3 schools, ~10 students each
n_schools <- 3L
n_per     <- 10L
school    <- rep(seq_len(n_schools), each = n_per)
x         <- cbind(1, rnorm(n_schools * n_per))  # intercept + covariate
b_true    <- matrix(c(5, 0.5, 3, -0.2, 7, 0.3), nrow = n_schools, byrow = TRUE)
sigma2    <- 1.5    # residual variance

y <- rowSums(x * b_true[school, ]) + rnorm(nrow(x), sd = sqrt(sigma2))

## Flat prior (large Sigma) shared across all schools; dispersion = sigma2
l1         <- ncol(x)
prior_list <- list(
  mu         = rep(0, l1),
  Sigma      = diag(100, l1),
  dispersion = sigma2,
  ddef       = FALSE
)

out <- block_rNormalReg(
  n          = 1L,
  y          = y,
  x          = x,
  block      = school,
  prior_list = prior_list
)

out$coefficients   ## k x l1 matrix: one row of b_j draws per school
out$coef.mode
