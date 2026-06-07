## One block_rNormalReg_update step: draw b_j given current fixed effects

set.seed(7)

## Simulate data: 3 schools, 10 students each
n_schools <- 3L
n_per     <- 10L
school    <- rep(seq_len(n_schools), each = n_per)
Z         <- cbind(1, rnorm(n_schools * n_per))  # within-school design
b_true    <- matrix(c(5, 0.5, 3, -0.2, 7, 0.3), nrow = n_schools, byrow = TRUE)
sigma2    <- 1.5

y <- rowSums(Z * b_true[school, ]) + rnorm(nrow(Z), sd = sqrt(sigma2))

## Precision matrix for the random-effect prior (shared across schools)
l1  <- ncol(Z)
P_b <- diag(0.01, l1)   # vague prior

## Current "fixed-effects" prior mean (one vector recycled to all blocks)
mu_current <- rep(0, l1)

out <- block_rNormalReg_update(
  mu_all     = mu_current,   ## recycled to all k blocks (prior_list path)
  P          = P_b,
  dispersion = sigma2,
  y          = y,
  x          = Z,
  block      = school
)

out$b_draws         ## k x l1 matrix of updated random effects
out$coefficients    ## same content (all columns)
