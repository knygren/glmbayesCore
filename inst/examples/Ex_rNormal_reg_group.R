## rNormal_reg_group() — independent Gaussian regressions by group (toy)

set.seed(42)
n_groups <- 3L
n_per    <- 10L
school   <- rep(seq_len(n_groups), each = n_per)
x        <- cbind(1, rnorm(n_groups * n_per))
colnames(x) <- c("(Intercept)", "X1")
b_true   <- matrix(c(5, 0.5, 3, -0.2, 7, 0.3), nrow = n_groups, byrow = TRUE)
sigma2   <- 1.5
y        <- rowSums(x * b_true[school, ]) + rnorm(nrow(x), sd = sqrt(sigma2))
l1       <- ncol(x)

prior_list <- list(
  mu         = rep(0, l1),
  Sigma      = diag(100, l1),
  dispersion = sigma2,
  ddef       = FALSE
)

out <- rNormal_reg_group(
  n = 1L, y = y, x = x, group = school, prior_list = prior_list
)

out$coefficients
out$coef.mode
