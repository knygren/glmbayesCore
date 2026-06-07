## One Gibbs block update: scalar intercept per observation (Poisson)

set.seed(42)

n <- 5L
y <- rpois(n, 2)
x <- matrix(1, n, 1)
mu <- rep(0, n)

upd <- block_rNormalGLM_update(
  mu_all = mu,
  sigma_theta_sq = 1,
  y = y,
  x = x,
  block = seq_len(n),
  family = poisson(),
  use_parallel = FALSE
)

upd$theta
upd$coefficients
