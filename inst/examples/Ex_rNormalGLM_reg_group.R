## rNormalGLM_reg_group() — independent Poisson GLMs by group (toy)

set.seed(1)
group <- factor(rep(c("A", "B"), each = 15L))
x <- cbind(1, rnorm(length(group)))
colnames(x) <- c("(Intercept)", "x1")
eta <- x %*% c(-0.5, 0.3) + ifelse(group == "A", 0, 0.4)
y <- rpois(length(group), exp(eta))

prior_list <- list(
  mu = rep(0, ncol(x)),
  P  = diag(10, ncol(x)),
  dispersion = 1
)

out <- rNormalGLM_reg_group(
  n = 1L,
  y = y,
  x = x,
  group = group,
  prior_list = prior_list,
  family = poisson(),
  use_parallel = FALSE
)

out$coefficients
out$coef.mode
