## two_block_rate() — Remark 8 contraction rate for a balanced
## random-intercept model, where lambda* has a closed form.
##
## For J groups of n_j observations each, p_re = q = 1, and a diffuse Block 2
## prior,  lambda* = sigma^2 / (sigma^2 + n_j * tau^2).

J      <- 12L
n_j    <- 8L
sigma2 <- 1.0
tau2   <- 0.5

group   <- factor(rep(seq_len(J), each = n_j))
Z       <- matrix(1, J * n_j, 1L, dimnames = list(NULL, "(Intercept)"))
x_hyper <- list("(Intercept)" = matrix(
  1, J, 1L, dimnames = list(levels(group), "(Intercept)")
))

rate <- two_block_rate(
  x                 = Z,
  group             = group,
  x_hyper           = x_hyper,
  prior_list_block1 = list(Sigma = matrix(tau2), dispersion = sigma2),
  prior_list_block2 = list("(Intercept)" = list(mu = 0, Sigma = matrix(100)))
)

c(computed = rate$lambda_star,
  closed_form = sigma2 / (sigma2 + n_j * tau2))

## print() adds the sweeps implied by each bound: the geometric (lambda*)^m
## proxy, the Corollary 1 envelope, and the exact Theorem 3 terms the sampler
## actually uses.
rate

## Larger groups carry more likelihood precision and mix faster; the rate does
## not depend on the number of groups.
lambda_for <- function(n_j) {
  g <- factor(rep(seq_len(J), each = n_j))
  two_block_rate(
    x       = matrix(1, J * n_j, 1L, dimnames = list(NULL, "(Intercept)")),
    group   = g,
    x_hyper = list("(Intercept)" = matrix(
      1, J, 1L, dimnames = list(levels(g), "(Intercept)")
    )),
    prior_list_block1 = list(Sigma = matrix(tau2), dispersion = sigma2),
    prior_list_block2 = list("(Intercept)" = list(mu = 0, Sigma = matrix(100)))
  )$lambda_star
}
round(vapply(c(2L, 8L, 32L), lambda_for, numeric(1)), 4)
