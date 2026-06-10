## Temp debug: locate the Poisson discrepancy between the C++ driver and the
## legacy R loop. Compare (a) single Block 1 call means given fixef_start,
## (b) one full inner step (m_convergence = 1) fixef means.
pkgload::load_all(export_all = FALSE)

set.seed(11)
J <- 8L
n_per <- 20L
grp <- factor(rep(sprintf("g%02d", seq_len(J)), each = n_per))
group_levels <- levels(grp)
w_j <- round(rnorm(J), 2)
z1 <- rnorm(J * n_per)
x_re <- cbind(`(Intercept)` = 1, slope = z1)
re_names <- colnames(x_re)
X_int <- cbind(1, w_j)
rownames(X_int) <- group_levels
colnames(X_int) <- c("(Intercept)", "w")
X_slp <- matrix(1, J, 1L, dimnames = list(NULL, "(Intercept)"))
x_hyper <- list(`(Intercept)` = X_int, slope = X_slp)
gamma_int <- c(1.0, 0.5)
gamma_slp <- 0.8
b_int <- as.numeric(X_int %*% gamma_int) + rnorm(J, sd = 0.4)
b_slp <- as.numeric(X_slp %*% gamma_slp) + rnorm(J, sd = 0.4)
eta <- b_int[as.integer(grp)] + b_slp[as.integer(grp)] * z1

prior_list_block2 <- list(
  `(Intercept)` = list(mu = c(0, 0), Sigma = diag(4, 2L), dispersion = 0.16),
  slope         = list(mu = 0, Sigma = diag(4, 1L), dispersion = 0.16)
)
fixef_start <- list(
  `(Intercept)` = stats::setNames(c(0, 0), colnames(X_int)),
  slope         = stats::setNames(0, colnames(X_slp))
)
set.seed(404); invisible(rnorm(1))
y_pois <- rpois(length(eta), exp(0.3 * eta))
pl1_pois <- list(Sigma = diag(0.25, 2L))

fam <- poisson()
meta <- glmbayesCore:::.two_block_validate_block1_prior(pl1_pois, family = fam)
mu_all0 <- glmbayesCore:::.two_block_mu_all(fixef_start, x_hyper, re_names, group_levels)
pl1 <- glmbayesCore:::.two_block_block1_prior_list(pl1_pois, mu_all0, meta)

## (a) single Block 1 call, fixed mu: R-wrapper path vs identical export call
n_rep <- 400L
acc_r <- matrix(0, J, 2L)
set.seed(1)
for (r in seq_len(n_rep)) {
  out <- block_rNormalGLM(
    n = 1L, y = y_pois, x = x_re, block = grp, prior_list = pl1,
    family = fam, use_parallel = FALSE
  )
  acc_r <- acc_r + unname(out$coefficients)
}
famfunc <- glmbayesCore::glmbfamfunc(fam)
acc_c <- matrix(0, J, 2L)
set.seed(2)
for (r in seq_len(n_rep)) {
  out <- glmbayesCore:::.block_rNormalGLM_cpp(
    n = 1L, y = y_pois, x = x_re, block = grp,
    prior_list = pl1, prior_lists = NULL,
    offset = rep(0, length(y_pois)), wt = rep(1, length(y_pois)),
    f2 = famfunc$f2, f3 = famfunc$f3,
    family = "poisson", link = "log",
    Gridtype = 2L, n_envopt = 1L,
    use_parallel = FALSE, use_opencl = FALSE, verbose = FALSE
  )
  acc_c <- acc_c + unname(out$coefficients)
}
cat("(a) Block1-only mean diff:", max(abs(acc_r / n_rep - acc_c / n_rep)), "\n")

## (b) one inner step: full driver (m=1, n=n_rep) vs legacy loop (m=1)
fit_cpp <- two_block_rNormal_reg(
  n = n_rep, y = y_pois, x = x_re, block = grp, x_hyper = x_hyper,
  prior_list_block1 = pl1_pois, prior_list_block2 = prior_list_block2,
  fixef_start = fixef_start, m_convergence = 1L, family = fam,
  use_parallel = FALSE, seed = 3L, progbar = FALSE
)

set.seed(4)
fx_sum <- list(`(Intercept)` = c(0, 0), slope = 0)
b_sum <- matrix(0, J, 2L)
for (r in seq_len(n_rep)) {
  fixef <- fixef_start
  mu_all <- glmbayesCore:::.two_block_mu_all(fixef, x_hyper, re_names, group_levels)
  pl1r <- glmbayesCore:::.two_block_block1_prior_list(pl1_pois, mu_all, meta)
  block_i <- block_rNormalGLM(
    n = 1L, y = y_pois, x = x_re, block = grp, prior_list = pl1r,
    family = fam, use_parallel = FALSE
  )
  b_i <- block_i$coefficients
  colnames(b_i) <- re_names
  fixef_draw <- multi_rNormal_reg(
    n = 1L, y = b_i, x = x_hyper, prior_list = prior_list_block2,
    family = gaussian(), progbar = FALSE
  )
  for (k in re_names) fx_sum[[k]] <- fx_sum[[k]] + fixef_draw[[k]]$coefficients[1L, ]
  b_sum <- b_sum + unname(b_i)
}
for (k in re_names) {
  cat("(b) fixef[", k, "] cpp:", colMeans(fit_cpp$fixef_draws[[k]]),
      " ref:", fx_sum[[k]] / n_rep, "\n")
}
b_cpp <- sapply(seq_along(re_names), function(jj) {
  tapply(fit_cpp$coefficients[[re_names[jj]]], fit_cpp$coefficients$grp, mean)
})
cat("(b) b mean diff:", max(abs(unname(b_cpp) - b_sum / n_rep)), "\n")
