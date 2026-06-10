## Tests for two_block_mode_weights(): per-observation likelihood precisions
## (IRLS/Fisher weights) evaluated at a supplied random-effects value.
##
## Checks:
##   1. gaussian: weights == 1/dispersion and two_block_rate(weights = w)
##      reproduces the exact (no-weights) spectrum
##   2. poisson-log / binomial-logit: info_total equals glmbfamfunc()$f7
##      evaluated at the same point; closed-form weights match
##   3. B_lik matches brute force per group; info_total = sum(B_lik)
##   4. probit: generic formula phi(eta)^2/(p(1-p)) (f7's logistic weights are
##      wrong for probit, so compare against the closed form directly)
##   5. plumbing: weights feed two_block_rate() and give a spectrum in [0, 1)
##   6. input validation
## Run: Rscript data-raw/test_two_block_mode_weights.R

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload.", call. = FALSE)
}
pkgload::load_all(export_all = FALSE)

## ---------------------------------------------------------------------------
## Fixture: same toy design as test_two_block_rate.R
## ---------------------------------------------------------------------------
set.seed(11)
J <- 8L
n_per <- 20L
grp <- factor(rep(sprintf("g%02d", seq_len(J)), each = n_per))
group_levels <- levels(grp)
w_j <- round(rnorm(J), 2)

z1 <- rnorm(J * n_per)
x_re <- cbind(`(Intercept)` = 1, slope = z1)
re_names <- colnames(x_re)
p_re <- length(re_names)

X_int <- cbind(1, w_j)
rownames(X_int) <- group_levels
colnames(X_int) <- c("(Intercept)", "w")
X_slp <- matrix(1, J, 1L, dimnames = list(NULL, "(Intercept)"))
x_hyper <- list(`(Intercept)` = X_int, slope = X_slp)

prior_list_block2 <- list(
  `(Intercept)` = list(mu = c(0, 0), Sigma = diag(4, 2L), dispersion = 0.16),
  slope         = list(mu = 0, Sigma = diag(4, 1L), dispersion = 0.16)
)

## A b_mode matrix to evaluate curvature at (any point is valid; weights are
## a function of the supplied point, not of any fit)
b_mode <- cbind(
  `(Intercept)` = 0.6 + 0.3 * w_j,
  slope         = rep(0.4, J)
)
rownames(b_mode) <- group_levels

disp <- 0.25
prior_b <- list(P = diag(c(2.5, 4.0)), dispersion = disp)

## ---------------------------------------------------------------------------
## 1. Gaussian: weights = 1/dispersion; rate(weights) == exact rate
## ---------------------------------------------------------------------------
mw_g <- two_block_mode_weights(
  x = x_re, block = grp, b_mode = b_mode,
  family = gaussian(), dispersion = disp, group_levels = group_levels
)
stopifnot(inherits(mw_g, "two_block_mode_weights"))
stopifnot(all(abs(mw_g$weights - 1 / disp) < 1e-15))
stopifnot(all(abs(mw_g$eta - (x_re %*% t(b_mode))[cbind(seq_along(grp), as.integer(grp))]) < 1e-12))

rate_exact <- two_block_rate(
  x = x_re, block = grp, x_hyper = x_hyper,
  prior_list_block1 = prior_b, prior_list_block2 = prior_list_block2,
  family = gaussian(), group_levels = group_levels
)
rate_w <- two_block_rate(
  x = x_re, block = grp, x_hyper = x_hyper,
  prior_list_block1 = prior_b, prior_list_block2 = prior_list_block2,
  weights = mw_g$weights, family = gaussian(), group_levels = group_levels
)
stopifnot(max(abs(rate_exact$eigenvalues - rate_w$eigenvalues)) < 1e-12)
cat("1. gaussian: weights = 1/dispersion; rate(weights) == exact rate\n")

## ---------------------------------------------------------------------------
## 2. f7 cross-check at the same point (branches where f7 is correct)
## ---------------------------------------------------------------------------
## glmbfamfunc()$f7 takes a single design matrix and one coefficient vector,
## so check on a single group with that group's b_mode.
j1 <- 1L
rows1 <- which(as.integer(grp) == j1)
Z_1 <- x_re[rows1, , drop = FALSE]
b_1 <- b_mode[j1, ]
y_dummy <- rep(1, length(rows1))   # f7 does not use y

## poisson-log: w_i = lambda_i
mw_p <- two_block_mode_weights(
  x = x_re, block = grp, b_mode = b_mode,
  family = poisson(), group_levels = group_levels
)
stopifnot(all(abs(mw_p$weights - exp(mw_p$eta)) < 1e-12))
f7_p <- glmbfamfunc(poisson())$f7(
  b = b_1, y = y_dummy, x = Z_1, mu = b_1 * 0, P = diag(0, p_re), wt = 1
)
stopifnot(max(abs(mw_p$B_lik[[group_levels[j1]]] - f7_p)) < 1e-10)

## binomial-logit: w_i = wt_i * p_i (1 - p_i), with trial counts wt
wt_tr <- rep(c(5, 10), length.out = length(grp))
mw_b <- two_block_mode_weights(
  x = x_re, block = grp, b_mode = b_mode,
  family = binomial(), wt = wt_tr, group_levels = group_levels
)
p_hat <- plogis(mw_b$eta)
stopifnot(all(abs(mw_b$weights - wt_tr * p_hat * (1 - p_hat)) < 1e-12))
f7_b <- glmbfamfunc(binomial())$f7(
  b = b_1, y = y_dummy, x = Z_1, mu = b_1 * 0, P = diag(0, p_re),
  wt = wt_tr[rows1]
)
stopifnot(max(abs(mw_b$B_lik[[group_levels[j1]]] - f7_b)) < 1e-10)
cat("2. poisson-log and binomial-logit match glmbfamfunc()$f7 at b_mode\n")

## ---------------------------------------------------------------------------
## 3. B_lik brute force and info_total = sum of blocks
## ---------------------------------------------------------------------------
for (j in seq_len(J)) {
  rows <- which(as.integer(grp) == j)
  Z_j <- x_re[rows, , drop = FALSE]
  B_bf <- crossprod(Z_j, Z_j * mw_p$weights[rows])
  stopifnot(max(abs(mw_p$B_lik[[group_levels[j]]] - B_bf)) < 1e-12)
}
stopifnot(max(abs(mw_p$info_total - Reduce(`+`, mw_p$B_lik))) < 1e-12)
cat("3. B_lik matches brute force; info_total = sum(B_lik)\n")

## ---------------------------------------------------------------------------
## 4. probit: generic formula (f7's logistic weights are wrong here)
## ---------------------------------------------------------------------------
mw_pr <- two_block_mode_weights(
  x = x_re, block = grp, b_mode = b_mode,
  family = binomial(link = "probit"), group_levels = group_levels
)
p_pr <- pnorm(mw_pr$eta)
w_pr_closed <- dnorm(mw_pr$eta)^2 / (p_pr * (1 - p_pr))
stopifnot(all(abs(mw_pr$weights - w_pr_closed) < 1e-12))
## and they differ from the logistic weights f7 would have used
stopifnot(max(abs(mw_pr$weights - p_pr * (1 - p_pr))) > 1e-3)
cat("4. probit: weights = phi(eta)^2 / (p(1-p)), not the logistic p(1-p)\n")

## ---------------------------------------------------------------------------
## 5. plumbing: weights -> two_block_rate gives a valid heuristic spectrum
## ---------------------------------------------------------------------------
rate_p <- two_block_rate(
  x = x_re, block = grp, x_hyper = x_hyper,
  prior_list_block1 = list(P = prior_b$P),
  prior_list_block2 = prior_list_block2,
  weights = mw_p$weights, family = poisson(), group_levels = group_levels
)
stopifnot(all(rate_p$eigenvalues >= 0), all(rate_p$eigenvalues < 1))
stopifnot(identical(rate_p$weights_source, "user"))
l_needed <- two_block_l_for_tv(rate_p, 1e-3)
stopifnot(l_needed >= 1L)
print(mw_p)
cat(sprintf("5. plumbing: heuristic lambda* = %.4f, l(1e-3) = %d\n",
            rate_p$lambda_star, l_needed))

## ---------------------------------------------------------------------------
## 6. input validation
## ---------------------------------------------------------------------------
expect_err <- function(expr, pat) {
  res <- tryCatch(expr, error = function(e) e)
  stopifnot(inherits(res, "error"), grepl(pat, conditionMessage(res)))
}
expect_err(
  two_block_mode_weights(x_re, grp, b_mode, family = gaussian(),
                         group_levels = group_levels),
  "dispersion"
)
expect_err(
  two_block_mode_weights(x_re, grp, b_mode[, 1L, drop = FALSE],
                         family = poisson(), group_levels = group_levels),
  "ncol\\(b_mode\\)"
)
expect_err(
  two_block_mode_weights(x_re, grp, b_mode[-1L, ], family = poisson(),
                         group_levels = group_levels),
  "missing rows"
)
expect_err(
  two_block_mode_weights(x_re, grp, b_mode, family = poisson(),
                         wt = rep(1, 3L), group_levels = group_levels),
  "length\\(wt\\)"
)
cat("6. input validation errors as expected\n")

cat("\ntest_two_block_mode_weights.R: all checks passed\n")
