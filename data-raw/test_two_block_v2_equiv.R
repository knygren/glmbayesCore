## Equivalence gate for the v2 (pfamily-based) two-block Gibbs driver
## (two_block_rNormal_reg_v2_cpp_export, src/twoBlockGibbs.cpp).
##
## v2 takes Block 2 priors as pfamily objects (dNormal /
## dIndependent_Normal_Gamma) instead of bare prior lists.  With dNormal
## priors throughout, v2 must produce draws IDENTICAL to v1 under the same
## seed for gaussian Block 1 (all randomness flows through R's RNG).  For
## non-Gaussian Block 1 the GLM envelope sampler uses its own RNG stream
## (std::mt19937 seeded from std::random_device), so only average
## coefficients are compared there.
## Run: Rscript data-raw/test_two_block_v2_equiv.R

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload.", call. = FALSE)
}
pkgload::load_all(export_all = FALSE)

tol_exact <- 1e-12

## ---------------------------------------------------------------------------
## Shared toy design: J groups, random intercept + slope, level-2 covariate
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

## v1 contract: bare Block 2 prior lists
prior_list_block2 <- list(
  `(Intercept)` = list(mu = c(0, 0), Sigma = diag(4, 2L), dispersion = 0.16),
  slope         = list(mu = 0, Sigma = diag(4, 1L), dispersion = 0.16)
)
## v2 contract: the same priors as dNormal pfamily objects
pfam_list <- list(
  `(Intercept)` = dNormal(mu = c(0, 0), Sigma = diag(4, 2L), dispersion = 0.16),
  slope         = dNormal(mu = 0, Sigma = diag(4, 1L), dispersion = 0.16)
)
fixef_start <- list(
  `(Intercept)` = stats::setNames(c(0, 0), colnames(X_int)),
  slope         = stats::setNames(0, colnames(X_slp))
)

## ---------------------------------------------------------------------------
## 1. Gaussian Block 1: same seed => identical draws (v1 vs v2)
## ---------------------------------------------------------------------------
y_gauss <- eta + rnorm(length(eta), sd = 0.5)
pl1_gauss <- list(Sigma = diag(0.25, 2L), dispersion = 0.25)

n_draw <- 50L
m_conv <- 3L

set.seed(515)
fit_v1 <- two_block_rNormal_reg(
  n = n_draw, y = y_gauss, x = x_re, block = grp,
  x_hyper = x_hyper,
  prior_list_block1 = pl1_gauss,
  prior_list_block2 = prior_list_block2,
  fixef_start = fixef_start,
  m_convergence = m_conv,
  family = gaussian(),
  progbar = FALSE
)

set.seed(515)
fit_v2 <- two_block_rNormal_reg_v2(
  n = n_draw, y = y_gauss, x = x_re, block = grp,
  x_hyper = x_hyper,
  prior_list_block1 = pl1_gauss,
  pfamily_list = pfam_list,
  fixef_start = fixef_start,
  m_convergence = m_conv,
  family = gaussian(),
  progbar = FALSE
)

stopifnot(inherits(fit_v2, "two_block_rNormal_reg_v2"))
stopifnot(inherits(fit_v2, "two_block_rNormal_reg"))

for (k in re_names) {
  d <- max(abs(fit_v1$fixef_draws[[k]] - fit_v2$fixef_draws[[k]]))
  if (!is.finite(d) || d > tol_exact) {
    stop("gaussian fixef_draws[", k, "] differ between v1 and v2: max = ", d)
  }
}
d_b <- max(abs(as.matrix(fit_v1$coefficients[, re_names]) -
               as.matrix(fit_v2$coefficients[, re_names])))
if (!is.finite(d_b) || d_b > tol_exact) {
  stop("gaussian b draws differ between v1 and v2: max = ", d_b)
}
stopifnot(identical(fit_v1$b_last, fit_v2$b_last))
cat("1. gaussian v1 vs v2: identical draws OK (max diff ",
    format(max(d_b), digits = 3), ")\n", sep = "")

## ---------------------------------------------------------------------------
## 2. dispersion_fixef_draws: present, constant at the dNormal dispersions
## ---------------------------------------------------------------------------
dd <- fit_v2$dispersion_fixef_draws
stopifnot(is.matrix(dd), nrow(dd) == n_draw, ncol(dd) == length(re_names))
stopifnot(identical(colnames(dd), re_names))
stopifnot(all(dd[, "(Intercept)"] == 0.16), all(dd[, "slope"] == 0.16))
## iters_fixef_draws: total Block 2 candidates per stored draw (summed over
## inner sweeps).  dNormal components are conjugate -- exactly 1 candidate
## per sweep, so every entry equals m_convergence.
it <- fit_v2$iters_fixef_draws
stopifnot(is.matrix(it), nrow(it) == n_draw, ncol(it) == length(re_names))
stopifnot(identical(colnames(it), re_names))
stopifnot(all(it == m_conv))
cat("2. dispersion_fixef_draws + iters_fixef_draws: OK\n")

## ---------------------------------------------------------------------------
## 3. Poisson Block 1: average-coefficient equivalence (RNG streams differ)
## ---------------------------------------------------------------------------
y_pois <- rpois(length(eta), lambda = exp(0.3 * eta))
pl1_pois <- list(Sigma = diag(0.25, 2L))

n_draw_p <- 200L
set.seed(616)
fit_v1p <- two_block_rNormal_reg(
  n = n_draw_p, y = y_pois, x = x_re, block = grp,
  x_hyper = x_hyper,
  prior_list_block1 = pl1_pois,
  prior_list_block2 = prior_list_block2,
  fixef_start = fixef_start,
  m_convergence = m_conv,
  family = poisson(),
  progbar = FALSE
)
set.seed(616)
fit_v2p <- two_block_rNormal_reg_v2(
  n = n_draw_p, y = y_pois, x = x_re, block = grp,
  x_hyper = x_hyper,
  prior_list_block1 = pl1_pois,
  pfamily_list = pfam_list,
  fixef_start = fixef_start,
  m_convergence = m_conv,
  family = poisson(),
  progbar = FALSE
)
for (k in re_names) {
  diff <- abs(colMeans(fit_v1p$fixef_draws[[k]]) -
              colMeans(fit_v2p$fixef_draws[[k]]))
  if (any(!is.finite(diff)) || any(diff > 0.15)) {
    stop("poisson fixef means differ between v1 and v2 [", k,
         "]: max = ", max(diff))
  }
}
cat("3. poisson v1 vs v2: average coefficients OK\n")

## ---------------------------------------------------------------------------
## 4. ING component: mixed dNormal + ING smoke run (full coverage in
##    data-raw/test_two_block_v2_ing.R).  ING hyperparameters are always
##    calibrated from a pwt_disp choice (pfamily_list() formulas), never
##    hand-picked: n_prior = J*pwt/(1-pwt), shape = (n_prior+1)/2 + p_k/2,
##    rate = d_k*(n_prior+p_k-1)/2 (= d_k*(shape-1), mean-matched),
##    disp_lower/disp_upper = central 98% prior-mass window
##    1/qgamma(0.99|0.01, shape, rate) (both required for sampling).
## ---------------------------------------------------------------------------
pwt_disp_s <- 0.5
n_prior_s <- J * pwt_disp_s / (1 - pwt_disp_s)
shape_s <- (n_prior_s + 1) / 2 + 1 / 2
rate_s  <- 0.16 * (shape_s - 1)   # d_k = 0.16 (the dNormal dispersion); b_0 = d_k*(shape-1)
dl_s    <- 1 / stats::qgamma(0.99, shape = shape_s, rate = rate_s)
du_s    <- 1 / stats::qgamma(0.01, shape = shape_s, rate = rate_s)
pfam_ing <- pfam_list
pfam_ing$slope <- dIndependent_Normal_Gamma(
  mu = 0, Sigma = diag(4, 1L), shape = shape_s, rate = rate_s,
  disp_lower = dl_s, disp_upper = du_s
)
fit_ing <- two_block_rNormal_reg_v2(
  n = 5L, y = y_gauss, x = x_re, block = grp,
  x_hyper = x_hyper,
  prior_list_block1 = pl1_gauss,
  pfamily_list = pfam_ing,
  fixef_start = fixef_start,
  m_convergence = 2L,
  family = gaussian(),
  progbar = FALSE
)
dd_ing <- fit_ing$dispersion_fixef_draws
stopifnot(all(is.finite(dd_ing)), all(dd_ing > 0))
stopifnot(all(dd_ing[, "(Intercept)"] == 0.16))  # dNormal: fixed
stopifnot(all(dd_ing[, "slope"] >= dl_s))        # ING: within fixed window
stopifnot(all(dd_ing[, "slope"] <= du_s))
## Candidate counts: dNormal column = m_convergence exactly; ING column is
## at least m_convergence (>= 1 envelope candidate per accepted draw).
it_ing <- fit_ing$iters_fixef_draws
stopifnot(all(it_ing[, "(Intercept)"] == 2L))
stopifnot(all(it_ing[, "slope"] >= 2L))
cat("4. mixed dNormal + ING smoke run: OK\n")

## Validation: ING without disp_lower is rejected up front
pfam_bad <- pfam_list
pfam_bad$slope <- dIndependent_Normal_Gamma(
  mu = 0, Sigma = diag(4, 1L), shape = 3, rate = 2
)
res2 <- tryCatch(
  two_block_rNormal_reg_v2(
    n = 2L, y = y_gauss, x = x_re, block = grp,
    x_hyper = x_hyper,
    prior_list_block1 = pl1_gauss,
    pfamily_list = pfam_bad,
    fixef_start = fixef_start,
    m_convergence = 1L,
    family = gaussian(),
    progbar = FALSE
  ),
  error = function(e) conditionMessage(e)
)
stopifnot(is.character(res2), grepl("disp_lower", res2))
cat("5. ING validation (missing disp_lower): OK\n")

## ---------------------------------------------------------------------------
## 6. two_block_rate_v2 matches two_block_rate on the unwrapped priors
## ---------------------------------------------------------------------------
r1 <- two_block_rate(
  x = x_re, block = grp, x_hyper = x_hyper,
  prior_list_block1 = pl1_gauss,
  prior_list_block2 = prior_list_block2,
  family = gaussian()
)
r2 <- two_block_rate_v2(
  x = x_re, block = grp, x_hyper = x_hyper,
  prior_list_block1 = pl1_gauss,
  pfamily_list = pfam_list,
  family = gaussian()
)
stopifnot(isTRUE(all.equal(r1$lambda_star, r2$lambda_star, tolerance = 1e-12)))
stopifnot(isTRUE(all.equal(r1$eigenvalues, r2$eigenvalues, tolerance = 1e-12)))
cat("6. two_block_rate_v2: OK (lambda* = ",
    format(r2$lambda_star, digits = 6), ")\n", sep = "")

cat("\nAll v2 equivalence tests passed.\n")
