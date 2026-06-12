## Tests for dIndependent_Normal_Gamma Block 2 priors in the v2 two-block
## Gibbs driver (two_block_rNormal_reg_v2_cpp_export, src/twoBlockGibbs.cpp).
##
## Block 2 ING components make a joint (gamma_k, tau2_k) draw via
## rIndepNormalGammaReg (the same likelihood-subgradient envelope sampler
## used by rglmb with an ING pfamily); the sampled tau2_k feeds back into
## the Block 1 prior precision on the next inner step.
##
## Design adapted from the schools example in test_block_rNormalReg_cpp.R:
## 20 schools, random intercept + slope, gaussian measurement model.
##
## The envelope sampler uses its own RNG stream, so checks are at the level
## of bounds, structure, and posterior means (not per-draw equality).
## Run: Rscript data-raw/test_two_block_v2_ing.R

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload.", call. = FALSE)
}
pkgload::load_all(export_all = FALSE)

## ---------------------------------------------------------------------------
## Schools toy data (adapted from test_block_rNormalReg_cpp.R, scaled up to
## J = 20 schools: the Block 2 hyper-regression has J observations, so J must
## comfortably exceed the number of hyper-predictors and leave residual
## degrees of freedom to inform tau2)
## ---------------------------------------------------------------------------
set.seed(42)
n_schools <- 20L
n_per     <- 10L
school    <- factor(rep(seq_len(n_schools), each = n_per))
x         <- cbind(1, rnorm(n_schools * n_per))
colnames(x) <- c("(Intercept)", "X1")
## School coefficients: intercepts ~ N(5, tau2 = 4), slopes ~ N(0.2, tau2 = 0.2)
b_true    <- cbind(
  rnorm(n_schools, mean = 5,   sd = 2),
  rnorm(n_schools, mean = 0.2, sd = sqrt(0.2))
)
sigma2    <- 1.5
y         <- rowSums(x * b_true[as.integer(school), ]) +
  rnorm(nrow(x), sd = sqrt(sigma2))
re_names  <- colnames(x)

## Group-level designs: intercept-only hyper regression for both components.
x_hyper <- list(
  "(Intercept)" = matrix(1, n_schools, 1L,
                         dimnames = list(NULL, "(Intercept)")),
  "X1"          = matrix(1, n_schools, 1L,
                         dimnames = list(NULL, "(Intercept)"))
)
fixef_start <- list(
  "(Intercept)" = stats::setNames(0, "(Intercept)"),
  "X1"          = stats::setNames(0, "(Intercept)")
)

## Block 1 prior: plug-in RE variances on the diagonal (overridden per sweep
## for ING components), measurement dispersion fixed at sigma2.
tau2_plug <- c(4, 0.2)
prior_b1 <- list(Sigma = diag(tau2_plug, 2L), dispersion = sigma2, ddef = FALSE)

## ING priors are ALWAYS calibrated from a pwt_disp choice, mirroring
## lmebayes::pfamily_list() on a Prior_Setup_lmebayes object:
##   n_prior  = J * pwt_disp / (1 - pwt_disp)
##   shape    = (n_prior + 1) / 2 + p_k / 2
##   rate     = d_k * (n_prior + p_k - 1) / 2   (d_k = dispersion guess tau2_k;
##              = d_k * (shape - 1), so the inv-Gamma prior mean is exactly d_k)
##   disp_lower = 1 / qgamma(0.99, shape, rate)   (0.01 quantile of inv-Gamma)
##   disp_upper = 1 / qgamma(0.01, shape, rate)   (0.99 quantile of inv-Gamma)
## Both bounds are required for sampling: the tau2_k truncation window is
## then the central 98% prior-mass interval, fixed across Gibbs sweeps.
## Hand-picked shape/rate are never used for sampling.
ing_pfamily <- function(d_k, pwt_disp, J, mu = 0, Sigma = diag(100, 1L)) {
  n_prior <- J * pwt_disp / (1 - pwt_disp)
  p_k <- length(mu)
  shape <- (n_prior + 1) / 2 + p_k / 2
  rate  <- d_k * (n_prior + p_k - 1) / 2
  dIndependent_Normal_Gamma(
    mu = mu, Sigma = Sigma, shape = shape, rate = rate,
    disp_lower = 1 / stats::qgamma(0.99, shape = shape, rate = rate),
    disp_upper = 1 / stats::qgamma(0.01, shape = shape, rate = rate)
  )
}

n_draw <- 100L
m_conv <- 2L

## ---------------------------------------------------------------------------
## 1. All-ING run: structure, truncation bounds, tau2 actually varies
##    pwt_disp = 0.5 (prior and the J school-level observations get equal
##    weight); dispersion guesses d_k = tau2_plug.
## ---------------------------------------------------------------------------
pwt_disp <- 0.5
pfam_ing <- list(
  "(Intercept)" = ing_pfamily(tau2_plug[1L], pwt_disp, n_schools),
  "X1"          = ing_pfamily(tau2_plug[2L], pwt_disp, n_schools)
)
set.seed(101)
fit_ing <- two_block_rNormal_reg_v2(
  n = n_draw, y = y, x = x, block = school,
  x_hyper = x_hyper,
  prior_list_block1 = prior_b1,
  pfamily_list = pfam_ing,
  fixef_start = fixef_start,
  m_convergence = m_conv,
  family = gaussian(),
  progbar = FALSE
)

stopifnot(inherits(fit_ing, "two_block_rNormal_reg_v2"))
dd <- fit_ing$dispersion_fixef_draws
stopifnot(is.matrix(dd), nrow(dd) == n_draw, ncol(dd) == 2L)
stopifnot(identical(colnames(dd), re_names))
stopifnot(all(is.finite(dd)), all(dd > 0))
## Both bounds are supplied, so every tau2 draw must lie inside the fixed
## truncation window [disp_lower, disp_upper] (renormalized inverse-CDF).
for (j in seq_along(re_names)) {
  pr_j <- pfam_ing[[re_names[j]]]$prior_list
  stopifnot(
    all(dd[, j] >= pr_j$disp_lower),
    all(dd[, j] <= pr_j$disp_upper)
  )
}
stopifnot(stats::sd(dd[, 1L]) > 0, stats::sd(dd[, 2L]) > 0)
stopifnot(all(is.finite(as.matrix(fit_ing$coefficients[, re_names]))))
for (k in re_names) {
  stopifnot(all(is.finite(fit_ing$fixef_draws[[k]])))
}
## Candidate counts: every stored draw needed at least one envelope candidate
## per inner sweep for each (all-ING) component.
it <- fit_ing$iters_fixef_draws
stopifnot(
  is.matrix(it), nrow(it) == n_draw, ncol(it) == 2L,
  identical(colnames(it), re_names),
  all(is.finite(it)), all(it >= m_conv)
)
cat("   mean candidates per accepted draw: ",
    paste(sprintf("%s = %.2f", re_names, colMeans(it) / m_conv),
          collapse = ", "), "\n", sep = "")
cat("1. all-ING run: structure + bounds OK (tau2 means: ",
    paste(sprintf("%s=%.3g", re_names, colMeans(dd)), collapse = ", "),
    ")\n", sep = "")

## Posterior sanity: gamma for the intercept component should sit near the
## average school intercept (b_true column 1 mean = 5).
g_int <- mean(fit_ing$fixef_draws[["(Intercept)"]])
if (abs(g_int - 5) > 1.5) {
  stop("ING intercept gamma mean far from truth: ", g_int)
}
cat("2. posterior location sane: gamma_int mean = ",
    format(g_int, digits = 4), "\n", sep = "")

## ---------------------------------------------------------------------------
## 3. Prior-vs-data guard: pwt_disp > 0.5 implies n_prior > J, which the
##    dispersion envelope cannot support (log-tilt capped at n_w/2 = J/2;
##    Remark 4.1.3 of the ING vignette).  Such calls must be rejected up
##    front.  [Historical note: pwt_disp = 0.999 at small J used to produce
##    biased gamma draws because the binding cap silently degraded the
##    envelope -- the guard makes that regime unreachable.]
## ---------------------------------------------------------------------------
tau2_star <- 0.16
pwt_disp_tight <- 0.999  # n_prior = 999 * J >> J: prior-dominated, illegal

pfam_tight <- list(
  "(Intercept)" = ing_pfamily(tau2_star, pwt_disp_tight, n_schools),
  "X1"          = ing_pfamily(tau2_star, pwt_disp_tight, n_schools)
)
pfam_norm <- list(
  "(Intercept)" = dNormal(mu = 0, Sigma = diag(100, 1L),
                          dispersion = tau2_star),
  "X1"          = dNormal(mu = 0, Sigma = diag(100, 1L),
                          dispersion = tau2_star)
)
prior_b1_t <- list(Sigma = diag(tau2_star, 2L), dispersion = sigma2,
                   ddef = FALSE)

err_t <- tryCatch(
  two_block_rNormal_reg_v2(
    n = 5L, y = y, x = x, block = school,
    x_hyper = x_hyper,
    prior_list_block1 = prior_b1_t,
    pfamily_list = pfam_tight,
    fixef_start = fixef_start,
    m_convergence = 2L,
    family = gaussian(),
    progbar = FALSE
  ),
  error = function(e) conditionMessage(e)
)
stopifnot(is.character(err_t), grepl("n_prior <= J", err_t, fixed = TRUE))

## Boundary case pwt_disp = 0.5 (n_prior = J) must remain legal: already
## exercised by sections 1-2 above (pfam_ing uses pwt_disp = 0.5).
cat("3. prior-vs-data guard: pwt_disp = 0.999 rejected (n_prior > J)\n")

## ---------------------------------------------------------------------------
## 4. Mixed priors: ING intercept + dNormal slope
## ---------------------------------------------------------------------------
pfam_mixed <- list(
  "(Intercept)" = pfam_ing[["(Intercept)"]],
  "X1"          = pfam_norm[["X1"]]
)
set.seed(303)
fit_mix <- two_block_rNormal_reg_v2(
  n = 50L, y = y, x = x, block = school,
  x_hyper = x_hyper,
  prior_list_block1 = prior_b1,
  pfamily_list = pfam_mixed,
  fixef_start = fixef_start,
  m_convergence = m_conv,
  family = gaussian(),
  progbar = FALSE
)
dd_m <- fit_mix$dispersion_fixef_draws
stopifnot(stats::sd(dd_m[, "(Intercept)"]) > 0)        # ING: varies
stopifnot(all(dd_m[, "X1"] == tau2_star))              # dNormal: fixed
cat("4. mixed ING + dNormal: OK\n")

## ---------------------------------------------------------------------------
## 5. two_block_rate_v2 with ING components uses the disp_lower plug-in
## ---------------------------------------------------------------------------
r_ing <- two_block_rate_v2(
  x = x, block = school, x_hyper = x_hyper,
  prior_list_block1 = prior_b1,
  pfamily_list = pfam_ing,
  family = gaussian()
)
stopifnot(is.finite(r_ing$lambda_star),
          r_ing$lambda_star >= 0, r_ing$lambda_star < 1)
cat("5. two_block_rate_v2 (ING plug-in): lambda* = ",
    format(r_ing$lambda_star, digits = 6), "\n", sep = "")

## ---------------------------------------------------------------------------
## 6. One-sided ING (disp_lower only) is rejected for sampling: without
##    disp_upper the envelope would fall back to a per-sweep surrogate
##    posterior window, making the truncation state-dependent.  The
##    calibration-only path (two_block_rate_v2) still accepts it.
## ---------------------------------------------------------------------------
pr1 <- pfam_ing[["(Intercept)"]]$prior_list
pf_onesided <- pfam_ing
pf_onesided[["(Intercept)"]] <- dIndependent_Normal_Gamma(
  mu = pr1$mu, Sigma = pr1$Sigma, shape = pr1$shape, rate = pr1$rate,
  disp_lower = as.numeric(pr1$disp_lower)
)
err_os <- tryCatch(
  two_block_rNormal_reg_v2(
    n = 5L, y = y, x = x, block = school,
    x_hyper = x_hyper,
    prior_list_block1 = prior_b1,
    pfamily_list = pf_onesided,
    fixef_start = fixef_start,
    m_convergence = 2L,
    family = gaussian(),
    progbar = FALSE
  ),
  error = function(e) conditionMessage(e)
)
stopifnot(is.character(err_os), grepl("disp_upper", err_os, fixed = TRUE))
r_os <- two_block_rate_v2(
  x = x, block = school, x_hyper = x_hyper,
  prior_list_block1 = prior_b1,
  pfamily_list = pf_onesided,
  family = gaussian()
)
stopifnot(isTRUE(all.equal(r_os$lambda_star, r_ing$lambda_star)))
cat("6. one-sided ING rejected for sampling, accepted for calibration: OK\n")

cat("\nAll v2 ING tests passed.\n")
