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
##   rate     = d_k * n_prior / 2          (d_k = dispersion guess tau2_k)
##   disp_lower = 1 / qgamma(0.99, shape, rate)   (0.01 quantile of inv-Gamma)
## Hand-picked shape/rate are never used for sampling.
ing_pfamily <- function(d_k, pwt_disp, J, mu = 0, Sigma = diag(100, 1L)) {
  n_prior <- J * pwt_disp / (1 - pwt_disp)
  p_k <- length(mu)
  shape <- (n_prior + 1) / 2 + p_k / 2
  rate  <- d_k * (n_prior / 2)
  dIndependent_Normal_Gamma(
    mu = mu, Sigma = Sigma, shape = shape, rate = rate,
    disp_lower = 1 / stats::qgamma(0.99, shape = shape, rate = rate)
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
dl_int <- pfam_ing[["(Intercept)"]]$prior_list$disp_lower
dl_slp <- pfam_ing[["X1"]]$prior_list$disp_lower

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
stopifnot(all(dd[, 1L] >= dl_int))
stopifnot(all(dd[, 2L] >= dl_slp))
stopifnot(stats::sd(dd[, 1L]) > 0, stats::sd(dd[, 2L]) > 0)
stopifnot(all(is.finite(as.matrix(fit_ing$coefficients[, re_names]))))
for (k in re_names) {
  stopifnot(all(is.finite(fit_ing$fixef_draws[[k]])))
}
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
## 3. Tight-ING vs dNormal equivalence: pwt_disp -> 1 concentrates the
##    inverse-Gamma at d_k = tau2*, so the run must reproduce the
##    dNormal(dispersion = tau2*) posterior
## ---------------------------------------------------------------------------
tau2_star <- 0.16
pwt_disp_tight <- 0.999  # n_prior = 2997: prior dominates, tau2 pinned at d_k

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

## Start near the posterior (school means of b_true) and use enough inner
## steps that both chains forget the start: replicate draws then come from
## the same stationary distribution and the means must agree within MC error.
fixef_start_t <- list(
  "(Intercept)" = stats::setNames(mean(b_true[, 1L]), "(Intercept)"),
  "X1"          = stats::setNames(mean(b_true[, 2L]), "(Intercept)")
)
m_conv_t <- 8L

n_draw_t <- 200L
set.seed(202)
fit_t_ing <- two_block_rNormal_reg_v2(
  n = n_draw_t, y = y, x = x, block = school,
  x_hyper = x_hyper,
  prior_list_block1 = prior_b1_t,
  pfamily_list = pfam_tight,
  fixef_start = fixef_start_t,
  m_convergence = m_conv_t,
  family = gaussian(),
  progbar = FALSE
)
set.seed(202)
fit_t_nrm <- two_block_rNormal_reg_v2(
  n = n_draw_t, y = y, x = x, block = school,
  x_hyper = x_hyper,
  prior_list_block1 = prior_b1_t,
  pfamily_list = pfam_norm,
  fixef_start = fixef_start_t,
  m_convergence = m_conv_t,
  family = gaussian(),
  progbar = FALSE
)

## tau2 draws pinned near tau2* by the tight prior
dd_t <- fit_t_ing$dispersion_fixef_draws
stopifnot(all(abs(dd_t - tau2_star) < 0.05))

for (k in re_names) {
  diff <- abs(colMeans(fit_t_ing$fixef_draws[[k]]) -
              colMeans(fit_t_nrm$fixef_draws[[k]]))
  if (any(!is.finite(diff)) || any(diff > 0.25)) {
    stop("tight-ING vs dNormal fixef means differ [", k,
         "]: max = ", max(diff))
  }
}
d_b <- abs(colMeans(as.matrix(fit_t_ing$coefficients[, re_names])) -
           colMeans(as.matrix(fit_t_nrm$coefficients[, re_names])))
stopifnot(all(d_b < 0.25))
cat("3. tight-ING vs dNormal: posterior means agree (max fixef diff < 0.25)\n")

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

cat("\nAll v2 ING tests passed.\n")
