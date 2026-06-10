## Regression tests for the C++ two-block Gibbs driver
## (two_block_rNormal_reg_cpp_export, src/twoBlockGibbs.cpp).
##
## The driver is a port-only migration of the former R loop in
## two_block_rNormal_reg(): mu_all -> Block 1 (block_rNormalReg /
## block_rNormalGLM C++ exports) -> Block 2 (rNormalReg per RE component).
##
## Equivalence policy: R and C++ random number generation are NOT the same
## (the envelope rejection sampler uses a thread-local std::mt19937 seeded
## from std::random_device, independent of R's seed), so individual draws
## can never be compared.  We reconstruct the legacy R loop inline (using
## the still-exported pieces) and compare AVERAGE coefficients across many
## draws with generous tolerances, for gaussian and Poisson Block 1.
## Run: Rscript data-raw/test_two_block_cpp.R

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload.", call. = FALSE)
}
pkgload::load_all(export_all = FALSE)

## ---------------------------------------------------------------------------
## Legacy R loop (verbatim copy of the pre-port two_block_rNormal_reg body)
## ---------------------------------------------------------------------------
ref_two_block <- function(n, y, x, block, x_hyper,
                          prior_list_block1, prior_list_block2,
                          fixef_start, re_names, group_levels,
                          m_convergence, family) {
  is_gaussian <- identical(family$family, "gaussian")
  meta <- glmbayesCore:::.two_block_validate_block1_prior(
    prior_list_block1, family = family
  )
  block1_fn <- if (is_gaussian) block_rNormalReg else block_rNormalGLM

  J <- length(group_levels)
  p_re <- length(re_names)
  fixef_draws <- stats::setNames(
    lapply(re_names, function(k) {
      matrix(NA_real_, nrow = n, ncol = length(fixef_start[[k]]))
    }),
    re_names
  )
  b_sum <- matrix(0, J, p_re)
  b_i <- NULL

  for (i in seq_len(n)) {
    fixef <- fixef_start
    for (m in seq_len(m_convergence)) {
      mu_all <- glmbayesCore:::.two_block_mu_all(
        fixef, x_hyper, re_names, group_levels
      )
      pl1 <- glmbayesCore:::.two_block_block1_prior_list(
        prior_list_block1, mu_all, meta
      )
      args <- list(
        n = 1L, y = y, x = x, block = block,
        prior_list = pl1, offset = NULL, weights = 1
      )
      if (!is_gaussian) {
        args <- c(args, list(family = family, use_parallel = FALSE))
      }
      block_i <- do.call(block1_fn, args)
      b_i <- block_i$coefficients
      colnames(b_i) <- re_names

      fixef_draw <- multi_rNormal_reg(
        n = 1L, y = b_i, x = x_hyper,
        prior_list = prior_list_block2,
        family = gaussian(), progbar = FALSE
      )
      fixef <- stats::setNames(
        lapply(re_names, function(k) fixef_draw[[k]]$coefficients[1L, ]),
        re_names
      )
    }
    for (k in re_names) fixef_draws[[k]][i, ] <- fixef[[k]]
    b_sum <- b_sum + unname(b_i)
  }

  list(
    fixef_means = lapply(fixef_draws, colMeans),
    fixef_vars  = lapply(fixef_draws, function(m) apply(m, 2L, stats::var)),
    b_means     = b_sum / n
  )
}

## Average-coefficient comparison with a generous fixed tolerance.
## (No per-draw or tight-SE checks: RNG streams differ across paths.
## Payload-level bugs shift means by > 0.3 in these setups.)
check_means <- function(m1, m2, tol, label) {
  diff <- abs(m1 - m2)
  if (any(!is.finite(diff)) || any(diff > tol)) {
    stop(label, ": average coefficients differ. max diff = ",
         max(diff), ", allowed = ", tol)
  }
  invisible(max(diff))
}

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

prior_list_block2 <- list(
  `(Intercept)` = list(mu = c(0, 0), Sigma = diag(4, 2L), dispersion = 0.16),
  slope         = list(mu = 0, Sigma = diag(4, 1L), dispersion = 0.16)
)
fixef_start <- list(
  `(Intercept)` = stats::setNames(c(0, 0), colnames(X_int)),
  slope         = stats::setNames(0, colnames(X_slp))
)

n_draw <- 300L
m_conv <- 3L

## ---------------------------------------------------------------------------
## 1. Gaussian Block 1: structure + mean equivalence vs legacy R loop
## ---------------------------------------------------------------------------
y_gauss <- eta + rnorm(length(eta), sd = 0.5)
pl1_gauss <- list(Sigma = diag(0.25, 2L), dispersion = 0.25)

set.seed(303)
fit_g <- two_block_rNormal_reg(
  n = n_draw, y = y_gauss, x = x_re, block = grp,
  x_hyper = x_hyper,
  prior_list_block1 = pl1_gauss,
  prior_list_block2 = prior_list_block2,
  fixef_start = fixef_start,
  m_convergence = m_conv,
  family = gaussian(),
  progbar = FALSE
)

stopifnot(inherits(fit_g, "two_block_rNormal_reg"))
stopifnot(identical(names(fit_g$fixef_draws), re_names))
stopifnot(all(dim(fit_g$fixef_draws[["(Intercept)"]]) == c(n_draw, 2L)))
stopifnot(identical(colnames(fit_g$fixef_draws[["(Intercept)"]]), colnames(X_int)))
stopifnot(nrow(fit_g$coefficients) == n_draw * J)
stopifnot(identical(colnames(fit_g$coefficients)[1:2], c("draw", "grp")))
stopifnot(all(dim(fit_g$b_last) == c(J, 2L)))
stopifnot(identical(rownames(fit_g$b_last), group_levels))
stopifnot(identical(dimnames(fit_g$mu_all_last), list(re_names, group_levels)))
stopifnot(all(is.finite(as.matrix(fit_g$coefficients[, re_names]))))
cat("1a. gaussian: structure OK\n")

set.seed(404)
ref_g <- ref_two_block(
  n = n_draw, y = y_gauss, x = x_re, block = grp, x_hyper = x_hyper,
  prior_list_block1 = pl1_gauss, prior_list_block2 = prior_list_block2,
  fixef_start = fixef_start, re_names = re_names,
  group_levels = group_levels, m_convergence = m_conv, family = gaussian()
)

for (k in re_names) {
  check_means(
    colMeans(fit_g$fixef_draws[[k]]), ref_g$fixef_means[[k]],
    tol = 0.15, label = paste0("gaussian fixef[", k, "]")
  )
}
b_means_cpp <- sapply(seq_along(re_names), function(jj) {
  tapply(fit_g$coefficients[[re_names[jj]]], fit_g$coefficients$grp, mean)
})
d_b <- check_means(unname(b_means_cpp), ref_g$b_means,
                   tol = 0.2, label = "gaussian b means")
cat("1b. gaussian: average fixef + b match legacy R loop (b max diff ",
    format(d_b, digits = 3), ")\n", sep = "")

## ---------------------------------------------------------------------------
## 2. Poisson Block 1 (GLM envelope path): structure + mean equivalence
## ---------------------------------------------------------------------------
y_pois <- rpois(length(eta), exp(0.3 * eta))
pl1_pois <- list(Sigma = diag(0.25, 2L))
n_draw_p <- 600L

set.seed(505)
fit_p <- two_block_rNormal_reg(
  n = n_draw_p, y = y_pois, x = x_re, block = grp,
  x_hyper = x_hyper,
  prior_list_block1 = pl1_pois,
  prior_list_block2 = prior_list_block2,
  fixef_start = fixef_start,
  m_convergence = m_conv,
  family = poisson(),
  use_parallel = FALSE,
  progbar = FALSE
)
stopifnot(inherits(fit_p, "two_block_rNormal_reg"))
stopifnot(nrow(fit_p$coefficients) == n_draw_p * J)
stopifnot(all(is.finite(as.matrix(fit_p$coefficients[, re_names]))))
cat("2a. poisson: structure OK\n")

set.seed(606)
ref_p <- ref_two_block(
  n = n_draw_p, y = y_pois, x = x_re, block = grp, x_hyper = x_hyper,
  prior_list_block1 = pl1_pois, prior_list_block2 = prior_list_block2,
  fixef_start = fixef_start, re_names = re_names,
  group_levels = group_levels, m_convergence = m_conv, family = poisson()
)

for (k in re_names) {
  check_means(
    colMeans(fit_p$fixef_draws[[k]]), ref_p$fixef_means[[k]],
    tol = 0.15, label = paste0("poisson fixef[", k, "]")
  )
}
b_means_cpp_p <- sapply(seq_along(re_names), function(jj) {
  tapply(fit_p$coefficients[[re_names[jj]]], fit_p$coefficients$grp, mean)
})
d_bp <- check_means(unname(b_means_cpp_p), ref_p$b_means,
                    tol = 0.2, label = "poisson b means")
cat("2b. poisson: average fixef + b match legacy R loop (b max diff ",
    format(d_bp, digits = 3), ")\n", sep = "")

## ---------------------------------------------------------------------------
## 3. Error paths: missing Block 2 dispersion; missing Block 1 P/Sigma
## ---------------------------------------------------------------------------
bad_pl2 <- prior_list_block2
bad_pl2[["slope"]]$dispersion <- NULL
err1 <- tryCatch(
  two_block_rNormal_reg(
    n = 1L, y = y_gauss, x = x_re, block = grp, x_hyper = x_hyper,
    prior_list_block1 = pl1_gauss, prior_list_block2 = bad_pl2,
    fixef_start = fixef_start, m_convergence = 1L, family = gaussian(),
    progbar = FALSE
  ),
  error = function(e) conditionMessage(e)
)
stopifnot(is.character(err1), grepl("dispersion", err1))

err2 <- tryCatch(
  two_block_rNormal_reg(
    n = 1L, y = y_gauss, x = x_re, block = grp, x_hyper = x_hyper,
    prior_list_block1 = list(dispersion = 0.25),
    prior_list_block2 = prior_list_block2,
    fixef_start = fixef_start, m_convergence = 1L, family = gaussian(),
    progbar = FALSE
  ),
  error = function(e) conditionMessage(e)
)
stopifnot(is.character(err2), grepl("'P' or 'Sigma'", err2))
cat("3. error paths: OK\n")

cat("\ntest_two_block_cpp.R: all checks passed\n")
