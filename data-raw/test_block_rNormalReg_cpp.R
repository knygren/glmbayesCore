## Regression tests for Gaussian block C++ path (block_rNormalReg_cpp_export).
## Run: Rscript data-raw/test_block_rNormalReg_cpp.R

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload.", call. = FALSE)
}
pkgload::load_all(export_all = FALSE)

tol <- 1e-10

## ---------------------------------------------------------------------------
## 1. Ex_block_rNormalReg toy example — coef.mode deterministic
## ---------------------------------------------------------------------------
set.seed(42)
n_schools <- 3L
n_per     <- 10L
school    <- rep(seq_len(n_schools), each = n_per)
x         <- cbind(1, rnorm(n_schools * n_per))
colnames(x) <- c("(Intercept)", "X1")
b_true    <- matrix(c(5, 0.5, 3, -0.2, 7, 0.3), nrow = n_schools, byrow = TRUE)
sigma2    <- 1.5
y         <- rowSums(x * b_true[school, ]) + rnorm(nrow(x), sd = sqrt(sigma2))
l1        <- ncol(x)
prior_list <- list(
  mu         = rep(0, l1),
  Sigma      = diag(100, l1),
  dispersion = sigma2,
  ddef       = FALSE
)

out <- block_rNormalReg(
  n = 1L, y = y, x = x, block = school, prior_list = prior_list
)
stopifnot(inherits(out, "block_rNormalReg"))
stopifnot(all(dim(out$coefficients) == c(n_schools, l1)))
stopifnot(all(is.finite(out$coef.mode)))
cat("1. Ex_block_rNormalReg toy: OK\n")

## ---------------------------------------------------------------------------
## 2. High-level C++ export vs low-level .rNormalRegBlocks_cpp (same payload)
## ---------------------------------------------------------------------------
block_info <- glmbayesCore::normalize_block(school, length(y))
k <- block_info$k
prior_block <- glmbayesCore:::normalize_prior_for_blocks(
  prior_list = prior_list, prior_lists = NULL,
  block_info = block_info, l1 = l1
)
prior_cpp <- glmbayesCore:::.prior_payload_for_rNormalGLMBlocks_cpp(prior_block, l1, k)

low <- glmbayesCore:::.rNormalRegBlocks_cpp(
  n = 1L, y = y, x = x,
  offset = rep(0, length(y)), wt = rep(1, length(y)),
  dispersion = prior_cpp$dispersion,
  mu = prior_cpp$mu,
  P_blocks = prior_cpp$P_blocks,
  prior_by_block = prior_cpp$prior_by_block,
  row_blocks = block_info$rows,
  f2 = glmbayesCore::glmbfamfunc(stats::gaussian())$f2,
  f3 = glmbayesCore::glmbfamfunc(stats::gaussian())$f3,
  Gridtype = 2L
)

diff_mode <- max(abs(out$coef.mode - low$coef.mode))
if (!is.finite(diff_mode) || diff_mode > tol) {
  stop("block_rNormalReg vs .rNormalRegBlocks_cpp coef.mode differ: max = ", diff_mode)
}
cat("2. High-level vs low-level blocks: OK (max diff ", format(diff_mode, digits = 3), ")\n", sep = "")

## ---------------------------------------------------------------------------
## 3. normalize_block equivalence on factor / integer / list inputs
## ---------------------------------------------------------------------------
l2 <- length(y)
blk_factor <- glmbayesCore::normalize_block(factor(school), l2)
blk_int    <- glmbayesCore::normalize_block(as.integer(school), l2)
blk_list   <- glmbayesCore::normalize_block(split(seq_len(l2), school), l2)
stopifnot(identical(blk_factor$k, blk_int$k))
stopifnot(identical(blk_factor$l2_blocks, blk_int$l2_blocks))
cat("3. normalize_block partition shapes: OK\n")

## ---------------------------------------------------------------------------
## 4. two_block_rNormal_reg smoke run (Gaussian) — deterministic with seed
## ---------------------------------------------------------------------------
x_hyper <- list("(Intercept)" = matrix(1, n_schools, 1), "X1" = matrix(0, n_schools, 1))
prior_b1 <- list(P = diag(0.01, l1), dispersion = sigma2, ddef = FALSE)
prior_b2 <- list(
  "(Intercept)" = list(mu = 0, Sigma = diag(100, 1L), dispersion = 1),
  "X1"          = list(mu = 0, Sigma = diag(100, 1L), dispersion = 1)
)
fixef_start <- list("(Intercept)" = 0, "X1" = 0)

run_once <- function(seed) {
  set.seed(seed)
  two_block_rNormal_reg(
    n = 2L, y = y, x = x, block = factor(school),
    x_hyper = x_hyper,
    prior_list_block1 = prior_b1,
    prior_list_block2 = prior_b2,
    fixef_start = fixef_start,
    m_convergence = 1L,
    family = gaussian(),
    progbar = FALSE,
    seed = NULL
  )
}

r1 <- run_once(123L)
r2 <- run_once(123L)
stopifnot(all(is.finite(r1$b_last)))
d_coef <- max(abs(r1$b_last - r2$b_last), na.rm = TRUE)
if (!is.finite(d_coef) || d_coef > tol) {
  stop("two_block_rNormal_reg not reproducible with fixed seed: max = ", d_coef)
}
cat("4. two_block reproducibility: OK (max diff ", format(d_coef, digits = 3), ")\n", sep = "")

## ---------------------------------------------------------------------------
## 5. big_word_club — one Block 1 step via block_rNormalReg (optional)
## ---------------------------------------------------------------------------
if (requireNamespace("bayesrules", quietly = TRUE) &&
    requireNamespace("lme4", quietly = TRUE)) {
  data(big_word_club, package = "bayesrules")
  dat <- big_word_club
  dat$school_id <- factor(dat$school_id)
  dat <- subset(
    dat,
    !is.na(score_ppvt) &
      !is.na(invalid_ppvt) & invalid_ppvt == 0L &
      complete.cases(dat[, c("score_ppvt", "distracted_a1", "distracted_ppvt",
                             "private_school", "title1", "free_reduced_lunch",
                             "school_id")])
  )
  if (requireNamespace("lmebayes", quietly = TRUE)) {
    form_lmer <- score_ppvt ~
      private_school + title1 + free_reduced_lunch +
      distracted_a1 + distracted_ppvt +
      free_reduced_lunch:distracted_a1 +
      (1 + distracted_ppvt + distracted_a1 || school_id)
    ctrl <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
    ps <- lmebayes::Prior_Setup_lmebayes(form_lmer, data = dat, pwt = 0.01)
    design <- lmebayes::model_setup(form_lmer, data = dat, control = ctrl)
    fixef <- lapply(ps$prior_list, `[[`, "mu_fixef")
    names(fixef) <- design$re_coef_names
    mu_all <- as.matrix(glmbayesCore::build_mu_all(design, fixef)$mu_all)
    set.seed(99)
    b1 <- block_rNormalReg(
      n = 1L,
      y = design$y,
      x = design$Z,
      block = design$groups,
      prior_list = list(
        mu = mu_all,
        Sigma = ps$Sigma_ranef,
        dispersion = ps$dispersion_ranef,
        ddef = FALSE
      )
    )
    stopifnot(nrow(b1$coefficients) == nlevels(design$groups))
    stopifnot(all(is.finite(b1$coef.mode)))
    cat("5. big_word_club Block 1: OK (", nrow(b1$coefficients), " schools)\n", sep = "")
  } else {
    cat("5. big_word_club Block 1: skipped (lmebayes not installed)\n")
  }
} else {
  cat("5. big_word_club Block 1: skipped (bayesrules/lme4 not installed)\n")
}

cat("\ntest_block_rNormalReg_cpp.R: all checks passed\n")
