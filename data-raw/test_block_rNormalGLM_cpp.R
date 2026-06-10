## Regression tests for GLM block C++ path (block_rNormalGLM_cpp_export).
## The export moves block partition + prior payload into C++; per-block
## sampling is the unchanged rNormalGLMBlocks -> rNormalGLM pipeline.
##
## Equivalence policy: posterior modes are deterministic (optim) and must
## match the legacy R-prep payload tightly.  Individual envelope draws are
## NOT expected to match across paths (R- vs C++-assembled payloads can
## differ in the last ulp, and rejection sampling is chaotic in those bits);
## instead, per-block draw MEANS from longer runs must agree within
## Monte Carlo error.
## Run: Rscript data-raw/test_block_rNormalGLM_cpp.R

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload.", call. = FALSE)
}
pkgload::load_all(export_all = FALSE)

tol <- 1e-10

## ---------------------------------------------------------------------------
## 1. Ex_block_rNormalGLM Dobson Poisson example — structure + finite mode
## ---------------------------------------------------------------------------
set.seed(42)
counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
outcome <- gl(3, 1, 9)
treatment <- gl(3, 3)
d.AD <- data.frame(outcome, treatment, counts)

ps <- glmbayesCore::Prior_Setup(counts ~ treatment, family = poisson(), data = d.AD)
y_pois <- ps$y
x_pois <- ps$x
l1 <- ncol(x_pois)

set.seed(101)
out_pois <- block_rNormalGLM(
  n = 1L,
  y = y_pois,
  x = x_pois,
  block = outcome,
  prior_list = list(mu = ps$mu, Sigma = ps$Sigma),
  family = poisson(),
  use_parallel = FALSE
)
stopifnot(inherits(out_pois, "block_rNormalGLM"))
stopifnot(all(dim(out_pois$coefficients) == c(3L, l1)))
stopifnot(all(dim(out_pois$coef.mode) == c(3L, l1)))
stopifnot(all(is.finite(out_pois$coef.mode)))
stopifnot(identical(out_pois$k, 3L))
cat("1. Ex_block_rNormalGLM Dobson Poisson: OK\n")

## ---------------------------------------------------------------------------
## 2. High-level C++ export vs legacy R-prep payload (.rNormalGLMBlocks_cpp)
##    (a) coef.mode is deterministic => tight match.
##    (b) draw means over n_rep repeated n=1 calls agree within MC error
##        (specific draws are not comparable across payload assembly paths).
## ---------------------------------------------------------------------------
prior_list_pois <- list(mu = ps$mu, Sigma = ps$Sigma)
block_info <- glmbayesCore::normalize_block(outcome, length(y_pois))
k <- block_info$k
prior_block <- glmbayesCore:::normalize_prior_for_blocks(
  prior_list = prior_list_pois, prior_lists = NULL,
  block_info = block_info, l1 = l1
)
prior_cpp <- glmbayesCore:::.prior_payload_for_rNormalGLMBlocks_cpp(prior_block, l1, k)
famfunc <- glmbayesCore::glmbfamfunc(poisson())

run_low <- function() {
  glmbayesCore:::.rNormalGLMBlocks_cpp(
    n = 1L, y = y_pois, x = x_pois,
    offset = rep(0, length(y_pois)), wt = rep(1, length(y_pois)),
    dispersion = prior_cpp$dispersion,
    mu = prior_cpp$mu,
    P_blocks = prior_cpp$P_blocks,
    prior_by_block = prior_cpp$prior_by_block,
    row_blocks = block_info$rows,
    f2 = famfunc$f2,
    f3 = famfunc$f3,
    family = "poisson",
    link = "log",
    Gridtype = 2L,
    n_envopt = 1L,
    use_parallel = FALSE,
    use_opencl = FALSE,
    verbose = FALSE
  )
}

low <- run_low()
diff_mode <- max(abs(unname(out_pois$coef.mode) - unname(low$coef.mode)))
if (!is.finite(diff_mode) || diff_mode > 1e-8) {
  stop("block_rNormalGLM vs legacy payload coef.mode differ: max = ", diff_mode)
}
cat("2a. coef.mode high-level vs legacy payload: OK (max diff ",
    format(diff_mode, digits = 3), ")\n", sep = "")

run_hi <- function() {
  block_rNormalGLM(
    n = 1L, y = y_pois, x = x_pois, block = outcome,
    prior_list = prior_list_pois, family = poisson(), use_parallel = FALSE
  )
}

n_rep <- 400L
acc_hi <- matrix(0, 3L, l1)
acc_low <- matrix(0, 3L, l1)
set.seed(2026)
for (r in seq_len(n_rep)) {
  acc_hi <- acc_hi + unname(run_hi()$coefficients)
}
set.seed(2027)
for (r in seq_len(n_rep)) {
  acc_low <- acc_low + unname(run_low()$coefficients)
}
mean_hi <- acc_hi / n_rep
mean_low <- acc_low / n_rep
diff_mean <- max(abs(mean_hi - mean_low))
## Conditional posterior sd per coordinate is O(0.1-0.5) here; with
## n_rep = 400 the MC standard error of a mean difference is ~ sd*sqrt(2/400).
## 0.15 gives ample slack while still catching payload-level errors
## (the Sigma->P transpose bug produced mean shifts > 0.3).
if (!is.finite(diff_mean) || diff_mean > 0.15) {
  stop("block_rNormalGLM vs legacy payload draw means differ: max = ", diff_mean)
}
cat("2b. Draw means over ", n_rep, " reps: OK (max diff ",
    format(diff_mean, digits = 3), ")\n", sep = "")

## ---------------------------------------------------------------------------
## 3. Family/link smoke: binomial (logit, probit, cloglog) and Gamma (log)
## ---------------------------------------------------------------------------
set.seed(7)
n_grp <- 3L
n_per <- 30L
grp <- factor(rep(seq_len(n_grp), each = n_per))
xb <- cbind(1, rnorm(n_grp * n_per))
colnames(xb) <- c("(Intercept)", "X1")
eta <- 0.3 + 0.8 * xb[, 2]

pl_smoke <- list(mu = rep(0, 2L), Sigma = diag(4, 2L))

for (lnk in c("logit", "probit", "cloglog")) {
  fam <- binomial(link = lnk)
  y_bin <- rbinom(length(eta), 1L, fam$linkinv(eta))
  set.seed(202)
  out_b <- block_rNormalGLM(
    n = 1L, y = y_bin, x = xb, block = grp,
    prior_list = pl_smoke, family = fam, use_parallel = FALSE
  )
  stopifnot(all(dim(out_b$coefficients) == c(n_grp, 2L)))
  stopifnot(all(is.finite(out_b$coef.mode)))
  cat("3. binomial/", lnk, ": OK\n", sep = "")
}

y_gam <- rgamma(length(eta), shape = 2, rate = 2 / exp(eta))
set.seed(303)
out_g <- block_rNormalGLM(
  n = 1L, y = y_gam, x = xb, block = grp,
  prior_list = c(pl_smoke, list(dispersion = 0.5)),
  family = Gamma(link = "log"), use_parallel = FALSE
)
stopifnot(all(dim(out_g$coefficients) == c(n_grp, 2L)))
stopifnot(all(is.finite(out_g$coef.mode)))
cat("3. Gamma/log: OK\n")

## ---------------------------------------------------------------------------
## 4. Per-block prior means (mu as l1 x k matrix) — the lmebayes Block 1 shape
## ---------------------------------------------------------------------------
mu_mat <- matrix(c(0, 0, 0.2, 0.1, -0.2, -0.1), nrow = 2L, ncol = n_grp)
y_pois2 <- rpois(length(eta), exp(eta))
set.seed(404)
out_pb <- block_rNormalGLM(
  n = 1L, y = y_pois2, x = xb, block = grp,
  prior_list = list(mu = mu_mat, Sigma = diag(4, 2L)),
  family = poisson(), use_parallel = FALSE
)
stopifnot(all(dim(out_pb$coefficients) == c(n_grp, 2L)))
stopifnot(length(out_pb$prior_lists) == n_grp)
mu_back <- vapply(out_pb$prior_lists, function(pl) pl$mu, numeric(2L))
stopifnot(max(abs(mu_back - mu_mat)) < tol)
cat("4. Per-block prior means (mu matrix): OK\n")

## ---------------------------------------------------------------------------
## 5. Error paths: gaussian family rejected; bad link rejected
## ---------------------------------------------------------------------------
err1 <- tryCatch(
  block_rNormalGLM(
    n = 1L, y = y_pois2, x = xb, block = grp,
    prior_list = pl_smoke, family = gaussian()
  ),
  error = function(e) conditionMessage(e)
)
stopifnot(is.character(err1), grepl("GLM envelope path", err1))

err2 <- tryCatch(
  block_rNormalGLM(
    n = 1L, y = y_pois2, x = xb, block = grp,
    prior_list = pl_smoke, family = poisson(link = "identity")
  ),
  error = function(e) conditionMessage(e)
)
stopifnot(is.character(err2), grepl("not available for family", err2))
cat("5. Error paths (gaussian, bad link): OK\n")

cat("\ntest_block_rNormalGLM_cpp.R: all checks passed\n")
