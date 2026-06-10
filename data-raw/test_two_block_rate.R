## Tests for two_block_rate(): Remark 8 eigenvalues for the two-block Gibbs
## sampler (Nygren 2020).
##
## Strategy: build the FULL dense joint precision over (gamma, b) by brute
## force, compute A = P11^{-1/2} P12 P22^{-1} P21 P11^{-1/2} densely, and
## check that the fast per-group accumulation in two_block_rate() reproduces
## the spectrum to 1e-10.  Also checks block-swap invariance of the nonzero
## spectrum and basic sanity of the weights (non-Gaussian heuristic) path.
## Run: Rscript data-raw/test_two_block_rate.R

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload.", call. = FALSE)
}
pkgload::load_all(export_all = FALSE)

## ---------------------------------------------------------------------------
## Fixture: same toy design as test_two_block_cpp.R
## J groups, random intercept + slope, level-2 covariate on the intercept
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

gamma_int <- c(1.0, 0.5)
gamma_slp <- 0.8
b_int <- as.numeric(X_int %*% gamma_int) + rnorm(J, sd = 0.4)
b_slp <- as.numeric(X_slp %*% gamma_slp) + rnorm(J, sd = 0.4)
eta <- b_int[as.integer(grp)] + b_slp[as.integer(grp)] * z1

prior_list_block2 <- list(
  `(Intercept)` = list(mu = c(0, 0), Sigma = diag(4, 2L), dispersion = 0.16),
  slope         = list(mu = 0, Sigma = diag(4, 1L), dispersion = 0.16)
)

## ---------------------------------------------------------------------------
## Dense brute-force reference
## ---------------------------------------------------------------------------
## Builds P11 (q x q), P22 (J p_re x J p_re), P12 (q x J p_re) explicitly and
## returns the sorted eigenvalues of A = P11^{-1/2} P12 P22^{-1} P21 P11^{-1/2}
## plus the dense blocks (for the block-swap check).
dense_rate <- function(x, grp, group_levels, x_hyper, P_b, w, V_inv) {
  J <- length(group_levels)
  p_re <- ncol(x)
  q_k <- vapply(x_hyper, ncol, integer(1L))
  q <- sum(q_k)
  cols <- split(seq_len(q), rep(seq_along(q_k), q_k))

  H <- lapply(seq_len(J), function(j) {
    Hj <- matrix(0, p_re, q)
    for (k in seq_len(p_re)) {
      X_k <- as.matrix(x_hyper[[k]])
      Hj[k, cols[[k]]] <- X_k[j, ]
    }
    Hj
  })

  P22 <- matrix(0, J * p_re, J * p_re)
  P12 <- matrix(0, q, J * p_re)
  P11 <- matrix(0, q, q)
  for (j in seq_len(J)) {
    rows <- which(as.integer(grp) == j)
    Z_j <- x[rows, , drop = FALSE]
    B_j <- crossprod(Z_j, Z_j * w[rows]) + P_b
    bc <- (j - 1L) * p_re + seq_len(p_re)
    P22[bc, bc] <- B_j
    P12[, bc] <- -t(H[[j]]) %*% P_b
    P11 <- P11 + t(H[[j]]) %*% P_b %*% H[[j]]
  }
  for (k in seq_along(cols)) {
    P11[cols[[k]], cols[[k]]] <- P11[cols[[k]], cols[[k]]] + V_inv[[k]]
  }

  inv_sqrt <- function(M) {
    e <- eigen(0.5 * (M + t(M)), symmetric = TRUE)
    e$vectors %*% diag(1 / sqrt(e$values), nrow(M)) %*% t(e$vectors)
  }

  P11_is <- inv_sqrt(P11)
  A <- P11_is %*% P12 %*% solve(P22, t(P12)) %*% P11_is
  ev <- sort(eigen(0.5 * (A + t(A)), symmetric = TRUE,
                   only.values = TRUE)$values, decreasing = TRUE)

  P22_is <- inv_sqrt(P22)
  A_swap <- P22_is %*% t(P12) %*% solve(P11, P12) %*% P22_is
  ev_swap <- sort(eigen(0.5 * (A_swap + t(A_swap)), symmetric = TRUE,
                        only.values = TRUE)$values, decreasing = TRUE)

  list(ev = ev, ev_swap = ev_swap, q = q)
}

V_inv_list <- lapply(prior_list_block2, function(pl) {
  chol2inv(chol(as.matrix(pl$Sigma)))
})

## ---------------------------------------------------------------------------
## 1. Gaussian, diagonal P_b: fast path == dense brute force
## ---------------------------------------------------------------------------
pl1_diag <- list(Sigma = diag(0.25, 2L), dispersion = 0.25)
P_b_diag <- chol2inv(chol(pl1_diag$Sigma))
w_gauss <- rep(1 / pl1_diag$dispersion, nrow(x_re))

rate1 <- two_block_rate(
  x = x_re, block = grp, x_hyper = x_hyper,
  prior_list_block1 = pl1_diag,
  prior_list_block2 = prior_list_block2,
  family = gaussian()
)
stopifnot(inherits(rate1, "two_block_rate"))
stopifnot(rate1$dims$q == 3L, rate1$dims$J == J, rate1$dims$p_re == p_re)
stopifnot(all(rate1$eigenvalues >= 0), all(rate1$eigenvalues < 1))
stopifnot(identical(rate1$weights_source, "dispersion"))

ref1 <- dense_rate(x_re, grp, group_levels, x_hyper, P_b_diag,
                   w_gauss, V_inv_list)
stopifnot(max(abs(rate1$eigenvalues - ref1$ev)) < 1e-10)
cat("1. gaussian, diagonal P_b: spectrum matches dense brute force (max diff ",
    format(max(abs(rate1$eigenvalues - ref1$ev)), digits = 3), ")\n", sep = "")

## Block-swap invariance: nonzero spectrum of the b-side matrix coincides
stopifnot(max(abs(ref1$ev_swap[seq_len(ref1$q)] - rate1$eigenvalues)) < 1e-10)
stopifnot(all(abs(ref1$ev_swap[-seq_len(ref1$q)]) < 1e-10))
cat("2. block-swap invariance: nonzero spectra coincide\n")

## ---------------------------------------------------------------------------
## 3. Gaussian, non-diagonal P_b (correlated REs): general path vs dense
## ---------------------------------------------------------------------------
Sigma_corr <- matrix(c(0.25, 0.10, 0.10, 0.25), 2L)
pl1_corr <- list(Sigma = Sigma_corr, dispersion = 0.25)
P_b_corr <- chol2inv(chol(Sigma_corr))

rate2 <- two_block_rate(
  x = x_re, block = grp, x_hyper = x_hyper,
  prior_list_block1 = pl1_corr,
  prior_list_block2 = prior_list_block2,
  family = gaussian()
)
ref2 <- dense_rate(x_re, grp, group_levels, x_hyper, P_b_corr,
                   w_gauss, V_inv_list)
stopifnot(max(abs(rate2$eigenvalues - ref2$ev)) < 1e-10)
cat("3. gaussian, correlated P_b: spectrum matches dense brute force\n")

## ---------------------------------------------------------------------------
## 4. Weights path (Poisson IRLS heuristic): runs, spectrum in [0, 1)
## ---------------------------------------------------------------------------
w_pois <- exp(0.3 * eta)   # IRLS weights at the "true" linear predictor
rate3 <- two_block_rate(
  x = x_re, block = grp, x_hyper = x_hyper,
  prior_list_block1 = list(Sigma = diag(0.25, 2L)),
  prior_list_block2 = prior_list_block2,
  weights = w_pois,
  family = poisson()
)
stopifnot(all(rate3$eigenvalues >= 0), all(rate3$eigenvalues < 1))
stopifnot(identical(rate3$weights_source, "user"))
ref3 <- dense_rate(x_re, grp, group_levels, x_hyper, P_b_diag,
                   w_pois, V_inv_list)
stopifnot(max(abs(rate3$eigenvalues - ref3$ev)) < 1e-10)
cat("4. weights path (poisson IRLS heuristic): matches dense brute force\n")

## Non-gaussian without weights must error
err <- tryCatch(
  two_block_rate(
    x = x_re, block = grp, x_hyper = x_hyper,
    prior_list_block1 = list(Sigma = diag(0.25, 2L)),
    prior_list_block2 = prior_list_block2,
    family = poisson()
  ),
  error = function(e) conditionMessage(e)
)
stopifnot(is.character(err), grepl("weights", err))

## ---------------------------------------------------------------------------
## 5. m_for_tol and print method
## ---------------------------------------------------------------------------
lam <- rate1$lambda_star
m3 <- rate1$m_for_tol(1e-3)
stopifnot(lam^m3 <= 1e-3, lam^(m3 - 1L) > 1e-3)
print(rate1)
cat("5. m_for_tol consistent with lambda*^m\n")

## ---------------------------------------------------------------------------
## 6. erf_n closed form: pchisq identity vs Remark 1 explicit series
## ---------------------------------------------------------------------------
## erf_{2m}(x)   = 1 - exp(-x^2) * sum_{j=0}^{m-1} x^{2j}/j!
## erf_{2m+1}(x) = erf_1(x) - exp(-x^2)/sqrt(pi) *
##                 sum_{j=1}^{m} (2x)^{2j-1} (j-1)! / (2j-1)!
erfn_series <- function(x, n) {
  if (n %% 2L == 0L) {
    m <- n %/% 2L
    1 - exp(-x^2) * sum(x^(2 * (0:(m - 1))) / factorial(0:(m - 1)))
  } else {
    m <- (n - 1L) %/% 2L
    base <- 2 * stats::pnorm(sqrt(2) * x) - 1     # classical erf(x)
    if (m == 0L) return(base)
    j <- seq_len(m)
    base - exp(-x^2) / sqrt(pi) *
      sum((2 * x)^(2 * j - 1) * factorial(j - 1) / factorial(2 * j - 1))
  }
}
for (n in 1:6) {
  for (x in c(0.1, 0.7, 1.5, 3)) {
    stopifnot(abs(glmbayesCore:::.two_block_erfn(x, n) -
                  erfn_series(x, n)) < 1e-12)
  }
}
cat("6. erf_n pchisq identity matches Remark 1 series (n = 1..6)\n")

## ---------------------------------------------------------------------------
## 7. TV bounds: Theorem 3 exact vs Corollary 1 envelope
## ---------------------------------------------------------------------------
l_grid <- c(1L, 2L, 3L, 5L, 8L, 12L, 20L)
b_t3 <- two_block_tv_bound(rate1, l_grid, method = "theorem3")
b_c1 <- two_block_tv_bound(rate1, l_grid, method = "corollary1")

stopifnot(all(b_t3 >= 0), all(b_t3 <= 1), all(b_c1 >= 0), all(b_c1 <= 1))
## Theorem 3 terms are exact; Corollary 1 is a term-wise relaxation
stopifnot(all(b_t3 <= b_c1 + 1e-12))
## both monotone decreasing in l, and -> 0
stopifnot(all(diff(b_t3) <= 1e-12), all(diff(b_c1) <= 1e-12))
stopifnot(two_block_tv_bound(rate1, 60L) < 1e-12)
cat("7. TV bounds: theorem3 <= corollary1, monotone decreasing, -> 0\n")

cat("\n   l    theorem3      corollary1\n")
for (ii in seq_along(l_grid)) {
  cat(sprintf("  %3d   %.6e  %.6e\n", l_grid[ii], b_t3[ii], b_c1[ii]))
}

## Mean term: D0 > 0 increases the bound; D0 = 0 default unchanged
b_d0 <- two_block_tv_bound(rate1, l_grid, method = "theorem3", D0 = 4)
stopifnot(all(b_d0 >= b_t3))
stopifnot(identical(two_block_tv_bound(rate1, 5L),
                    two_block_tv_bound(rate1, 5L, D0 = 0)))

## ---------------------------------------------------------------------------
## 8. l_for_tv: inversion consistency for both methods
## ---------------------------------------------------------------------------
for (mth in c("theorem3", "corollary1")) {
  for (tol in c(1e-2, 1e-3, 1e-6)) {
    l_star <- two_block_l_for_tv(rate1, tol, method = mth)
    stopifnot(two_block_tv_bound(rate1, l_star, method = mth) <= tol)
    if (l_star > 1L) {
      stopifnot(two_block_tv_bound(rate1, l_star - 1L, method = mth) > tol)
    }
  }
}
cat("8. l_for_tv inverts the bound exactly (both methods)\n")

cat("\ntest_two_block_rate.R: all checks passed\n")
