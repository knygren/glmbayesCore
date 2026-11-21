## Validate glmbayesCore::dbinom_raw() against stats::dbinom() on integer
## inputs, and check its behavior on non-integer and invalid inputs.
## Run after Rcpp::compileAttributes() and reinstalling glmbayesCore:
##   Rscript data-raw/check_dbinom_raw.R

library(glmbayesCore)

rel_diff <- function(a, b) {
  same <- (a == b) | (is.na(a) & is.na(b))
  d <- abs(a - b) / pmax(1, abs(b))
  d[same] <- 0
  d
}

## 1. Integer grid: must match stats::dbinom (density and log density)
ns <- c(0, 1, 2, 5, 50, 1e3, 1e6)
ps <- c(0, 1e-12, 0.01, 0.5, 0.99, 1 - 1e-9, 1)
grid <- do.call(rbind, lapply(ns, function(n) {
  xs <- unique(c(-1, 0, 1, n - 1, floor(n / 2), n, n + 1))
  expand.grid(x = xs, n = n, p = ps)
}))

max_diff <- 0
for (lg in c(FALSE, TRUE)) {
  ref <- stats::dbinom(grid$x, grid$n, grid$p, log = lg)
  new <- dbinom_raw(grid$x, grid$n, grid$p, log = lg)
  d <- rel_diff(new, ref)
  max_diff <- max(max_diff, d)
  if (any(d > 1e-12)) {
    bad <- cbind(grid, ref = ref, new = new)[d > 1e-12, ]
    print(utils::head(bad, 20))
    stop("dbinom_raw differs from stats::dbinom (log = ", lg, ")")
  }
}
message("1. integer grid: ", nrow(grid), " points x 2, max rel diff = ",
        signif(max_diff, 3))

## 2. Directly computed q near p = 1
eta <- c(10, 20, 30, 35)
d_direct <- dbinom_raw(3, 5, plogis(eta), plogis(-eta), log = TRUE)
d_default <- dbinom_raw(3, 5, plogis(eta), log = TRUE)
message("2. q near p = 1 (direct vs 1 - p):")
print(data.frame(eta = eta, direct_q = d_direct, one_minus_p = d_default))

## 3. Non-integer x and n: finite, and continuous between integers
n0 <- 10
x_fine <- seq(0, n0, by = 0.01)
lp <- dbinom_raw(x_fine, n0, 0.3, log = TRUE)
stopifnot(all(is.finite(lp)))
at_int <- x_fine %in% 0:n0
stopifnot(max(rel_diff(lp[at_int], stats::dbinom(0:n0, n0, 0.3, log = TRUE))) < 1e-12)
stopifnot(max(abs(diff(lp))) < 0.5)
lp_n <- dbinom_raw(3, seq(3, 6, by = 0.01), 0.3, log = TRUE)
stopifnot(all(is.finite(lp_n)), max(abs(diff(lp_n))) < 0.1)
message("3. non-integer x and n: finite and smooth")

## 4. Non-integer values agree with the gamma-function form
##    lgamma(n+1) - lgamma(x+1) - lgamma(n-x+1) + x log p + (n-x) log q
xs <- c(0.3, 2.25, 4.5, 7.75); nn <- 9.5; p <- 0.35
lg_form <- lgamma(nn + 1) - lgamma(xs + 1) - lgamma(nn - xs + 1) +
  xs * log(p) + (nn - xs) * log1p(-p)
stopifnot(max(rel_diff(dbinom_raw(xs, nn, p, log = TRUE), lg_form)) < 1e-10)
message("4. non-integer values match the lgamma form")

## 5. Invalid inputs: NaN with a warning, like stats::dbinom
expect_nan_warning <- function(expr) {
  w <- NULL
  val <- withCallingHandlers(expr, warning = function(cnd) {
    w <<- conditionMessage(cnd)
    invokeRestart("muffleWarning")
  })
  stopifnot(is.nan(val), identical(w, "NaNs produced"))
}
expect_nan_warning(dbinom_raw(1, -1, 0.5))
expect_nan_warning(dbinom_raw(1, 2, 1.5))
expect_nan_warning(dbinom_raw(1, 2, -0.1))
expect_nan_warning(dbinom_raw(1, 2, 0.5, q = 0.4))
stopifnot(is.na(dbinom_raw(NA, 2, 0.5)), is.nan(suppressWarnings(stats::dbinom(1, -1, 0.5))))
stopifnot(dbinom_raw(-1, 2, 0.5) == 0, dbinom_raw(3, 2, 0.5) == 0,
          dbinom_raw(3, 2, 0.5, log = TRUE) == -Inf)
stopifnot(identical(dbinom_raw(numeric(0), 2, 0.5), numeric(0)))
message("5. invalid inputs: NaN + warning; out-of-range x gives 0")

message("check_dbinom_raw.R: OK")
