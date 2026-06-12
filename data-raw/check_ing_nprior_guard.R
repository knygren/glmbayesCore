## Check the n_prior <= n_w guard in rindepNormalGamma_reg.
## Run: Rscript data-raw/check_ing_nprior_guard.R
pkgload::load_all("C:/Rpackages/glmbayesCore", export_all = FALSE, quiet = TRUE)

set.seed(7)
n  <- 20L
x  <- cbind(1, rnorm(n))
y  <- as.numeric(x %*% c(2, 0.5) + rnorm(n))
mu <- c(0, 0)
Sigma <- diag(100, 2L)

## 1. Violating prior: n_prior = 1000 >> n_w = 20 (shape = (n_prior+1+p)/2)
shape_bad <- (1000 + 1 + 2) / 2
pf_bad <- dIndependent_Normal_Gamma(mu = mu, Sigma = Sigma,
                                    shape = shape_bad, rate = shape_bad)
err <- tryCatch(
  rlmb(n = 10, y = y, x = x, pfamily = pf_bad, progbar = FALSE),
  error = function(e) conditionMessage(e)
)
stopifnot(is.character(err), grepl("n_prior <= n_w", err, fixed = TRUE))
cat("1. violating prior (n_prior = 1000, n_w = 20) rejected: OK\n")

## 2. Boundary-compliant prior: n_prior = n_w = 20 must pass the guard
shape_ok <- (20 + 1 + 2) / 2
pf_ok <- dIndependent_Normal_Gamma(mu = mu, Sigma = Sigma,
                                   shape = shape_ok, rate = shape_ok)
fit <- rlmb(n = 25, y = y, x = x, pfamily = pf_ok, progbar = FALSE)
stopifnot(is.matrix(fit$coefficients), nrow(fit$coefficients) == 25L)
cat("2. compliant prior (n_prior = n_w = 20) samples: OK\n")

cat("\nAll n_prior guard checks passed.\n")
