## Stage 2 check for the switch of the binomial likelihoods to dbinom_raw():
## compare glmbfamfunc() f1/f2/f4 against the previous stats::dbinom formulas,
## then fit glmb() models for each link.
## Run after reinstalling glmbayesCore:
##   Rscript data-raw/check_binomial_dbinom_raw_switch.R

library(glmbayesCore)

## Previous R formulas (before the switch), for reference
old_f1 <- list(
  logit = function(b, y, x, alpha = 0, wt = 1) {
    lambda <- t(exp(alpha + x %*% b)); p <- lambda / (1 + lambda)
    -sum(dbinom(round(wt * y), round(wt), p, log = TRUE))
  },
  probit = function(b, y, x, alpha = 0, wt = 1) {
    p <- pnorm(alpha + x %*% b)
    -sum(dbinom(round(wt * y), round(wt), p, log = TRUE))
  },
  cloglog = function(b, y, x, alpha = 0, wt = 1) {
    p <- 1 - exp(-exp(alpha + x %*% b))
    -sum(dbinom(round(wt * y), round(wt), p, log = TRUE))
  }
)
old_f4_sat <- function(y, wt) -sum(dbinom(round(wt * y), round(wt), y, log = TRUE))

## Exact negative log-likelihood with log p and log q computed on the log scale
log_pq <- list(
  logit   = function(eta) cbind(plogis(eta, log.p = TRUE), plogis(-eta, log.p = TRUE)),
  probit  = function(eta) cbind(pnorm(eta, log.p = TRUE), pnorm(-eta, log.p = TRUE)),
  cloglog = function(eta) cbind(log(-expm1(-exp(eta))), -exp(eta))
)
exact_f1 <- function(lnk, b, y, x, wt) {
  s <- round(wt * y); n <- round(wt)
  lpq <- log_pq[[lnk]](as.vector(x %*% b))
  ll <- lchoose(n, s) + ifelse(s > 0, s * lpq[, 1], 0) + ifelse(n > s, (n - s) * lpq[, 2], 0)
  -sum(ll)
}

## 1. New closures vs the exact reference (must agree); old formulas for information
##    (old cloglog loses precision once exp(-exp(eta)) < ~1e-12, and is Inf beyond eta ~ 3.6)
n_obs <- 40
x <- cbind(1, rnorm(n_obs), runif(n_obs))
wt <- sample(1:12, n_obs, replace = TRUE)
y <- rbinom(n_obs, wt, 0.4) / wt
mu <- matrix(0, 3, 1); P <- diag(0.1, 3)

for (lnk in names(old_f1)) {
  ff <- glmbfamfunc(binomial(link = lnk))
  res <- replicate(200, {
    b <- rnorm(3, sd = 0.5)
    ex <- exact_f1(lnk, b, y, x, wt)
    new1 <- ff$f1(b, y, x, wt = wt)
    c(rel_f1 = abs(new1 - ex) / max(1, abs(ex)),
      rel_f2 = abs(ff$f2(b, y, x, mu, P, wt = wt) -
                     (ex + 0.5 * t(b - mu) %*% P %*% (b - mu))) / max(1, abs(ex)),
      rel_f4 = abs(ff$f4(b, y, x, wt = wt) - (2 * ex - 2 * old_f4_sat(y, wt))) /
        max(1, abs(ex)),
      old_diff = abs(old_f1[[lnk]](b, y, x, wt = wt) - new1))
  })
  rel <- max(res[c("rel_f1", "rel_f2", "rel_f4"), ])
  message(sprintf("1. %-7s max relative |new - exact| = %.3g | old vs new: median %.3g, max %.3g",
                  lnk, rel, median(res["old_diff", ]), max(res["old_diff", ])))
  stopifnot(rel < 1e-10)
}

## 2. Extreme linear predictors: new closures stay finite where p or q -> 0
x1 <- matrix(1, 3, 1); y1 <- c(0, 0.5, 1); wt1 <- c(2, 2, 2)
for (lnk in names(old_f1)) {
  ff <- glmbfamfunc(binomial(link = lnk))
  for (b in c(-40, 40, 800)) {
    message(sprintf("2. %-7s b = %4g: new f1 = %-12.6g old f1 = %.6g", lnk, b,
                    ff$f1(b, y1, x1, wt = wt1),
                    suppressWarnings(old_f1[[lnk]](b, y1, x1, wt = wt1))))
  }
}

## 3. Model fits (optim mode, envelope, accept/reject) with a diffuse prior:
##    posterior mode should be close to glm(); compare with a run before the switch
colnames(x) <- c("(Intercept)", "x1", "x2")
for (lnk in names(old_f1)) {
  fam <- binomial(link = lnk)
  fit <- rglmb(n = 2000, y = y, x = x, family = fam,
               pfamily = dNormal(mu = rep(0, 3), Sigma = diag(100, 3)),
               weights = wt)
  ref <- glm.fit(x, y, weights = wt, family = fam)
  message(sprintf("3. %-7s coef.mode - glm coef: %s | mean iters/draw = %.3f",
                  lnk,
                  paste(signif(fit$coef.mode - ref$coefficients, 3), collapse = ", "),
                  mean(fit$iters)))
  print(summary(fit))
}

message("check_binomial_dbinom_raw_switch.R: done")
