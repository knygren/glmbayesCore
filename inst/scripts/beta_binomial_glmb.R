#!/usr/bin/env Rscript

## Beta–Binomial conjugate: Beta prior on the probability θ of a Binomial likelihood
## (intercept-only, identity link on θ).
##
## Model:
##   Y_i | θ ~ Binomial(n_i, θ),  identity link → coefficient = θ directly
##   θ ~ Beta(α, β)                                              (prior)
##   θ | y ~ Beta(α + total_successes, β + total_failures)      (posterior)
##
## This is the same update as the Beta-Binomial scalar example in Chapter 02,
## but here we fit it via glmb() + dBeta() and verify the draws match the
## analytic posterior exactly.

library(glmbayes)

## =============================================================================
## Example A — Bayes Rules! style: binary outcomes, Beta(2,2) prior
## =============================================================================

cat("\n=== Example A: binary outcomes, Beta(2, 2) prior ===\n\n")

set.seed(42)
n_obs   <- 25L
theta_true <- 0.28
## Individual binary outcomes; weights = 1 per row makes the trial count explicit.
y_A     <- rbinom(n_obs, size = 1, prob = theta_true)

alpha0_A <- 2;  beta0_A <- 2    ## Beta(2, 2) prior: mean = 0.5

## Analytic conjugate posterior:
##   shape1_post = alpha0 + sum(y) = 2 + sum(y_A)
##   shape2_post = beta0  + (n - sum(y)) = 2 + (25 - sum(y_A))
shape1_post_A <- alpha0_A + sum(y_A)
shape2_post_A <- beta0_A  + (n_obs - sum(y_A))

cat(sprintf("Observed data: n = %d, successes = %d, p_hat = %.4f\n",
            n_obs, sum(y_A), mean(y_A)))
cat(sprintf("Prior:    Beta(shape1 = %g, shape2 = %g)  -->  prior mean = %.4f\n",
            alpha0_A, beta0_A, alpha0_A / (alpha0_A + beta0_A)))
cat(sprintf("Posterior: Beta(shape1 = %.1f, shape2 = %.1f)\n",
            shape1_post_A, shape2_post_A))
cat(sprintf("  Posterior mean theta = %.6f\n",
            shape1_post_A / (shape1_post_A + shape2_post_A)))
cat(sprintf("  90%% CI for theta: [%.4f, %.4f]\n\n",
            qbeta(0.05, shape1_post_A, shape2_post_A),
            qbeta(0.95, shape1_post_A, shape2_post_A)))

## Fit with glmb()
beta_init_A <- matrix(alpha0_A / (alpha0_A + beta0_A), nrow = 1L, ncol = 1L)
colnames(beta_init_A) <- "(Intercept)"

pf_A <- dBeta(shape1 = alpha0_A, shape2 = beta0_A, beta = beta_init_A)

data_A <- data.frame(y = y_A)
set.seed(2026)
fit_A <- glmb(
  n       = 20000,
  y ~ 1,
  data    = data_A,
  weights = rep(1L, n_obs),    ## explicit: each row is one Bernoulli trial
  family  = binomial(link = "identity"),
  pfamily = pf_A
)

cat("glmb() summary (coefficient = Binomial probability theta):\n")
print(summary(fit_A))

## Verify: draws should match the analytic posterior
smp_A <- fit_A$coefficients[, 1L]
cat(sprintf("glmb draw mean = %.6f  |  analytic mean = %.6f\n",
            mean(smp_A),
            shape1_post_A / (shape1_post_A + shape2_post_A)))
cat(sprintf("glmb draw SD   = %.6f  |  analytic SD   = %.6f\n",
            sd(smp_A),
            sqrt(shape1_post_A * shape2_post_A /
                   ((shape1_post_A + shape2_post_A)^2 *
                      (shape1_post_A + shape2_post_A + 1)))))


## =============================================================================
## Example B — Prior_Setup() calibration for the same data
## =============================================================================

cat("\n=== Example B: same data via Prior_Setup() ===\n\n")

ps_A <- Prior_Setup(
  y ~ 1,
  data    = data_A,
  weights = rep(1L, n_obs),
  family  = binomial(link = "identity"),
  pwt     = 0.01
)
print(ps_A)

cb    <- ps_A$conj_binomial
pf_ps <- dBeta(shape1 = cb$shape1, shape2 = cb$shape2, beta = cb$beta)

set.seed(2026)
fit_ps <- glmb(
  n       = 20000,
  y ~ 1,
  data    = data_A,
  weights = rep(1L, n_obs),
  family  = binomial(link = "identity"),
  pfamily = pf_ps
)

cat("glmb() summary (Prior_Setup calibrated Beta prior, pwt = 0.01):\n")
print(summary(fit_ps))

message("\nFinished Beta-Binomial conjugate demos.")
