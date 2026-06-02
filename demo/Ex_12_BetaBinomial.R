## Beta-Binomial Conjugacy: dBeta() prior with glmb()
##
## Mirrors vignette Chapter-02-S03 (Beta-Binomial Conjugacy for One Proportion).
## Requires the 'bayesrules' package (in Suggests).
##
## Two examples:
##   Part 1  - Analgesic trial scenario from Bayes Rules!
##             Visualise prior/likelihood/posterior with plot_beta_binomial().
##   Part 2  - Bechdel test (bayesrules::bechdel, n ~ 1794 films)
##             Fit with glmb() using dBeta(); compare to analytic posterior.
##
## Run: demo("Ex_12_BetaBinomial", package = "glmbayes")
## See also: vignette("Chapter-02-S03", package = "glmbayes")

library(glmbayes)

if (!requireNamespace("bayesrules", quietly = TRUE))
  stop("This demo requires the 'bayesrules' package. ",
       "Install it with: install.packages('bayesrules')")

library(bayesrules)

## ---- Part 1: Analgesic trial (Bayes Rules! scenario) -------------------------

cat("\n=== Part 1: Analgesic trial (Beta(4,6) prior, n=30, y=14) ===\n\n")

## Prior: Beta(4, 6)  =>  mean = 0.4, effective prior sample size = 10
## Data: 14 responders out of 30 patients
a0 <- 4;  b0 <- 6
y  <- 14; n  <- 30

## Analytic posterior: Beta(a0 + y, b0 + n - y)
post_a <- a0 + y          ## 4 + 14 = 18
post_b <- b0 + (n - y)    ## 6 + 16 = 22

cat(sprintf("Prior:     Beta(%d, %d)  =>  mean = %.3f\n", a0, b0,
            a0 / (a0 + b0)))
cat(sprintf("Data:      y = %d, n = %d  =>  freq = %.3f\n", y, n, y / n))
cat(sprintf("Posterior: Beta(%d, %d)  =>  mean = %.4f,  SD = %.4f\n\n",
            post_a, post_b,
            post_a / (post_a + post_b),
            sqrt(post_a * post_b / ((post_a + post_b)^2 * (post_a + post_b + 1)))))

## Overlay prior, likelihood (scaled), and posterior densities
plot_beta_binomial(
  alpha      = a0,
  beta       = b0,
  y          = y,
  n          = n,
  prior      = TRUE,
  likelihood = TRUE,
  posterior  = TRUE
)

readline("press any key to continue")

## Tabular summary of prior, likelihood, and posterior
summarize_beta_binomial(alpha = a0, beta = b0, y = y, n = n)

readline("press any key to continue")

## ---- Part 2: Bechdel test (bayesrules::bechdel) ------------------------------

cat("\n=== Part 2: Bechdel test (Beta(9,11) prior, n ~ 1794 films) ===\n\n")

## Binary outcome: 1 = PASS, 0 = FAIL
pass   <- as.integer(bechdel[["binary"]] == "PASS")
n_bech <- length(pass)
y_bech <- sum(pass)

cat(sprintf("n = %d,  y (PASS) = %d,  proportion = %.3f\n\n",
            n_bech, y_bech, y_bech / n_bech))

## Prior: Beta(9, 11)  =>  mean = 0.45, effective prior sample size = 20
a0_b <- 9;  b0_b <- 11

## Analytic posterior
post_a_b <- a0_b + y_bech
post_b_b <- b0_b + (n_bech - y_bech)

cat(sprintf("Prior:     Beta(%d, %d)  =>  mean = %.3f\n", a0_b, b0_b,
            a0_b / (a0_b + b0_b)))
cat(sprintf("Posterior: Beta(%d, %d)  =>  mean = %.4f,  SD = %.4f\n",
            post_a_b, post_b_b,
            post_a_b / (post_a_b + post_b_b),
            sqrt(post_a_b * post_b_b /
                   ((post_a_b + post_b_b)^2 * (post_a_b + post_b_b + 1)))))
cat("95% credible interval: ")
print(round(qbeta(c(0.025, 0.975), post_a_b, post_b_b), 4))
cat("\n")

## Overlay prior, likelihood (scaled), and posterior densities
plot_beta_binomial(
  alpha      = a0_b,
  beta       = b0_b,
  y          = y_bech,
  n          = n_bech,
  prior      = TRUE,
  likelihood = TRUE,
  posterior  = TRUE
)

readline("press any key to continue")

## Tabular summary
summarize_beta_binomial(alpha = a0_b, beta = b0_b, y = y_bech, n = n_bech)

readline("press any key to continue")

## ---- Fitting with glmb() using dBeta() ---------------------------------------

cat("\n--- Fitting the Bechdel test with glmb() + dBeta() ---\n\n")
cat("family = binomial(link = 'identity'),  pfamily = dBeta()\n\n")

df_bech  <- data.frame(y = pass)
bech_beta <- matrix(a0_b / (a0_b + b0_b), nrow = 1L, ncol = 1L,
                    dimnames = list(NULL, "(Intercept)"))
bech_pf  <- dBeta(shape1 = a0_b, shape2 = b0_b, beta = bech_beta)

set.seed(2026)
fit_bech <- glmb(
  n       = 20000,
  y ~ 1,
  data    = df_bech,
  weights = rep(1L, n_bech),
  family  = binomial(link = "identity"),
  pfamily = bech_pf
)
print(fit_bech)

readline("press any key to continue")

## ---- Analytic vs. glmb() comparison ------------------------------------------

analytic_mean_b <- post_a_b / (post_a_b + post_b_b)
analytic_sd_b   <- sqrt(post_a_b * post_b_b /
                          ((post_a_b + post_b_b)^2 * (post_a_b + post_b_b + 1)))

glmb_mean <- fit_bech$coef.means["(Intercept)"]
glmb_sd   <- sd(fit_bech$coefficients[, "(Intercept)", drop = TRUE])

cat("\nAnalytic vs. glmb() posterior (Bechdel test):\n")
cmp <- data.frame(
  Posterior      = sprintf("Beta(%d, %d)", post_a_b, post_b_b),
  Analytic.Mean  = round(analytic_mean_b, 5),
  Analytic.SD    = round(analytic_sd_b,   5),
  glmb.Post.Mean = round(glmb_mean, 5),
  glmb.Post.Sd   = round(glmb_sd,   5),
  check.names    = FALSE
)
print(cmp, row.names = FALSE)
cat("\nglmb Post.Mean and Post.Sd should match the analytic values to Monte Carlo error.\n")
cat("dBeta() implements the exact conjugate update.\n\n")

cat("See vignette('Chapter-02-S03', package = 'glmbayes') for full derivations.\n")
cat("For Binomial regression with covariates use family = binomial(logit) + dNormal().\n")
