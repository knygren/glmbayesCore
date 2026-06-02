## Gamma-Poisson Conjugacy: dGamma(Inv_Dispersion = FALSE) with glmb()
##
## Mirrors vignette Chapter-02-S04 (Gamma-Poisson Conjugacy for One Count Rate).
## Part 1 requires the 'bayesrules' package (in Suggests).
## Appendix A requires the 'LearnBayes' package (in Suggests).
##
## Two examples:
##   Part 1      - Daily bicycle rental counts (Bayes Rules! scenario)
##                 Visualise prior/likelihood/posterior with plot_gamma_poisson().
##                 Fit with glmb() using dGamma(Inv_Dispersion = FALSE).
##   Appendix A  - Heart transplant mortality (Albert 2009 / LearnBayes)
##                 Exposure-weighted Poisson; fit two hospitals with glmb().
##
## Run: demo("Ex_13_GammaPoisson", package = "glmbayes")
## See also: vignette("Chapter-02-S04", package = "glmbayes")

library(glmbayes)

if (!requireNamespace("bayesrules", quietly = TRUE))
  stop("Part 1 of this demo requires the 'bayesrules' package. ",
       "Install it with: install.packages('bayesrules')")

library(bayesrules)

## ---- Part 1: Daily bicycle counts (Bayes Rules! scenario) --------------------

cat("\n=== Part 1: Daily bike counts (Gamma(3,4) prior, n=9, sum(y)=3) ===\n\n")

## Data: 9 days, first three had 1 rental each, remaining six had none
br_y_df  <- data.frame(y = c(1, 1, 1, rep(0, 6)))
br_n     <- nrow(br_y_df)
br_shape <- 3
br_rate  <- 4

## Analytic posterior: Gamma(alpha + sum(y), beta + n)
post_shape_br <- br_shape + sum(br_y_df$y)   ## 3 + 3 = 6
post_rate_br  <- br_rate  + br_n             ## 4 + 9 = 13

cat(sprintf("Prior:     Gamma(%d, %d)  =>  mean = %.4f\n", br_shape, br_rate,
            br_shape / br_rate))
cat(sprintf("Data:      n = %d, sum(y) = %d  =>  mean(y) = %.4f\n",
            br_n, sum(br_y_df$y), mean(br_y_df$y)))
cat(sprintf("Posterior: Gamma(%d, %d)  =>  mean = %.4f,  SD = %.4f\n\n",
            post_shape_br, post_rate_br,
            post_shape_br / post_rate_br,
            sqrt(post_shape_br) / post_rate_br))

## Overlay prior, likelihood (scaled), and posterior densities
plot_gamma_poisson(
  shape      = br_shape,
  rate       = br_rate,
  sum_y      = sum(br_y_df$y),
  n          = br_n,
  prior      = TRUE,
  likelihood = TRUE,
  posterior  = TRUE
)

readline("press any key to continue")

## Tabular summary of prior, likelihood, and posterior
summarize_gamma_poisson(shape = br_shape, rate = br_rate,
                        sum_y = sum(br_y_df$y), n = br_n)

readline("press any key to continue")

## ---- Fitting with glmb() using dGamma(Inv_Dispersion = FALSE) ----------------

cat("\n--- Fitting bike-counts with glmb() + dGamma(Inv_Dispersion = FALSE) ---\n\n")
cat("family = poisson(link = 'identity'),  pfamily = dGamma(Inv_Dispersion = FALSE)\n\n")

gp_beta <- matrix(br_shape / br_rate, 1L, 1L, dimnames = list(NULL, "(Intercept)"))
gp_pf   <- dGamma(shape = br_shape, rate = br_rate,
                  beta = gp_beta, Inv_Dispersion = FALSE)

set.seed(2026)
fit_gp <- glmb(
  n       = 20000,
  y ~ 1,
  data    = br_y_df,
  weights = rep(1L, br_n),
  family  = poisson(link = "identity"),
  pfamily = gp_pf
)
print(fit_gp)

readline("press any key to continue")

## ---- Analytic vs. glmb() comparison ------------------------------------------

glmb_mean_gp <- fit_gp$coef.means["(Intercept)"]
glmb_sd_gp   <- sd(fit_gp$coefficients[, "(Intercept)", drop = TRUE])

cat("\nAnalytic vs. glmb() posterior (daily bike counts):\n")
cmp_gp <- data.frame(
  Posterior      = sprintf("Gamma(%d, %d)", post_shape_br, post_rate_br),
  Analytic.Mean  = round(post_shape_br / post_rate_br, 5),
  Analytic.SD    = round(sqrt(post_shape_br) / post_rate_br, 5),
  glmb.Post.Mean = round(glmb_mean_gp, 5),
  glmb.Post.Sd   = round(glmb_sd_gp,   5),
  check.names    = FALSE
)
print(cmp_gp, row.names = FALSE)
cat("\nglmb Post.Mean and Post.Sd should match the analytic values to Monte Carlo error.\n")
cat("dGamma(Inv_Dispersion = FALSE) implements the exact conjugate update.\n\n")

readline("press any key to continue")

## ---- Appendix A: Heart transplant mortality (Albert 2009 / LearnBayes) -------

if (!requireNamespace("LearnBayes", quietly = TRUE)) {
  cat("Appendix A skipped: 'LearnBayes' package not available.\n")
  cat("Install with: install.packages('LearnBayes')\n")
} else {
  library(LearnBayes)
  data("hearttransplants")

  cat("=== Appendix A: Heart transplant mortality (Albert 2009) ===\n\n")
  cat(sprintf("%d hospitals  |  sum(y) = %d  |  sum(e) = %.0f\n\n",
              nrow(hearttransplants), sum(hearttransplants$y),
              sum(hearttransplants$e)))

  ## Albert (Ch. 3.2): Gamma(16, 15174) prior on lambda
  ## Focus on two specific hospitals
  alpha0_ht <- 16L;  beta0_ht <- 15174L
  ex_A      <- 66L;  yobs_A   <- 1L
  ex_B      <- 1767L; yobs_B  <- 4L

  cat("Albert's prior: Gamma(16, 15174)  =>  mean ~= 0.001054 (national baseline)\n\n")

  cat("Analytic posteriors:\n")
  ht_analytic <- data.frame(
    Hospital = c("A", "B"),
    e        = c(ex_A, ex_B),
    y        = c(yobs_A, yobs_B),
    Posterior = c(
      sprintf("Gamma(%d, %d)", alpha0_ht + yobs_A, beta0_ht + ex_A),
      sprintf("Gamma(%d, %d)", alpha0_ht + yobs_B, beta0_ht + ex_B)
    ),
    Mean = c(
      (alpha0_ht + yobs_A) / (beta0_ht + ex_A),
      (alpha0_ht + yobs_B) / (beta0_ht + ex_B)
    ),
    SD = c(
      sqrt(alpha0_ht + yobs_A) / (beta0_ht + ex_A),
      sqrt(alpha0_ht + yobs_B) / (beta0_ht + ex_B)
    ),
    check.names = FALSE
  )
  print(ht_analytic, digits = 6, row.names = FALSE)

  readline("press any key to continue")

  ## Fit each hospital separately; exposure enters as weights
  ht_beta <- matrix(alpha0_ht / beta0_ht, 1L, 1L,
                    dimnames = list(NULL, "(Intercept)"))
  ht_pf   <- dGamma(shape = alpha0_ht, rate = beta0_ht,
                    beta = ht_beta, Inv_Dispersion = FALSE)

  cat("\n--- Hospital A (e = 66, y = 1) ---\n")
  set.seed(2026)
  fit_A <- glmb(n = 20000, y ~ 1, data = data.frame(y = yobs_A),
                weights = ex_A, family = poisson(link = "identity"),
                pfamily = ht_pf)
  print(fit_A)

  readline("press any key to continue")

  cat("\n--- Hospital B (e = 1767, y = 4) ---\n")
  set.seed(2026)
  fit_B <- glmb(n = 20000, y ~ 1, data = data.frame(y = yobs_B),
                weights = ex_B, family = poisson(link = "identity"),
                pfamily = ht_pf)
  print(fit_B)

  readline("press any key to continue")

  cat("\nAlbert analytic vs. glmb() comparison:\n")
  ht_compare <- data.frame(
    Hospital       = c("A", "B"),
    Posterior      = ht_analytic$Posterior,
    Albert.Mean    = round(ht_analytic$Mean, 6),
    Albert.SD      = round(ht_analytic$SD,   6),
    glmb.Post.Mean = round(c(fit_A$coef.means["(Intercept)"],
                              fit_B$coef.means["(Intercept)"]), 6),
    glmb.Post.Sd   = round(c(sd(fit_A$coefficients[, "(Intercept)", drop = TRUE]),
                              sd(fit_B$coefficients[, "(Intercept)", drop = TRUE])), 6),
    check.names = FALSE
  )
  print(ht_compare, row.names = FALSE)
  cat("\nEach fit is for one hospital; (Intercept) = posterior draw for that hospital's rate.\n")
}

cat("\nSee vignette('Chapter-02-S04', package = 'glmbayes') for full derivations.\n")
cat("For Poisson regression with covariates use family = poisson(log) + dNormal().\n")
