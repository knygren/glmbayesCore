#!/usr/bin/env Rscript

## Gamma–Poisson conjugate (Gamma prior on lambda, Poisson likelihood, intercept-only identity link):
##
## **Part 1** — Matches the Gamma–Poisson illustration in Bayes Rules! vignette:
## `vignette("conjugate-families", package = "bayesrules")` (`plot_gamma_poisson(shape = 3, rate = 4,
## sum_y = 3, n = 9)`).  Any nonnegative integer outcome vector `y` of length 9 summing to 3 induces
## the same conjugate update.
##
## **Part 2** — Uses real daily rental counts `bayesrules::bikes$rides` (first week) together with the
## same `bayesrules::plot_gamma_poisson()` and `glmbayes::dGamma_Conjugate()` + `glmbayes::rglmb()`
## sampler as a second illustration (still intercept-only Lambda; pedagogical shorthand only).
##
## Suggested packages: `bayesrules`, `ggplot2` (already listed under Suggests in `glmbayes`).
## Run after installing the package, e.g.:
##    Rscript inst/scripts/bayesrules_gamma_poisson_glmb.R

require_suggested <- function(pkg, install_tip) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      sprintf("missing suggested package `%s` — %s", pkg, install_tip),
      call. = FALSE
    )
  }
}

require_suggested(
  "bayesrules",
  "install.packages(\"bayesrules\", repos = \"https://cloud.r-project.org\")"
)
require_suggested(
  "ggplot2",
  "install.packages(\"ggplot2\", repos = \"https://cloud.r-project.org\")"
)

if (!("package:glmbayes" %in% search())) {
  suppressPackageStartupMessages(library(glmbayes))
}
suppressPackageStartupMessages({
  library(bayesrules)
  library(ggplot2)
})

## --- Example A (vignette numbers) -----------------------------------------------

shape_a <- 3
rate_a <- 4
n_a <- 9L
y_a <- integer(n_a)
y_a[seq_len(3L)] <- 1L
stopifnot(sum(y_a) == 3)

print(summarize_gamma_poisson(shape = shape_a, rate = rate_a, sum_y = sum(y_a), n = n_a))

plot_gamma_poisson(
  shape = shape_a,
  rate = rate_a,
  sum_y = sum(y_a),
  n = n_a,
  prior = TRUE,
  likelihood = TRUE,
  posterior = TRUE
)

xa <- matrix(1, nrow = n_a, ncol = 1L)
colnames(xa) <- "(Intercept)"
beta_a <- matrix(shape_a / rate_a, nrow = 1L, ncol = 1L)
colnames(beta_a) <- "(Intercept)"

pf_a <- dGamma_Conjugate(shape = shape_a, rate = rate_a, beta = beta_a)

set.seed(2026)
fit_a <- rglmb(
  n = 20000,
  y = as.numeric(y_a),
  x = xa,
  family = poisson(link = "identity"),
  pfamily = pf_a,
  weights = rep(1, n_a)
)

smp_a <- fit_a$coefficients[, 1L]
post_row_a <- summarize_gamma_poisson(shape = shape_a, rate = rate_a, sum_y = sum(y_a), n = n_a)
post_row_a <- post_row_a[post_row_a$model == "posterior", ]

message(sprintf(
  "Vignette example — analytic posterior mean: %.6f ; glmbayes::rglmb() mean: %.6f",
  post_row_a$mean,
  mean(smp_a)
))

grid_a <- seq(1e-4, stats::qgamma(0.999, shape_a + sum(y_a), rate_a + n_a), length.out = 400)

plot_post_a <- ggplot(data.frame(lambda = smp_a), aes(lambda)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50, fill = "steelblue", alpha = 0.45, color = NA) +
  geom_line(
    data = data.frame(lambda = grid_a, y = stats::dgamma(grid_a, shape_a + sum(y_a), rate_a + n_a)),
    aes(lambda, y),
    inherit.aes = FALSE,
    linewidth = 1,
    colour = "black",
    linetype = "dashed"
  ) +
  labs(
    title = "glmbayes draws versus conjugate Gamma density",
    subtitle = paste0(
      "Synthetic y aligned with vignette(\"conjugate-families\", ",
      'package = "bayesrules") Gamma–Poisson example'
    ),
    x = expression(lambda),
    y = "density"
  ) +
  theme_bw()

print(plot_post_a)

## --- Example B (`bikes`: first-week ride counts) ------------------------------

data(bikes, package = "bayesrules")
n_b <- 7L
y_b <- bikes$rides[seq_len(n_b)]

## Mild conjugate Gamma prior centred near the empirical mean (illustrative only).
shape_b <- 2
rate_b <- 2 / mean(y_b)

print(summarize_gamma_poisson(shape = shape_b, rate = rate_b, sum_y = sum(y_b), n = n_b))

plot_gamma_poisson(
  shape = shape_b,
  rate = rate_b,
  sum_y = sum(y_b),
  n = n_b,
  prior = TRUE,
  likelihood = TRUE,
  posterior = TRUE
)

xb <- matrix(1, nrow = n_b, ncol = 1L)
colnames(xb) <- "(Intercept)"
beta_b <- matrix(shape_b / rate_b, nrow = 1L, ncol = 1L)
colnames(beta_b) <- "(Intercept)"

pf_b <- dGamma_Conjugate(shape = shape_b, rate = rate_b, beta = beta_b)

set.seed(7)
fit_b <- rglmb(
  n = 20000,
  y = as.numeric(y_b),
  x = xb,
  family = poisson(link = "identity"),
  pfamily = pf_b,
  weights = rep(1, n_b)
)

smp_b <- fit_b$coefficients[, 1L]
post_row_b <- summarize_gamma_poisson(shape = shape_b, rate = rate_b, sum_y = sum(y_b), n = n_b)
post_row_b <- post_row_b[post_row_b$model == "posterior", ]

message(sprintf(
  "bikes$rides — analytic posterior mean: %.6f ; glmbayes::rglmb() mean: %.6f",
  post_row_b$mean,
  mean(smp_b)
))

grid_b <- seq(1e-4, stats::qgamma(0.999, shape_b + sum(y_b), rate_b + n_b), length.out = 400)

plot_post_b <- ggplot(data.frame(lambda = smp_b), aes(lambda)) +
  geom_histogram(aes(y = after_stat(density)), bins = 55, fill = "tan2", alpha = 0.45, color = NA) +
  geom_line(
    data = data.frame(lambda = grid_b, y = stats::dgamma(grid_b, shape_b + sum(y_b), rate_b + n_b)),
    aes(lambda, y),
    inherit.aes = FALSE,
    linewidth = 1,
    colour = "brown4",
    linetype = "dashed"
  ) +
  labs(
    title = "glmbayes draws versus conjugate Gamma density",
    subtitle = "`bayesrules::bikes$rides`: first-week counts (intercept-only illustrative fit)",
    x = expression(lambda),
    y = "density"
  ) +
  theme_bw()

print(plot_post_b)

message("Finished bayesrules + glmbayes Gamma–Poisson conjugate demos.")
