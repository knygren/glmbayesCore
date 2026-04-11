## Non-Zellner prior: Independent Normal-Gamma with a strongly shrunk diagonal
## covariance (0.001 * diag(diag(ps$Sigma))) instead of full Prior_Setup Sigma.
## Formerly in inst/examples/Ex_lmb.R as a "temporary non zellner test".

test_that("lmb: Independent Normal-Gamma with scaled diagonal Sigma (non-Zellner)", {
  ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
  trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
  group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
  weight <- c(ctl, trt)

  ps <- Prior_Setup(weight ~ group, gaussian())
  p <- ncol(ps$x)
  Sigma_non_zellner <- 0.001 * diag(diag(ps$Sigma))

  fit <- lmb(
    weight ~ group,
    dIndependent_Normal_Gamma(
      ps$mu,
      Sigma_non_zellner,
      shape = ps$shape + p / 2,
      rate  = ps$rate
    ),
    n = 500L,
    verbose = FALSE,
    use_parallel = FALSE
  )

  expect_s3_class(fit, "lmb")
  expect_equal(nrow(fit$coefficients), 500L)
  expect_equal(ncol(fit$coefficients), nrow(ps$mu))
  expect_true(all(is.finite(fit$coef.means)))
})

test_that("Prior_Setup Gaussian calibration: shape from n_prior, E[sigma^2|y] = dispersion", {
  ctl <- c(4.17, 5.58)
  trt <- c(4.81, 4.17)
  group <- gl(2, 2, 4)
  weight <- c(ctl, trt)
  ps <- Prior_Setup(
    weight ~ group,
    gaussian(),
    pwt = 0.01
  )
  p <- ncol(ps$x)
  n_w <- ps$PriorSettings$n_effective
  n_prior <- ps$PriorSettings$n_prior
  S_marg <- ps$dispersion * (n_w - p)
  expect_equal(ps$shape, (n_prior + 1) / 2)
  post_mean_sigma2 <- (ps$rate + S_marg / 2) / (ps$shape + n_w / 2 - 1)
  expect_equal(post_mean_sigma2, ps$dispersion)
})
