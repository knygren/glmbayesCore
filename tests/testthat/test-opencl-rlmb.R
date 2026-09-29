test_that("rlmb Gaussian with use_opencl = TRUE (plant weights)", {
  skip_if_no_opencl()
  skip_on_cran()

  ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
  trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
  group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
  weight <- c(ctl, trt)
  d <- data.frame(weight, group)

  ps <- Prior_Setup(weight ~ group, gaussian(), data = d)
  fit <- rlmb(
    n = 50L,
    y = ps$y,
    x = as.matrix(ps$x),
    pfamily = dNormal(mu = ps$mu, Sigma = ps$Sigma, dispersion = ps$dispersion),
    use_opencl = TRUE,
    use_parallel = FALSE,
    verbose = FALSE
  )

  expect_s3_class(fit, "rlmb")
  expect_equal(nrow(fit$coefficients), 50L)
})

test_that("rlmb Gaussian OpenCL path matches CPU output structure", {
  skip_if_no_opencl()
  skip_on_cran()

  set.seed(7)
  ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
  trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
  group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
  weight <- c(ctl, trt)
  d <- data.frame(weight, group)

  ps <- Prior_Setup(weight ~ group, gaussian(), data = d)
  pf <- dNormal(mu = ps$mu, Sigma = ps$Sigma, dispersion = ps$dispersion)
  common <- list(
    n = 30L,
    y = ps$y,
    x = as.matrix(ps$x),
    pfamily = pf,
    use_parallel = FALSE,
    verbose = FALSE,
    Gridtype = 2L
  )

  fit_cpu <- do.call(rlmb, c(common, list(use_opencl = FALSE)))
  fit_gpu <- do.call(rlmb, c(common, list(use_opencl = TRUE)))

  expect_equal(dim(fit_cpu$coefficients), dim(fit_gpu$coefficients))
  expect_true(all(is.finite(fit_gpu$coefficients)))
})
