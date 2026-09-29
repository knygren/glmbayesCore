# OpenCL path for matrix fitters (rglmb -> envelope C++).

test_that("rglmb Poisson with use_opencl = TRUE (Dobson RCT)", {
  skip_if_no_opencl()
  skip_on_cran()

  set.seed(333)
  counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
  outcome <- gl(3, 1, 9)
  treatment <- gl(3, 3)
  d.AD <- data.frame(treatment, outcome, counts)

  ps <- Prior_Setup(counts ~ outcome + treatment, family = poisson(), data = d.AD)
  fit <- rglmb(
    n = 40L,
    y = ps$y,
    x = as.matrix(ps$x),
    pfamily = dNormal(mu = ps$mu, Sigma = ps$Sigma),
    family = poisson(),
    use_opencl = TRUE,
    use_parallel = FALSE,
    verbose = FALSE
  )

  expect_s3_class(fit, "rglmb")
  expect_equal(nrow(fit$coefficients), 40L)
})

test_that("rglmb Poisson OpenCL path matches CPU output structure", {
  skip_if_no_opencl()
  skip_on_cran()

  set.seed(42)
  counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
  d.AD <- data.frame(
    outcome = gl(3, 1, 9),
    treatment = gl(3, 3),
    counts = counts
  )
  ps <- Prior_Setup(counts ~ outcome + treatment, family = poisson(), data = d.AD)

  common <- list(
    n = 25L,
    y = ps$y,
    x = as.matrix(ps$x),
    pfamily = dNormal(mu = ps$mu, Sigma = ps$Sigma),
    family = poisson(),
    use_parallel = FALSE,
    verbose = FALSE,
    Gridtype = 2L
  )

  fit_cpu <- do.call(rglmb, c(common, list(use_opencl = FALSE)))
  fit_gpu <- do.call(rglmb, c(common, list(use_opencl = TRUE)))

  expect_equal(dim(fit_cpu$coefficients), dim(fit_gpu$coefficients))
  expect_true(all(is.finite(fit_gpu$coefficients)))
  expect_equal(ncol(fit_gpu$coefficients), length(ps$mu))
})
