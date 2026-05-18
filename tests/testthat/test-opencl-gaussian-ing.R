# inst/examples/Ex_Boston_centered.R (Independent Normal-Gamma block);
# glmb as in tests/testthat/test-boston.R; only use_opencl = TRUE added

test_that("OpenCL f2_f3_gaussian ING (Ex_Boston_centered)", {
  skip_if_no_opencl()
  skip_on_cran()

  data("Boston_centered")

  form <- medv ~
    crim + zn + indus + chas + nox + age + dis + rad + tax + ptratio + black + lstat + rm

  ps <- Prior_Setup(form, gaussian(), data = Boston_centered)

  fit <- glmb(
    n = 1000,
    form,
    data = Boston_centered,
    family = gaussian(),
    pfamily = dIndependent_Normal_Gamma(
      ps$mu,
      ps$Sigma,
      shape = ps$shape_ING,
      rate = ps$rate
    ),
    use_opencl = TRUE
  )

  expect_s3_class(fit, "glmb")
})
