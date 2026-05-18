# OpenCL binomial tests follow inst/examples/Ex_Cleveland.R (the package's
# documented OpenCL logistic example). Ex_glmb.R menarche uses cbind(success, failure)
# with trial weights; the current f2_f3 binomial OpenCL kernels expect the Cleveland
# style (0/1 response, weight 1) and fail PLSD on menarche even with use_opencl alone.
#
# Probit and cloglog: same Ex_Cleveland.R setup; link from Ex_glmb.R menarche block.

test_that("OpenCL f2_f3_binomial_logit (Ex_Cleveland)", {
  skip_if_no_opencl()
  skip_on_cran() # parallel/OpenCL: avoids R CMD check NOTE on CPU vs elapsed time

  data("Cleveland")

  ps <- Prior_Setup(
    hd ~ age + sex + cp + trestbps + chol +
      fbs + restecg + thalach + exang + oldpeak + slope + ca + thal,
    family = binomial(logit),
    data = Cleveland
  )

  fit <- glmb(
    hd ~ age + sex + cp + trestbps + chol +
      fbs + restecg + thalach + exang + oldpeak + slope + ca + thal,
    family       = binomial(link = "logit"),
    pfamily      = dNormal(mu = ps$mu, Sigma = ps$Sigma),
    data         = Cleveland,
    n            = 1000,
    Gridtype     = 2,
    use_parallel = TRUE,
    use_opencl   = TRUE,
    verbose      = FALSE
  )

  expect_s3_class(fit, "glmb")
})

test_that("OpenCL f2_f3_binomial_probit (Ex_Cleveland setup, Ex_glmb link)", {
  skip_if_no_opencl()
  skip_on_cran()

  data("Cleveland")

  ps <- Prior_Setup(
    hd ~ age + sex + cp + trestbps + chol +
      fbs + restecg + thalach + exang + oldpeak + slope + ca + thal,
    family = binomial(probit),
    data = Cleveland
  )

  fit <- glmb(
    hd ~ age + sex + cp + trestbps + chol +
      fbs + restecg + thalach + exang + oldpeak + slope + ca + thal,
    family       = binomial(link = "probit"),
    pfamily      = dNormal(mu = ps$mu, Sigma = ps$Sigma),
    data         = Cleveland,
    n            = 1000,
    Gridtype     = 2,
    use_parallel = TRUE,
    use_opencl   = TRUE,
    verbose      = FALSE
  )

  expect_s3_class(fit, "glmb")
})

test_that("OpenCL f2_f3_binomial_cloglog (Ex_Cleveland setup, Ex_glmb link)", {
  skip_if_no_opencl()
  skip_on_cran()

  data("Cleveland")

  fit <- suppressWarnings({
    ps <- Prior_Setup(
      hd ~ age + sex + cp + trestbps + chol +
        fbs + restecg + thalach + exang + oldpeak + slope + ca + thal,
      family = binomial(cloglog),
      data = Cleveland
    )

    glmb(
      hd ~ age + sex + cp + trestbps + chol +
        fbs + restecg + thalach + exang + oldpeak + slope + ca + thal,
      family       = binomial(link = "cloglog"),
      pfamily      = dNormal(mu = ps$mu, Sigma = ps$Sigma),
      data         = Cleveland,
      n            = 1000,
      Gridtype     = 2,
      use_parallel = TRUE,
      use_opencl   = TRUE,
      verbose      = FALSE
    )
  })

  expect_s3_class(fit, "glmb")
})
