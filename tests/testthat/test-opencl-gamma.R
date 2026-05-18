# inst/examples/Ex_glmb.R (carinsca Gamma block); only use_opencl = TRUE added

test_that("OpenCL f2_f3_gamma (Ex_glmb carinsca)", {
  skip_if_no_opencl()
  skip_on_cran()

  data(carinsca)
  carinsca$Merit <- ordered(carinsca$Merit)
  carinsca$Class <- factor(carinsca$Class)
  oldopt <- options(contrasts = c("contr.treatment", "contr.treatment"))
  on.exit(options(oldopt), add = TRUE)

  Claims <- carinsca$Claims
  Merit <- carinsca$Merit
  Class <- carinsca$Class
  Cost <- carinsca$Cost

  out <- glm(
    Cost / Claims ~ Merit + Class,
    family = Gamma(link = "log"),
    weights = Claims,
    x = TRUE
  )
  disp <- gamma.dispersion(out)
  ps <- Prior_Setup(
    Cost / Claims ~ Merit + Class,
    family = Gamma(link = "log"),
    weights = Claims
  )
  mu <- ps$mu
  V <- ps$Sigma

  out3 <- glmb(
    Cost / Claims ~ Merit + Class,
    family = Gamma(link = "log"),
    pfamily = dNormal(mu = mu, Sigma = V, dispersion = disp),
    weights = Claims,
    use_opencl = TRUE
  )

  expect_s3_class(out3, "glmb")
})
