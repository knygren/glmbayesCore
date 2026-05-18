# inst/examples/Ex_glmb.R (Dobson Poisson block); only use_opencl = TRUE added

test_that("OpenCL f2_f3_poisson (Ex_glmb Dobson RCT)", {
  skip_if_no_opencl()
  skip_on_cran() # OpenCL: avoids R CMD check NOTE on CPU vs elapsed time

  set.seed(333)
  counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
  outcome <- gl(3, 1, 9)
  treatment <- gl(3, 3)
  d.AD <- data.frame(treatment, outcome, counts)

  ps <- Prior_Setup(counts ~ outcome + treatment, family = poisson())
  mu <- ps$mu
  V <- ps$Sigma

  glmb.D93 <- glmb(
    counts ~ outcome + treatment,
    family = poisson(),
    pfamily = dNormal(mu = mu, Sigma = V),
    use_opencl = TRUE
  )

  expect_s3_class(glmb.D93, "glmb")
})
