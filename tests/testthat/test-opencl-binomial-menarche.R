# Reproduces OpenCL envelope evaluation with MASS::menarche in **aggregated
# binomial** form: `cbind(successes, failures) ~ ...`, so each row has trial
# counts (not binary 0/1 with weight 1).
#
# Mirrors `inst/examples/Ex_glmb.R` (menarche logit / probit / cloglog blocks)
# with `use_opencl = TRUE` and explicit `n`, `Gridtype`, `use_parallel`, `verbose`.
#
# The Cleveland-style OpenCL tests use a Bernoulli-style response; the package
# comments in test-opencl-binomial.R note that the current f2_f3 binomial
# OpenCL path can disagree with this menarche-style setup. This file exercises
# all three links under testthat so failures surface with a clear stack trace.

menarche_opencl_ag_fit <- function(link) {
  data("menarche", package = "MASS")
  menarche$Age2 <- menarche$Age - 13
  fam <- binomial(link = link)
  ps <- Prior_Setup(
    cbind(Menarche, Total - Menarche) ~ Age2,
    family = fam,
    data = menarche
  )
  glmb(
    cbind(Menarche, Total - Menarche) ~ Age2,
    family       = fam,
    pfamily      = dNormal(mu = ps$mu, Sigma = ps$Sigma),
    data         = menarche,
    n            = 200,
    Gridtype     = 2,
    use_parallel = TRUE,
    use_opencl   = TRUE,
    verbose      = FALSE
  )
}

test_that("OpenCL binomial logit with MASS menarche (cbind successes, failures)", {
  skip_if_no_opencl()
  skip_on_cran() # parallel/OpenCL: avoids R CMD check NOTE on CPU vs elapsed time

  fit <- menarche_opencl_ag_fit("logit")
  expect_s3_class(fit, "glmb")
})

test_that("OpenCL binomial probit with MASS menarche (cbind successes, failures)", {
  skip_if_no_opencl()
  skip_on_cran()

  fit <- menarche_opencl_ag_fit("probit")
  expect_s3_class(fit, "glmb")
})

test_that("OpenCL binomial cloglog with MASS menarche (cbind successes, failures)", {
  skip_if_no_opencl()
  skip_on_cran()

  fit <- suppressWarnings(menarche_opencl_ag_fit("cloglog"))
  expect_s3_class(fit, "glmb")
})
