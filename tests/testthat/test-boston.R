test_that("Bayesian Gaussian regression with Independent Normal-Gamma prior — OpenCL", {
  
  skip_if(!has_opencl(), "OpenCL not available")
  skip_on_cran()
  
  library(MASS)
  
  # Original dataset
  data("Boston")
  
  # Identify predictors (exclude response)
  predictors <- setdiff(names(Boston), "medv")
  
  # Mean-center predictors (subtract column means, keep original units)
  Boston_centered <- Boston
  Boston_centered[predictors] <- scale(Boston[predictors], center = TRUE, scale = FALSE)
  
  # --------------------------------
  # Explicit regression formula with all predictors
  # -------------------------------
  form <- medv ~     crim + zn +     indus + chas + nox + age +  dis + rad + tax + ptratio + black + lstat+ rm 
  
  # --------------------------------
  # Conjugate Normal prior (fixed dispersion)
  # --------------------------------
  ps <- Prior_Setup(form, gaussian(), data = Boston_centered)
  p <- ncol(ps$x)

  lmb.boston <- lmb(form,data     = Boston_centered,
    pfamily  = dNormal(mu = ps$mu,
                       Sigma = ps$Sigma,
                       dispersion = ps$dispersion)
  )

  # --------------------------------
  # Conjugate Normal–Gamma prior
  # --------------------------------
  lmb.boston_v2 <- lmb(
    form,
    data     = Boston_centered,
    pfamily  = dNormal_Gamma(mu     = ps$mu,
                             Sigma_0 = ps$Sigma_0,
                             shape  = ps$shape,
                             rate   = ps$rate)
  )
  

  # --------------------------------
  # Independent Normal–Gamma prior (non-conjugate, uses envelope)
  # --------------------------------
  lmb.boston_v3 <- glmb(n = 1000,
                        form,
                        data       = Boston_centered,
                        family     = gaussian(),
                        pfamily    = dIndependent_Normal_Gamma(ps$mu, ps$Sigma,
                                                               shape = ps$shape + p / 2,
                                                               rate  = ps$rate
                        ),
                        use_parallel = TRUE,
                        use_opencl = TRUE,
                        verbose    = FALSE
  )
  
  
  # Basic acceptance diagnostic:
  # Expect mean candidates per acceptance to be reasonably small (< 10)
  avg_candidates <- mean(lmb.boston_v3$iters)
  expect_true(avg_candidates < 400)
  
})