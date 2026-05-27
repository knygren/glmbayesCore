## Tests for dBeta pfamily + rBeta_reg + glmbfamfunc(binomial(identity))
## -----------------------------------------------------------------------
## Covers:
##   1. dBeta() constructor validation
##   2. Conjugate draw mean/SD match analytic posterior
##   3. summary(), logLik(), DIC_Info() run without error
##   4. Prior_Setup() produces non-NULL conj_binomial

test_that("dBeta() rejects bad inputs", {
  b <- matrix(0.5, nrow = 1, ncol = 1, dimnames = list(NULL, "(Intercept)"))
  expect_error(dBeta(shape1 = -1, shape2 = 2,  beta = b), "positive")
  expect_error(dBeta(shape1 =  2, shape2 = -1, beta = b), "positive")
  expect_error(dBeta(shape1 =  0, shape2 =  2, beta = b), "positive")
  expect_error(dBeta(shape1 = "a", shape2 = 2, beta = b), "non-numeric")
  expect_error(dBeta(shape1 = c(1, 2), shape2 = 2, beta = b), "single")
})

test_that("dBeta() constructor returns valid pfamily object", {
  b  <- matrix(0.5, nrow = 1, ncol = 1, dimnames = list(NULL, "(Intercept)"))
  pf <- dBeta(shape1 = 2, shape2 = 2, beta = b)
  expect_s3_class(pf, "pfamily")
  expect_equal(pf$pfamily, "dBeta")
  expect_true("binomial" %in% pf$okfamilies)
  expect_equal(pf$plinks(binomial(link = "identity")), "identity")
  expect_null(pf$plinks(poisson()))
  pl <- pf$prior_list
  expect_equal(as.numeric(pl$mu), 0.5)          ## Beta(2,2) mean
  expect_true(as.numeric(pl$Sigma) > 0)
})

test_that("glmb() with dBeta draws match analytic Beta posterior", {
  skip_if_not_installed("glmbayes")

  set.seed(101)
  n_obs     <- 40L
  theta_true <- 0.3
  ## Individual binary outcomes: each row is one Bernoulli trial (n_i = 1).
  ## Passing weights = 1L per row makes this explicit to glmb() / glm().
  y_dat  <- rbinom(n_obs, size = 1, prob = theta_true)

  alpha0 <- 3;  beta0 <- 7   ## Beta(3,7) prior: mean = 0.3

  ## Analytic posterior
  s1_post <- alpha0 + sum(y_dat)
  s2_post <- beta0  + (n_obs - sum(y_dat))
  post_mean_analytic <- s1_post / (s1_post + s2_post)
  post_sd_analytic   <- sqrt(s1_post * s2_post /
                               ((s1_post + s2_post)^2 * (s1_post + s2_post + 1)))

  b_init <- matrix(alpha0 / (alpha0 + beta0), nrow = 1L, ncol = 1L,
                   dimnames = list(NULL, "(Intercept)"))
  pf <- dBeta(shape1 = alpha0, shape2 = beta0, beta = b_init)

  set.seed(2026)
  fit <- glmb(
    n       = 30000,
    y ~ 1,
    data    = data.frame(y = y_dat),
    weights = rep(1L, n_obs),
    family  = binomial(link = "identity"),
    pfamily = pf
  )

  smp <- fit$coefficients[, 1L]

  expect_equal(mean(smp), post_mean_analytic, tolerance = 0.006)
  expect_equal(sd(smp),   post_sd_analytic,   tolerance = 0.006)
})

test_that("summary(), logLik(), DIC_Info() run without error for dBeta fit", {
  skip_if_not_installed("glmbayes")

  set.seed(7)
  ## 8 successes, 12 failures — supply as a two-column cbind response so the
  ## number of trials is unambiguous (one aggregated row, n = 20 trials).
  y_mat  <- cbind(success = 8L, failure = 12L)
  b_init <- matrix(0.5, nrow = 1, ncol = 1, dimnames = list(NULL, "(Intercept)"))
  pf     <- dBeta(shape1 = 2, shape2 = 2, beta = b_init)

  set.seed(2026)
  fit <- glmb(
    n       = 5000,
    y_mat ~ 1,
    family  = binomial(link = "identity"),
    pfamily = pf
  )

  sm <- summary(fit)
  expect_no_error(print(sm))
  expect_no_error(logLik(fit))
  ## DIC is computed inside summary(); verify the slot is present and finite
  expect_true(is.finite(sm$DIC))
})

test_that("Prior_Setup() produces non-NULL conj_binomial for binomial(identity)", {
  skip_if_not_installed("glmbayes")

  ## Represent as 25 individual binary trials with explicit weights = 1 per trial,
  ## so R knows these are single-trial Binomial outcomes (not ambiguous proportions).
  y_dat <- c(rep(1, 7), rep(0, 18))    ## 7 successes, 18 failures
  df    <- data.frame(y = y_dat)

  ps <- Prior_Setup(y ~ 1, data = df,
                    weights = rep(1L, nrow(df)),
                    family  = binomial(link = "identity"),
                    pwt     = 0.05)

  expect_false(is.null(ps$conj_binomial))
  cb <- ps$conj_binomial
  expect_true(cb$shape1 > 0)
  expect_true(cb$shape2 > 0)
  ## Prior mean should equal the weighted data proportion
  expect_equal(cb$shape1 / (cb$shape1 + cb$shape2),
               cb$weighted_mean_prop, tolerance = 1e-10)
})
