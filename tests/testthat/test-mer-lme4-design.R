test_that("extract_re_hyper_matrices matches lmebayesCore (sleepstudy)", {
  skip_if_not_installed("lmebayesCore")

  dat <- lme4::sleepstudy
  f <- Reaction ~ Days + (Days || Subject)

  old <- lmebayesCore:::extract_re_hyper_matrices(f, data = dat)
  new <- glmbayesCore:::extract_re_hyper_matrices(f, data = dat)

  expect_identical(old$group_name, new$group_name)
  expect_identical(old$groupef.names, new$groupef.names)
  expect_equal(old$y, new$y)
  expect_equal(old$weights, new$weights)
  expect_equal(old$offset, new$offset)
  expect_equal(old$D, new$D)
  expect_equal(old$W, new$W)
  expect_equal(old$popef.moderation, new$popef.moderation)
  expect_identical(old$group, new$group)
})

test_that("is_single_factor_model agrees with lmebayesCore", {
  skip_if_not_installed("lmebayesCore")

  dat <- lme4::sleepstudy
  f_ok <- Reaction ~ Days + (Days || Subject)
  f_bad <- Reaction ~ Days + (Days || Subject) + (1 | Dummy)

  expect_identical(
    lmebayesCore:::is_single_factor_model(f_ok, data = dat),
    glmbayesCore:::is_single_factor_model(f_ok, data = dat)
  )
  expect_identical(
    lmebayesCore:::is_single_factor_model(f_bad, data = dat),
    glmbayesCore:::is_single_factor_model(f_bad, data = dat)
  )
})

test_that("model_setup matches lmebayesCore (sleepstudy)", {
  skip_if_not_installed("lmebayesCore")

  dat <- lme4::sleepstudy
  f <- Reaction ~ Days + (Days || Subject)
  args <- list(formula = f, data = dat, REML = TRUE)

  old <- do.call(lmebayesCore::model_setup, args)
  new <- do.call(glmbayesCore::model_setup, args)

  expect_equal(old$y, new$y)
  expect_equal(old$D, new$D)
  expect_equal(old$W, new$W)
  expect_equal(old$Psi, new$Psi, tolerance = 1e-6)
  expect_equal(old$dispersion, new$dispersion, tolerance = 1e-6)
  expect_identical(old$groupef.names, new$groupef.names)
  expect_identical(old$groupef.rank, new$groupef.rank)
  expect_identical(old$popef.rank_ok, new$popef.rank_ok)
})

test_that("model_setup: dispformula ~ Subject matches lmebayesCore (glmmTMB)", {
  skip_if_not_installed("lmebayesCore")

  dat <- lme4::sleepstudy
  f <- Reaction ~ Days + (Days || Subject)
  args <- list(
    formula = f,
    data = dat,
    dispformula = ~Subject
  )

  old <- do.call(lmebayesCore::model_setup, args)
  new <- do.call(glmbayesCore::model_setup, args)

  expect_false(is.null(old[["glmmTMB_fit"]]))
  expect_false(is.null(new[["glmmTMB_fit"]]))
  expect_equal(
    stats::fitted(old[["glmmTMB_fit"]]),
    stats::fitted(new[["glmmTMB_fit"]]),
    tolerance = 1e-5
  )
})

test_that("check_identifiability matches lmebayesCore on model_setup design", {
  skip_if_not_installed("lmebayesCore")

  ms <- glmbayesCore::model_setup(
    Reaction ~ Days + (Days || Subject),
    data = lme4::sleepstudy
  )
  args <- list(
    y = ms$y,
    D = ms$D,
    group = ms$group,
    W = ms$W,
    family = gaussian(),
    group_name = ms$group_name
  )
  old <- do.call(lmebayesCore::check_identifiability, args)
  new <- do.call(glmbayesCore::check_identifiability, args)

  expect_identical(old$groupef.rank, new$groupef.rank)
  expect_identical(old$popef.rank_ok, new$popef.rank_ok)
})
