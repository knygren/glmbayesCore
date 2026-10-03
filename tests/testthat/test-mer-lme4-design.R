test_that("extract_re_hyper_matrices: sleepstudy design structure", {
  dat <- lme4::sleepstudy
  f <- Reaction ~ Days + (Days || Subject)

  out <- glmbayesCore:::extract_re_hyper_matrices(f, data = dat)

  expect_identical(out$group_name, "Subject")
  expect_identical(out$groupef.names, c("(Intercept)", "Days"))
  expect_length(out$y, nrow(dat))
  expect_equal(nrow(out$D), length(out$y))
  J <- nlevels(out$group)
  expect_equal(nrow(out$W[[1L]]), J)
  expect_equal(nrow(out$W[[2L]]), J)
  expect_identical(as.character(out$group), as.character(dat$Subject))
  expect_identical(attr(out$group, "group_name"), "Subject")
  expect_length(out$weights, length(out$y))
  expect_length(out$offset, length(out$y))
})

test_that("is_single_factor_model: single vs multiple grouping factors", {
  dat <- lme4::sleepstudy
  f_ok <- Reaction ~ Days + (Days || Subject)
  f_bad <- Reaction ~ Days + (Days || Subject) + (1 | Dummy)

  expect_true(glmbayesCore:::is_single_factor_model(f_ok, data = dat))
  expect_false(glmbayesCore:::is_single_factor_model(f_bad, data = dat))
})

test_that("model_setup: sleepstudy design and reference fit", {
  dat <- lme4::sleepstudy
  f <- Reaction ~ Days + (Days || Subject)

  ms <- glmbayesCore::model_setup(formula = f, data = dat, REML = TRUE)

  expect_s4_class(ms$lmer, "merMod")
  expect_identical(ms$group_name, "Subject")
  expect_identical(ms$groupef.names, c("(Intercept)", "Days"))
  expect_true(ms$popef.rank_ok)
  expect_type(ms$Psi, "double")
  expect_true(all(is.finite(ms$Psi)))
})

test_that("model_setup: dispformula ~ Subject stores glmmTMB reference fit", {
  skip_if_not_installed("glmmTMB")
  dat <- lme4::sleepstudy
  f <- Reaction ~ Days + (Days || Subject)

  ms <- glmbayesCore::model_setup(
    formula = f,
    data = dat,
    dispformula = ~Subject
  )

  expect_false(is.null(ms$glmmTMB_fit))
  expect_s3_class(ms$glmmTMB_fit, "glmmTMB")
  expect_s4_class(ms$lmer, "merMod")
})

test_that("check_identifiability on model_setup design (sleepstudy)", {
  ms <- glmbayesCore::model_setup(
    Reaction ~ Days + (Days || Subject),
    data = lme4::sleepstudy
  )
  id <- glmbayesCore::check_identifiability(
    y = ms$y,
    D = ms$D,
    group = ms$group,
    W = ms$W,
    family = gaussian(),
    group_name = ms$group_name
  )

  expect_length(id$groupef.rank, nlevels(ms$group))
  expect_true(all(id$groupef.rank))
  expect_true(id$popef.rank_ok)
})
