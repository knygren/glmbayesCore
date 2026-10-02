# Step 4: Prior_Setup_GLMM tests (glmbayesCore; parity vs lmebayesCore optional).

test_that("Prior_Setup_GLMM: pooled group.dispersion is A12 3.3.4 S_marg center", {
  dat <- lme4::sleepstudy

  ps <- glmbayesCore::Prior_Setup_GLMM(
    Reaction ~ Days + (Days || Subject),
    data = dat,
    pop.pwt = 0.01
  )

  expect_identical(attr(ps$group.dispersion, "source"), "pooled S_marg (A12 3.3.4)")
  expect_true(is.finite(attr(ps$group.dispersion, "classical")))
  expect_equal(
    ps$group.ing_prior$sigma2_hat,
    as.numeric(ps$group.dispersion),
    ignore_attr = TRUE
  )
  p_re <- length(ps$design$groupef.names)
  n_prior <- ps$group.ing_prior$n_prior
  expect_equal(
    ps$group.ing_prior$rate,
    as.numeric(ps$group.dispersion) * (n_prior + p_re - 1) / 2
  )
})

test_that("Prior_Setup_GLMM: per-group group.dispersion follows sigma2_hat", {
  skip_if_not_installed("glmmTMB")
  dat <- lme4::sleepstudy
  J <- nlevels(dat$Subject)

  ps <- glmbayesCore::Prior_Setup_GLMM(
    Reaction ~ Days + (Days || Subject),
    data = dat,
    pop.pwt = 0.01,
    group.dispersion.pwt = 0.1,
    group.alpha_target = NULL,
    dispformula = ~Subject
  )

  expect_length(ps$group.dispersion, J)
  expect_named(ps$group.dispersion, levels(dat$Subject))
  expect_identical(
    attr(ps$group.dispersion, "source"),
    "per-group S_marg (A12 3.3.4)"
  )
  for (lev in levels(dat$Subject)) {
    expect_equal(
      unname(ps$group.dispersion[[lev]]),
      ps$group.ing_prior[[lev]]$sigma2_hat
    )
  }
})

test_that("model_setup: subset/na.action align design rows with reference fit", {
  dat <- lme4::sleepstudy
  dat$Reaction[1:5] <- NA_real_

  ms <- glmbayesCore::model_setup(
    Reaction ~ Days + (Days || Subject),
    data = dat,
    na.action = stats::na.omit
  )
  expect_equal(length(ms$y), length(stats::fitted(ms$lmer)))
  expect_equal(length(ms$y), nrow(dat) - 5L)
  expect_equal(nlevels(ms$group), nlevels(dat$Subject))

  ms_sub <- glmbayesCore::model_setup(
    Reaction ~ Days + (Days || Subject),
    data = lme4::sleepstudy,
    subset = Days >= 5
  )
  expect_equal(length(ms_sub$y), length(stats::fitted(ms_sub$lmer)))
  expect_true(length(ms_sub$y) < nrow(lme4::sleepstudy))
  expect_equal(length(ms_sub$y), sum(lme4::sleepstudy$Days >= 5))
})

test_that("model_setup: contrasts align design W with reference fixef coding", {
  set.seed(1)
  dat <- data.frame(
    y = rnorm(40),
    g = gl(10, 4),
    f = factor(rep(c("A", "B"), each = 20))
  )
  ## f must be group-constant (level-2) for glmbayesCore::model_setup().
  dat$f <- factor(ifelse(as.integer(dat$g) %% 2L == 0L, "A", "B"))

  ms_tr <- glmbayesCore::model_setup(y ~ f + (1 | g), data = dat)
  ms_sum <- glmbayesCore::model_setup(
    y ~ f + (1 | g),
    data = dat,
    contrasts = list(f = "contr.sum")
  )
  expect_equal(length(ms_tr$y), length(stats::fitted(ms_tr$lmer)))
  expect_equal(length(ms_sum$y), length(stats::fitted(ms_sum$lmer)))
  ## Treatment vs sum contrasts change the fixed-effect column names on W
  ## and on fixef() of the reference fit -- and those names must agree.
  expect_identical(
    colnames(ms_tr$W[["(Intercept)"]]),
    names(lme4::fixef(ms_tr$lmer))
  )
  expect_identical(
    colnames(ms_sum$W[["(Intercept)"]]),
    names(lme4::fixef(ms_sum$lmer))
  )
  expect_false(
    identical(
      colnames(ms_tr$W[["(Intercept)"]]),
      colnames(ms_sum$W[["(Intercept)"]])
    )
  )
})

test_that("Prior_Setup_GLMM: forwards REML / weights to model_setup", {
  dat <- lme4::sleepstudy

  ps_reml <- glmbayesCore::Prior_Setup_GLMM(
    Reaction ~ Days + (Days || Subject),
    data = dat,
    pop.pwt = 0.01,
    REML = TRUE
  )
  ps_ml <- glmbayesCore::Prior_Setup_GLMM(
    Reaction ~ Days + (Days || Subject),
    data = dat,
    pop.pwt = 0.01,
    REML = FALSE
  )
  ## ML vs REML changes the classical reference used for calibration.
  expect_false(isTRUE(all.equal(ps_reml$group.Sigma, ps_ml$group.Sigma)))

  ## weights/offset are stored on design (Phase 1); reference fit uses them.
  w <- rep(1, nrow(dat))
  w[1:10] <- 2
  off <- rep(0, nrow(dat))
  off[1:5] <- 1
  ps_w <- glmbayesCore::Prior_Setup_GLMM(
    Reaction ~ Days + (Days || Subject),
    data = dat,
    pop.pwt = 0.01,
    weights = w,
    offset = off
  )
  expect_s3_class(ps_w, "Prior_Setup_GLMM")
  expect_true(is.numeric(ps_w$group.dispersion))
  expect_equal(ps_w$design$weights, w)
  expect_equal(ps_w$design$offset, off)
  expect_equal(length(ps_w$design$weights), length(ps_w$design$y))

  ms_default <- glmbayesCore::model_setup(
    Reaction ~ Days + (Days || Subject),
    data = dat
  )
  expect_equal(ms_default$weights, rep(1, length(ms_default$y)))
  expect_equal(ms_default$offset, rep(0, length(ms_default$y)))

  ## subset expression is NSE-forwarded (same rows as model_setup); weights
  ## are subset with the model frame.
  w_sub <- rep(1, nrow(dat))
  w_sub[dat$Days >= 5] <- 2
  ps_sub <- glmbayesCore::Prior_Setup_GLMM(
    Reaction ~ Days + (Days || Subject),
    data = dat,
    pop.pwt = 0.01,
    subset = Days >= 5,
    weights = w_sub
  )
  expect_equal(length(ps_sub$design$y), sum(dat$Days >= 5))
  expect_equal(
    length(ps_sub$design$y),
    length(stats::fitted(ps_sub$design$lmer))
  )
  expect_equal(ps_sub$design$weights, rep(2, length(ps_sub$design$y)))
  expect_error(
    .lmebayes_stop_if_nondefault_weights_offset(
      ps_sub$design$weights,
      ps_sub$design$offset,
      where = "test"
    ),
    "non-unit 'weights'"
  )
})

test_that("Prior_Setup_GLMM: dispformula = ~1 keeps the pooled lme4 reference", {
  dat <- lme4::sleepstudy

  ps <- glmbayesCore::Prior_Setup_GLMM(
    Reaction ~ Days + (Days || Subject),
    data = dat,
    pop.pwt  = 0.01
  )

  expect_identical(ps$calibration_source, "lme4")
  expect_s4_class(ps$fit_ref, "merMod")
  expect_identical(ps$fit_ref, ps$mer_fit)
  expect_identical(ps$fit_ref, ps$design$lmer)
  expect_null(ps$group.dispersion.fit)
  expect_null(ps$group.dispersion.ref)
  expect_false(is.null(ps$group.ing_prior))
  expect_false(lmebayesCore:::.lmebayes_ing_prior_is_grouped(ps$group.ing_prior))
})

test_that("Prior_Setup_GLMM: dispformula = ~group routes calibration through glmmTMB", {
  skip_if_not_installed("glmmTMB")
  dat <- lme4::sleepstudy

  ps <- glmbayesCore::Prior_Setup_GLMM(
    Reaction ~ Days + (Days || Subject),
    data            = dat,
    pop.pwt         = 0.01,
    group.dispersion.pwt       = 0.1,
    dispformula     = ~Subject
  )

  expect_identical(ps$calibration_source, "glmmTMB")
  expect_s3_class(ps$fit_ref, "glmmTMB")
  expect_s4_class(ps$mer_fit, "merMod")
  expect_identical(ps$group.dispersion.fit, ps$fit_ref)
  expect_identical(ps$mer_fit, ps$design$lmer)

  ## Per-group dispformula: group.dispersion is the named vector of
  ## per-group S_marg centers (not the pooled lmer residual).
  expect_identical(
    attr(ps$group.dispersion, "source"),
    "per-group S_marg (A12 3.3.4)"
  )
  expect_length(ps$group.dispersion, nlevels(dat$Subject))
  expect_named(ps$group.dispersion, levels(dat$Subject))

  ## group.ing_prior holds either the pooled (flat) or per-group (named
  ## list of lists) calibration, selected by dispformula: only the
  ## per-group shape is calibrated here.
  expect_false(is.null(ps$group.ing_prior))
  expect_true(lmebayesCore:::.lmebayes_ing_prior_is_grouped(ps$group.ing_prior))
  expect_length(ps$group.ing_prior, nlevels(dat$Subject))
  for (lev in levels(dat$Subject)) {
    expect_equal(
      unname(ps$group.dispersion[[lev]]),
      ps$group.ing_prior[[lev]]$sigma2_hat
    )
  }

  expect_type(ps$group.dispersion.ref, "double")
  expect_named(ps$group.dispersion.ref, levels(dat$Subject))

  disp_pf <- glmbayesCore::dGamma_list(ps, warn_asymmetric = FALSE)
  expect_length(disp_pf, nlevels(dat$Subject))
  expect_true(all(vapply(disp_pf, inherits, logical(1L), "pfamily")))
  expect_identical(attr(disp_pf, "group.dispersion.fit"), ps$group.dispersion.fit)
  expect_identical(attr(disp_pf, "calibration_source"), "glmmTMB")
})

test_that("Prior_Setup_GLMM: pop.mu supports partial (NULL-entry) overrides", {
  dat <- lme4::sleepstudy

  ps_default <- glmbayesCore::Prior_Setup_GLMM(
    Reaction ~ Days + (Days || Subject),
    data = dat
  )

  ps_mu <- glmbayesCore::Prior_Setup_GLMM(
    Reaction ~ Days + (Days || Subject),
    data = dat,
    pop.mu = list("(Intercept)" = c("(Intercept)" = 250), Days = NULL)
  )

  ## Overridden component: mu is exactly the supplied value.
  expect_equal(
    unname(ps_mu$pop.prior_list[["(Intercept)"]]$mu),
    250
  )

  ## Non-overridden component (NULL entry): unchanged from the usual
  ## pop.intercept_source/pop.effects_source-derived default.
  expect_identical(
    ps_mu$pop.prior_list[["Days"]]$mu,
    ps_default$pop.prior_list[["Days"]]$mu
  )

  ## Sigma is untouched by pop.mu.
  expect_identical(
    ps_mu$pop.prior_list[["(Intercept)"]]$Sigma,
    ps_default$pop.prior_list[["(Intercept)"]]$Sigma
  )
})

test_that("Prior_Setup_GLMM: pop.sd converts to the expected Sigma scaling", {
  dat <- lme4::sleepstudy

  sd_val <- 50
  ps_sd <- glmbayesCore::Prior_Setup_GLMM(
    Reaction ~ Days + (Days || Subject),
    data = dat,
    pop.sd = list("(Intercept)" = sd_val, Days = NULL)
  )
  ps_default <- glmbayesCore::Prior_Setup_GLMM(
    Reaction ~ Days + (Days || Subject),
    data = dat
  )

  ## For a single-column component, Sigma = V_ii * (1-w)/w with
  ## w = V_ii / (V_ii + sd^2) collapses to exactly sd^2, independent of V_ii.
  expect_equal(
    ps_sd$pop.prior_list[["(Intercept)"]]$Sigma[1, 1],
    sd_val^2
  )

  ## Non-overridden component (NULL entry): unchanged pop.pwt-derived scaling.
  expect_identical(
    ps_sd$pop.prior_list[["Days"]]$Sigma,
    ps_default$pop.prior_list[["Days"]]$Sigma
  )

  ## mu is untouched by pop.sd.
  expect_identical(
    ps_sd$pop.prior_list[["(Intercept)"]]$mu,
    ps_default$pop.prior_list[["(Intercept)"]]$mu
  )
})

test_that("Prior_Setup_GLMM: pop.pwt, pop.sd, and pop.nprior are mutually exclusive", {
  dat <- lme4::sleepstudy

  expect_error(
    glmbayesCore::Prior_Setup_GLMM(
      Reaction ~ Days + (Days || Subject),
      data = dat,
      pop.pwt = 0.02,
      pop.sd  = list("(Intercept)" = 50, Days = NULL)
    ),
    "at most one of 'pop.pwt', 'pop.sd', and 'pop.nprior'"
  )

  expect_error(
    glmbayesCore::Prior_Setup_GLMM(
      Reaction ~ Days + (Days || Subject),
      data = dat,
      pop.pwt    = 0.02,
      pop.nprior = list("(Intercept)" = 5, Days = NULL)
    ),
    "at most one of 'pop.pwt', 'pop.sd', and 'pop.nprior'"
  )

  expect_error(
    glmbayesCore::Prior_Setup_GLMM(
      Reaction ~ Days + (Days || Subject),
      data = dat,
      pop.sd     = list("(Intercept)" = 50, Days = NULL),
      pop.nprior = list("(Intercept)" = 5, Days = NULL)
    ),
    "at most one of 'pop.pwt', 'pop.sd', and 'pop.nprior'"
  )

  ## Relying on the pop.pwt default (not explicitly supplied) does not
  ## conflict with an explicit pop.sd/pop.nprior.
  expect_no_error(
    glmbayesCore::Prior_Setup_GLMM(
      Reaction ~ Days + (Days || Subject),
      data = dat,
      pop.sd = list("(Intercept)" = 50, Days = NULL)
    )
  )
  expect_no_error(
    glmbayesCore::Prior_Setup_GLMM(
      Reaction ~ Days + (Days || Subject),
      data = dat,
      pop.nprior = list("(Intercept)" = 5, Days = NULL)
    )
  )
})

test_that("Prior_Setup_GLMM: pop.nprior converts to the expected Sigma scaling", {
  dat <- lme4::sleepstudy
  J   <- nlevels(dat$Subject)

  n_val <- 5
  ps_nprior <- glmbayesCore::Prior_Setup_GLMM(
    Reaction ~ Days + (Days || Subject),
    data = dat,
    pop.nprior = list("(Intercept)" = n_val, Days = NULL)
  )
  ps_default <- glmbayesCore::Prior_Setup_GLMM(
    Reaction ~ Days + (Days || Subject),
    data = dat
  )

  ## n_i = J*w_i/(1-w_i) <=> w_i = n_i/(n_i+J): an equivalent pop.pwt list
  ## must reproduce exactly the same Sigma for the overridden component.
  w_equiv <- n_val / (n_val + J)
  ps_pwt_equiv <- glmbayesCore::Prior_Setup_GLMM(
    Reaction ~ Days + (Days || Subject),
    data = dat,
    pop.pwt = list("(Intercept)" = w_equiv, Days = 0.01)
  )
  expect_equal(
    ps_nprior$pop.prior_list[["(Intercept)"]]$Sigma,
    ps_pwt_equiv$pop.prior_list[["(Intercept)"]]$Sigma
  )

  ## Non-overridden component (NULL entry): unchanged pop.pwt-derived scaling.
  expect_identical(
    ps_nprior$pop.prior_list[["Days"]]$Sigma,
    ps_default$pop.prior_list[["Days"]]$Sigma
  )

  ## mu is untouched by pop.nprior.
  expect_identical(
    ps_nprior$pop.prior_list[["(Intercept)"]]$mu,
    ps_default$pop.prior_list[["(Intercept)"]]$mu
  )

  ## pop.nprior output reproduces the supplied value for the overridden
  ## component, and the pop.pwt-derived default for the other.
  expect_equal(unname(ps_nprior$pop.nprior[["(Intercept)"]]), n_val)
  expect_equal(ps_nprior$pop.nprior[["Days"]], ps_default$pop.nprior[["Days"]])
})

test_that("Prior_Setup_GLMM: pop.mu/pop.sd/pop.nprior outputs are always populated and round-trip", {
  dat <- lme4::sleepstudy
  J   <- nlevels(dat$Subject)

  ps_default <- glmbayesCore::Prior_Setup_GLMM(
    Reaction ~ Days + (Days || Subject),
    data = dat
  )

  ## Always present (derived defaults), one named vector per RE component,
  ## matching pop.prior_list's mu / sqrt(diag(Sigma)).
  expect_named(ps_default$pop.mu, c("(Intercept)", "Days"))
  expect_named(ps_default$pop.sd, c("(Intercept)", "Days"))
  expect_named(ps_default$pop.nprior, c("(Intercept)", "Days"))
  for (k in c("(Intercept)", "Days")) {
    expect_identical(ps_default$pop.mu[[k]], ps_default$pop.prior_list[[k]]$mu)
    sd_k <- sqrt(diag(ps_default$pop.prior_list[[k]]$Sigma, names = FALSE))
    expect_equal(ps_default$pop.sd[[k]], sd_k, ignore_attr = TRUE)
    ## n_i = J*w_i/(1-w_i) with w_i from the low-d adaptive default (0.01
    ## for sleepstudy: both components have p_k < 14).
    expect_equal(
      ps_default$pop.nprior[[k]],
      stats::setNames(J * 0.01 / (1 - 0.01), names(ps_default$pop.mu[[k]])),
      ignore_attr = TRUE
    )
  }

  ## Round-trip: feeding the default pop.mu/pop.sd straight back in as
  ## overrides reproduces the same prior (calling once without
  ## pop.mu/pop.sd to inspect the defaults, then supplying them back in).
  ## Sigma is rebuilt as sd^2 after extracting sd = sqrt(diag(Sigma)), so
  ## expect_equal (not identical) -- MacOS x86_64 can differ by 1 ULP.
  ps_roundtrip <- glmbayesCore::Prior_Setup_GLMM(
    Reaction ~ Days + (Days || Subject),
    data = dat,
    pop.mu = ps_default$pop.mu,
    pop.sd = ps_default$pop.sd
  )
  expect_equal(ps_roundtrip$pop.prior_list, ps_default$pop.prior_list)
  expect_equal(ps_roundtrip$pop.mu, ps_default$pop.mu)
  expect_equal(ps_roundtrip$pop.sd, ps_default$pop.sd)
  expect_equal(ps_roundtrip$pop.nprior, ps_default$pop.nprior)

  ## Same round-trip via pop.nprior instead of pop.sd.
  ps_roundtrip_n <- glmbayesCore::Prior_Setup_GLMM(
    Reaction ~ Days + (Days || Subject),
    data = dat,
    pop.mu     = ps_default$pop.mu,
    pop.nprior = ps_default$pop.nprior
  )
  expect_equal(ps_roundtrip_n$pop.prior_list, ps_default$pop.prior_list)
  expect_identical(ps_roundtrip_n$pop.nprior, ps_default$pop.nprior)

  ## Editing a copy of the default pop.mu before feeding it back in changes
  ## only the edited component/entry.
  mu_edit <- ps_default$pop.mu
  mu_edit[["(Intercept)"]]["(Intercept)"] <- 999
  ps_edited <- glmbayesCore::Prior_Setup_GLMM(
    Reaction ~ Days + (Days || Subject),
    data = dat,
    pop.mu = mu_edit
  )
  expect_equal(unname(ps_edited$pop.mu[["(Intercept)"]]), 999)
  expect_identical(ps_edited$pop.mu[["Days"]], ps_default$pop.mu[["Days"]])
})

## Minimal fixture: (Intercept) Block~2 model has p_k = 14 group-level
## predictors (intercept + n_cov covariates); optional random slope on
## 'x1' adds a second RE component with p_k = 1 (low-d by default).
high_d_prior_setup_fixture <- function(n_groups = 60L, n_cov = 13L, slope = "x1") {
  g <- factor(rep(seq_len(n_groups), each = 5L))
  n <- length(g)
  dat <- data.frame(g = g)
  for (j in seq_len(n_cov)) {
    nm <- paste0("x", j)
    dat[[nm]] <- rep(stats::rnorm(n_groups), each = 5L)
  }
  b0 <- stats::rnorm(1)
  b <- stats::rnorm(n_cov)
  X <- as.matrix(dat[, paste0("x", seq_len(n_cov)), drop = FALSE])
  dat$y <- as.numeric(X %*% b) + b0 +
    stats::rnorm(n_groups)[as.integer(g)] + stats::rnorm(n, sd = 0.5)
  covars <- paste(paste0("x", seq_len(n_cov)), collapse = " + ")
  re_part <- if (is.null(slope)) {
    "(1 | g)"
  } else {
    paste0("(1 + ", slope, " || g)")
  }
  form <- stats::as.formula(paste("y ~", covars, "+", re_part))
  list(data = dat, formula = form)
}

test_that("Prior_Setup_GLMM: dimension-adaptive pop.pwt defaults per RE component", {
  dat <- lme4::sleepstudy

  ps_low <- suppressMessages(
    glmbayesCore::Prior_Setup_GLMM(
      Reaction ~ Days + (Days || Subject),
      data = dat
    )
  )
  expect_true(all(unlist(ps_low$pop.pwt) == 0.01))

  fix <- high_d_prior_setup_fixture()
  ps_high <- suppressMessages(
    glmbayesCore::Prior_Setup_GLMM(fix$formula, data = fix$data)
  )
  expect_equal(unique(unname(ps_high$pop.pwt[["(Intercept)"]])), 0.05)
  expect_equal(unique(unname(ps_high$pop.pwt[["x1"]])), 0.01)

  ps_explicit <- glmbayesCore::Prior_Setup_GLMM(
    fix$formula,
    data = fix$data,
    pop.pwt = 0.01
  )
  expect_true(all(unlist(ps_explicit$pop.pwt) == 0.01))

  ps_custom <- suppressMessages(
    glmbayesCore::Prior_Setup_GLMM(
      fix$formula,
      data = fix$data,
      pop.pwt_default_low  = 0.02,
      pop.pwt_default_high = 0.08
    )
  )
  expect_equal(unique(unname(ps_custom$pop.pwt[["(Intercept)"]])), 0.08)
  expect_equal(unique(unname(ps_custom$pop.pwt[["x1"]])), 0.02)

  for (k in names(ps_high$pop.pwt)) {
    expect_equal(
      ps_high$pop.dispersion.pwt[[k]],
      mean(ps_high$pop.pwt[[k]])
    )
  }

  expect_no_error(
    glmbayesCore::Prior_Setup_GLMM(
      Reaction ~ Days + (Days || Subject),
      data = dat,
      pop.sd = list("(Intercept)" = 50, Days = NULL)
    )
  )
})
