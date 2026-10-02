## Prior_Setup_GLMM() — Block~1 / Block~2 prior calibration (Ex_13b data)
##
## Same formula and cleaning as
## demo("Ex_13b_rLMM_estimated_dispersion_known_vcov_BigWordClub_PartVI").
## Requires suggested package bayesrules.

if (requireNamespace("bayesrules", quietly = TRUE)) {
  data(big_word_club, package = "bayesrules")

  dat <- big_word_club
  dat$school_id <- factor(dat$school_id)
  dat <- subset(
    dat,
    !is.na(score_ppvt) &
      !is.na(invalid_ppvt) & invalid_ppvt == 0L &
      complete.cases(dat[, c(
        "score_ppvt", "distracted_a1", "distracted_ppvt",
        "private_school", "title1", "free_reduced_lunch", "school_id"
      )])
  )

  form_lmer <- score_ppvt ~
    private_school + title1 + free_reduced_lunch +
    distracted_ppvt + distracted_a1 +
    free_reduced_lunch:distracted_a1 +
    (1 + distracted_ppvt + distracted_a1 || school_id)

  ## Keep full-rank schools (sampler / per-group ING path)
  design0 <- model_setup(form_lmer, data = dat)
  dat <- subset(dat, school_id %in% names(design0$groupef.rank)[design0$groupef.rank])
  dat$school_id <- droplevels(dat$school_id)

  ## Pooled measurement dispersion (dispformula = ~1, default)
  ps_pooled <- Prior_Setup_GLMM(
    form_lmer,
    data = dat,
    pop.pwt = 0.01
  )
  print(ps_pooled)
  names(ps_pooled$pop.prior_list)

  ## Per-group Block~1 calibration (Ex_13b-style dispformula).
  ## group.alpha_target = NULL skips the ellipsoid pwt search (faster for examples);
  ## Ex_13b uses group.alpha_target = 0.01 with group.dispersion.pwt = 0.1 as floor.
  ps_group <- Prior_Setup_GLMM(
    form_lmer,
    data = dat,
    pop.pwt = 0.01,
    dispformula = ~school_id,
    group.max_disp_perc = 0.8,
    group.dispersion.pwt = 0.1,
    group.alpha_target = NULL
  )
  print(ps_group)
  stopifnot(is.list(ps_group$group.ing_prior[[1L]]))
  length(ps_group$group.ing_prior)
}
