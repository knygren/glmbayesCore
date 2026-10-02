## dGamma_list() — per-group dGamma from Prior_Setup_GLMM (Ex_13b data)
##
## Requires dispformula = ~school_id so group.ing_prior is per-group.
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

  design0 <- model_setup(form_lmer, data = dat)
  dat <- subset(dat, school_id %in% names(design0$groupef.rank)[design0$groupef.rank])
  dat$school_id <- droplevels(dat$school_id)

  ps <- Prior_Setup_GLMM(
    form_lmer,
    data = dat,
    pop.pwt = 0.01,
    dispformula = ~school_id,
    group.max_disp_perc = 0.8,
    group.dispersion.pwt = 0.1,
    group.alpha_target = NULL
  )

  disp_pf <- dGamma_list(ps)

  ## Group names (school_id levels)
  g <- names(disp_pf)
  utils::head(g, 10L)

  ## Print a subset (recommended when J is large)
  print(disp_pf, groups = c("2", "18", "47"))

  ## Variation: rebuild truncation windows; print one group
  disp_pf_tight <- dGamma_list(ps, max_disp_perc_measurement = 0.9)
  print(disp_pf_tight, groups = "2")

  ## Downstream attributes (not shown by print.dGamma_list)
  names(attributes(disp_pf))
  utils::head(attr(disp_pf, "window_diagnostics"), 3L)
}
