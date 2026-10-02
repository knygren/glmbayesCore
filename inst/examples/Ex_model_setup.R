## model_setup() — design + reference lmer (Ex_13b / big_word_club)
##
## Same formula and data cleaning as
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

  ## Default: design matrices + reference lmer + identifiability
  design <- model_setup(form_lmer, data = dat)
  print(design)

  ## Full-rank filter (same as Ex_13b before prior setup / sampling)
  full_rank <- names(design$groupef.rank)[design$groupef.rank]
  dat_fr <- subset(dat, school_id %in% full_rank)
  dat_fr$school_id <- droplevels(dat_fr$school_id)
  design_fr <- model_setup(form_lmer, data = dat_fr)
  stopifnot(all(design_fr$groupef.rank))
  print(design_fr)

  ## Variation: design only (no reference merMod fit)
  design_only <- model_setup(form_lmer, data = dat_fr, fit_mer = FALSE)
  stopifnot(is.null(design_only$lmer))
  names(design_only$W)
}
