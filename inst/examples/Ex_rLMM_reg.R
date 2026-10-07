## rLMMNormal_reg_known_vcov_iid() — exact iid draws on big_word_club
##
## Same data / formula path as the Ex_10 / Ex_13b demos; uses only the
## fixed-dispersion, known-vcov iid sampler. Requires suggested package
## bayesrules.

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

  design <- model_setup(form_lmer, data = dat)
  ps <- Prior_Setup_GLMM(form_lmer, data = dat, pop.pwt = 0.01)
  pf <- pfamily_list(ps)

  ## Exact iid draws (fixed sigma^2, known tau^2_k / all dNormal).
  ## design$group already carries attr(*, "group_name").
  set.seed(1)
  fit <- rLMMNormal_reg_known_vcov_iid(
    n = 20L,
    y = design$y,
    D = design$D,
    group = design$group,
    W = design$W,
    pfamily_list = pf,
    dispprior_list = list(dispersion = ps$group.dispersion),
    progbar = FALSE,
    verbose = FALSE
  )

  ## Exact-iid route tag, then the first few draws for the first few groups.
  fit$convergence_info$draw_engine
  fit$convergence_info$sim_method_used
  print_groupef(fit, draws = 1:5, groups = 1:3)
}
