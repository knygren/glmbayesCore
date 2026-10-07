## rlmerb() — matrix-level Gaussian LMM sampler on lme4::sleepstudy
##
## Reaction time versus centered days of sleep deprivation, with a random
## intercept and slope by Subject.  lmebayesCore requires uncorrelated
## random-effect terms, hence the `||` form.  For the formula interface see
## lmebayes::lmerb(), which wraps this function.

data(sleepstudy, package = "lme4")
dat <- sleepstudy
dat$Days_c <- dat$Days - mean(dat$Days)

form <- Reaction ~ Days_c + (1 + Days_c || Subject)

design <- model_setup(form, data = dat)
ps     <- Prior_Setup_GLMM(form, data = dat, pop.pwt = 0.01)

## Point estimates only: exact joint Gaussian posterior mean, no sampling.
pt <- rlmerb(
  design          = design,
  pfamily_list    = pfamily_list(ps),
  dispprior_list  = list(dispersion = ps$group.dispersion),
  simulate        = FALSE,
  verbose         = FALSE,
  print_icm_table = FALSE
)
pt$popef.mode

## Draws: fixed sigma^2 and known tau^2_k, so the joint posterior is exactly
## multivariate normal and sim_method = "DEFAULT" samples it directly (iid).
set.seed(1)
fit <- rlmerb(
  n               = 20L,
  design          = design,
  pfamily_list    = pfamily_list(ps),
  dispprior_list  = list(dispersion = ps$group.dispersion),
  progbar         = FALSE,
  verbose         = FALSE,
  print_icm_table = FALSE
)

fit$convergence_info$sim_method_used
print_groupef(fit, draws = 1:5, groups = 1:3)
