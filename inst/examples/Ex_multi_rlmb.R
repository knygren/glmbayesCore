## Multi-response samplers (iris): shared y and x via multi_prior_setup
## y: four numeric columns; x: intercept + Species (p = 3)

set.seed(42)

n_draw <- 150L
common <- list(
  n = n_draw,
  family = gaussian(),
  use_parallel = FALSE,
  progbar = FALSE
)

ps_multi <- glmbayesCore::multi_prior_setup(
  cbind(Sepal.Length, Sepal.Width, Petal.Length, Petal.Width) ~ Species,
  data = iris,
  family = gaussian()
)

x <- ps_multi[[1L]]$x
y <- do.call(cbind, lapply(ps_multi, function(ps) ps$y))
colnames(y) <- names(ps_multi)

## --- multi_rlmb (pfamily_list) ---------------------------------------------
pfamily_list <- lapply(ps_multi, function(ps) {
  glmbayesCore::dNormal_Gamma(
    mu = ps$mu, Sigma_0 = ps$Sigma_0, shape = ps$shape, rate = ps$rate
  )
})
out_rlmb <- multi_rlmb(
  n = n_draw, y = y, x = x, pfamily_list = pfamily_list,
  use_parallel = FALSE, progbar = FALSE
)
summary(out_rlmb)

## --- multi_rNormal_reg (shared-x path) -------------------------------------
prior_list_normal <- lapply(ps_multi, function(ps) {
  list(mu = as.numeric(ps$mu), Sigma = ps$Sigma, dispersion = ps$dispersion)
})
out_normal <- do.call(
  multi_rNormal_reg, c(common, list(y = y, x = x, prior_list = prior_list_normal))
)
summary(out_normal)

## --- multi_rNormalGamma_reg ------------------------------------------------
prior_list_ng <- lapply(ps_multi, function(ps) {
  list(
    mu = as.numeric(ps$mu), Sigma = ps$Sigma_0,
    shape = ps$shape, rate = ps$rate
  )
})
out_ng <- do.call(
  multi_rNormalGamma_reg, c(common, list(y = y, x = x, prior_list = prior_list_ng))
)
summary(out_ng)

## --- multi_rindepNormalGamma_reg -------------------------------------------
prior_list_ing <- lapply(ps_multi, function(ps) {
  list(
    mu = as.numeric(ps$mu), Sigma = ps$Sigma,
    shape = ps$shape_ING, rate = ps$rate
  )
})
out_ing <- do.call(
  multi_rindepNormalGamma_reg,
  c(common, list(y = y, x = x, prior_list = prior_list_ing))
)
summary(out_ing)
