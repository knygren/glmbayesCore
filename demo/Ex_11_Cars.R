## mtcars: Prior_Setup() - dNormal, dNormal_Gamma, dIndependent_Normal_Gamma
##
## Uses `Prior_Setup()` with **`pwt = 0.001`** (weak prior; override defaults). Prior
## mean **mu** uses the **full model** MLE for intercept
## and effects (`intercept_source` / `effects_source` = `"full_model"`) so the
## prior-mean penalty in `S_marg` is small and posterior variance is easier to
## compare to **(1-pwt)*vcov(lm)** (scalar-pwt Zellner shrinkage of the likelihood
## covariance; see Chapter A12).
##
## `dNormal()` uses `ps$Sigma` and **default** `dispersion = ps$dispersion`. 
## `dIndependent_Normal_Gamma()` uses `shape = ps$shape + p/2` with `p = ncol(ps$x)`
## (same `ps$rate` as `dNormal_Gamma()`); see `?Prior_Setup`.
##
## Run: demo(Ex_11_Cars, package = "glmbayes")

library(glmbayes)

vcov_congruence <- function(V_post, V_base) {
  U <- chol(V_base)
  Ui <- backsolve(U, diag(nrow(U)))
  M <- crossprod(Ui, V_post %*% Ui)
  ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  list(M = M, eigenvalues = ev)
}

congruence_line <- function(V_post, label, V_base, base_label) {
  cq <- vcov_congruence(V_post, V_base)
  ev <- cq$eigenvalues
  sprintf(
    "%s vs %s: eigenvalues min / mean / max = %s / %s / %s",
    label,
    base_label,
    format(min(ev), digits = 6),
    format(mean(ev), digits = 6),
    format(max(ev), digits = 6)
  )
}

data("mtcars", package = "datasets")
mt <- mtcars
mt$c_wt  <- as.numeric(scale(mtcars$wt, center = TRUE, scale = FALSE))
mt$c_cyl <- as.numeric(scale(mtcars$cyl, center = TRUE, scale = FALSE))

form <- mpg ~ c_wt + c_cyl
fit_lm <- lm(form, data = mt, x = TRUE, y = TRUE)
V_lm <- vcov(fit_lm)

pwt=0.001

# Explicit full-model prior means (must match lmb() below - same ps)
ps <- Prior_Setup(
  form,
  gaussian(),
  data = mt,
  pwt = pwt,
  intercept_source = "full_model",
  effects_source = "full_model"
)
p <- ncol(ps$x)
pwt <- ps$PriorSettings$pwt
V_lm_shrunk <- (1 - pwt) * V_lm

cat("\n======== Prior_Setup (pwt = 0.001, Gaussian) =========\n")
print(ps$PriorSettings)
cat("p =", p, "  ps$shape =", ps$shape, "  ps$shape + p/2 (ING) =", ps$shape + p / 2, "\n")

n_mc <- 100000L

fit_dn <- lmb(
  form,
  data = mt,
  pfamily = dNormal(
    mu = ps$mu,
    Sigma = ps$Sigma,
    dispersion = ps$dispersion
  ),
  n = n_mc
)

fit_ng <- lmb(
  form,
  data = mt,
  pfamily = dNormal_Gamma(
    ps$mu,
    Sigma_0 = ps$Sigma_0,
    shape = ps$shape,
    rate = ps$rate
  ),
  n = n_mc
)

fit_ing <- lmb(
  form,
  data = mt,
  pfamily = dIndependent_Normal_Gamma(
    ps$mu,
    ps$Sigma,
    shape = ps$shape + p / 2,
    rate = ps$rate
  ),
  n = n_mc
)


cat("\n======== Classical vs. Posterior vcov(lmb) =========\n")
cat("Classical Scaled:\n")
print((1 - pwt) * vcov(fit_lm))
cat("\ndNormal (dispersion = ps$dispersion):\n")
print(vcov(fit_dn))
cat("\ndNormal_Gamma:\n")
print(vcov(fit_ng))
cat("\ndIndependent_Normal_Gamma:\n")
print(vcov(fit_ing))

cat("\n======== Congruence vs (1-pwt)*vcov(lm) =========\n")
cat("pwt =", pwt, "\n")
cat(congruence_line(vcov(fit_dn), "dNormal", V_lm_shrunk, "(1-pwt)*vcov(lm)"), "\n", sep = "")
cat(congruence_line(vcov(fit_ng), "NG", V_lm_shrunk, "(1-pwt)*vcov(lm)"), "\n", sep = "")
cat(congruence_line(vcov(fit_ing), "ING", V_lm_shrunk, "(1-pwt)*vcov(lm)"), "\n", sep = "")

cat("\n======== Dispersion (classical vs posterior means) =========\n")
disp_classical <- summary(fit_lm)$sigma^2
cat("lm sigma^2 (RSS/(n-p))   =", disp_classical, "\n")
cat("ps$dispersion (dNormal)  =", ps$dispersion, "\n")
cat("mean(lmb dNormal dispersion) =", mean(fit_dn$dispersion), "\n")
cat("mean(lmb NG dispersion)  =", mean(fit_ng$dispersion), "\n")
cat("mean(lmb ING dispersion) =", mean(fit_ing$dispersion), "\n")

invisible(NULL)
