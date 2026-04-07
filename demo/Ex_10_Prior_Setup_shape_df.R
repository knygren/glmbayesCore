## Prior_Setup(shape_df) demo: matching lm()-style variance (n - p) vs n-counting
##
## Prior means = full-model MLE (intercept_source / effects_source = "full_model")
## so comparisons focus on variance, not location.  n = 20 PlantGrowth-style rows.
##
## Idea (weak prior on coefficients, see ?Prior_Setup and shape_df):
##   * lm()-like (residual df n - p): use dNormal_Gamma with default shape_df
##     (shape ~ n_prior/2) together with dIndependent_Normal_Gamma and shape_df
##     "n_prior+p" (and the same rate rule as in the package).
##   * Variance analogous to counting n (not n - p): use dIndependent_Normal_Gamma
##     with default shape_df together with dNormal_Gamma and shape_df "n_prior-p".
##
## Calibrations used here:
##   * For the n-p style (section 1): pwt = 0.01 for both Prior_Setup calls; do not
##     pass n_prior (it is implied by pwt).  First call: default shape_df for NG;
##     second: shape_df = "n_prior+p" for ING.
##     pwt = 0.01 is not a flat prior: implied n_prior = (pwt/(1-pwt))*n_effective
##     (about 0.2 when n = 20), plus a Gamma prior on residual precision.  Posterior
##     variance can sit slightly *below* vcov(lm) (ratio < 1 on the diagonal), not
##     exactly equal; that is expected unless pwt -> 0.
##   * For the n style (section 2): n_prior_sec2 = 3 for both Prior_Setup calls
##     (requires n_prior > p for "n_prior-p").  Section 2 replicates each row 5 times
##     (n_effective = 100).  Gamma *rates* for lmb: ING has shape = n_prior/2 = 1.5,
##     so rate = (shape - 1) * d_np with d_np = summary(lm)$sigma^2 (RSS/(n-p)).
##     NG with shape_df "n_prior-p" has shape = (3-2)/2 = 0.5 < 1, so E[sigma^2]
##     calibration does not apply; section 2 uses Prior_Setup default rate for NG.
##
## MC draws: n_mc = 100000 for stable vcov(lmb).  Early output includes a "scaled
## vcov(lm)" = vcov(lm) * (1 - pwt): heuristic when posterior beta precision scales
## like 1/(1-pwt) under the Zellner setup (section 1: pwt = 0.01; section 2: pwt from
## n_prior_sec2 and n_effective = 100 after 5-fold replication).  Congruence summaries use
## M = U^{-T} vcov(lmb) U^{-1} with vcov(lm) = U'U (chol); eigenvalues of M replace
## elementwise ratios for a matrix-consistent scale comparison.
##
## Run: demo(Ex_10_Prior_Setup_shape_df, package = "glmbayes")

library(glmbayes)

## Matrix comparison of covariances: V_base = crossprod(U) with U = chol(V_base)
## upper-triangular.  M = U^{-T} V_post U^{-1} = t(Ui) %*% V_post %*% Ui with
## Ui = U^{-1}.  Eigenvalues of M are generalized eigenvalues of (V_post, V_base);
## if V_post = k * V_base then M = k * I (all eigenvalues k).
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

ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
weight <- c(ctl, trt)

## Section 2 only: replicate each observation 5 times (n = 100) for smaller implied pwt.
rep_sec2 <- 5L
idx_s2 <- rep(seq_along(weight), each = rep_sec2)
weight_s2 <- weight[idx_s2]
group_s2 <- group[idx_s2]
n_eff_s2 <- length(weight_s2)
n_prior_sec2 <- 3L

n_mc <- 100000L

fit_lm <- lm(weight ~ group, x = TRUE, y = TRUE)
V_lm <- vcov(fit_lm)

fit_lm_s2 <- lm(weight_s2 ~ group_s2, x = TRUE, y = TRUE)
V_lm_s2 <- vcov(fit_lm_s2)

pwt_sec1 <- 0.01
V_lm_scaled_s1 <- V_lm * (1 - pwt_sec1)
pwt_sec2 <- n_prior_sec2 / (n_prior_sec2 + n_eff_s2)
V_lm_scaled_s2 <- V_lm_s2 * (1 - pwt_sec2)

cat("\n======== Reference: vcov(lm) and scaled vcov (Zellner variance heuristic) =========\n")
cat(
  "Heuristic: if beta precision scales by 1/(1-pwt), then Var ~ (1-pwt)*vcov(lm).\n",
  "So scaled vcov(lm) = vcov(lm) * (1-pwt).  Compare posterior vcov(lmb) to these.\n\n",
  sep = ""
)
cat("--- Ordinary least squares ---\n")
print(stats::coef(fit_lm))
cat("\nvcov(lm):\n")
print(V_lm)
cat(
  "\n--- Section 1 (pwt = ", pwt_sec1, "): scaled = vcov(lm) * (1-pwt) = vcov(lm) * ",
  format(1 - pwt_sec1, digits = 6), " ---\n",
  sep = ""
)
print(V_lm_scaled_s1)
cat(
  "\n--- Section 2 (5x replication: n_effective = ", n_eff_s2,
  "; n_prior = ", n_prior_sec2, " => pwt = ", format(pwt_sec2, digits = 6),
  "): lm on replicated data; scaled = vcov(lm_s2) * (1-pwt) = vcov(lm_s2) * ",
  format(1 - pwt_sec2, digits = 6), " ---\n",
  sep = ""
)
cat("\nvcov(lm) on replicated data (same coef as n=20, smaller SE^2):\n")
print(V_lm_s2)
print(V_lm_scaled_s2)

## ----- Section 1: lm-like (n - p) variance -----
## Prior_Setup twice: pwt = 0.01, no n_prior; NG uses default shape_df, ING uses n_prior+p.

cat("\n")
cat("######################################################################\n")
cat("# Section 1: priors aimed at lm()-style (n - p) residual counting   #\n")
cat("#   dNormal_Gamma: default shape_df (n_prior)                       #\n")
cat("#   dIndependent_Normal_Gamma: shape_df = n_prior+p, pwt = 0.01    #\n")
cat("######################################################################\n")

p1_ng <- Prior_Setup(
  weight ~ group,
  pwt = 0.01,
  shape_df = "n_prior",
  intercept_source = "full_model",
  effects_source = "full_model"
)

p1_ing <- Prior_Setup(
  weight ~ group,
  pwt = 0.01,
  shape_df = "n_prior+p",
  intercept_source = "full_model",
  effects_source = "full_model"
)

fit_p1_ng <- lmb(
  weight ~ group,
  pfamily = dNormal_Gamma(
    p1_ng$mu,
    p1_ng$Sigma / p1_ng$dispersion,
    shape = p1_ng$shape,
    rate = p1_ng$rate
  ),
  n = n_mc
)

fit_p1_ing <- lmb(
  weight ~ group,
  pfamily = dIndependent_Normal_Gamma(
    p1_ing$mu,
    p1_ing$Sigma,
    shape = p1_ing$shape,
    rate = p1_ing$rate
  ),
  n = n_mc
)

cat("\n======== Section 1: location (lm vs prior means) =========\n")
print(stats::coef(fit_lm))
print(cbind(lm = coef(fit_lm), ng_mu = p1_ng$mu[, 1], ing_mu = p1_ing$mu[, 1]))

cat("\n======== Section 1: posterior vcov(lmb), dNormal_Gamma (default shape_df) =========\n")
print(vcov(fit_p1_ng))

cat("\n======== Section 1: posterior vcov(lmb), dIndependent_Normal_Gamma (n_prior+p) =========\n")
print(vcov(fit_p1_ing))

cat("\n======== Section 1: ratios vcov(lmb) / vcov(lm) =========\n")
R1_ng <- vcov(fit_p1_ng) / V_lm
R1_ing <- vcov(fit_p1_ing) / V_lm
cat("dNormal_Gamma (default n_prior):\n")
print(R1_ng)
cat("dIndependent_Normal_Gamma (n_prior+p):\n")
print(R1_ing)
cat(sprintf(
  "Mean diagonal ratio vs lm:  NG (default) %5.3f;  ING (n_prior+p) %5.3f\n",
  mean(diag(R1_ng)),
  mean(diag(R1_ing))
))
R1_ng_s <- vcov(fit_p1_ng) / V_lm_scaled_s1
R1_ing_s <- vcov(fit_p1_ing) / V_lm_scaled_s1
cat("\nSection 1: vcov(lmb) / scaled vcov(lm)  (near 1 if Zellner scalar heuristic fits):\n")
cat("dNormal_Gamma (default n_prior):\n")
print(R1_ng_s)
cat("dIndependent_Normal_Gamma (n_prior+p):\n")
print(R1_ing_s)
cat(sprintf(
  "Mean diagonal vs scaled:  NG %5.3f;  ING (n_prior+p) %5.3f\n",
  mean(diag(R1_ng_s)),
  mean(diag(R1_ing_s))
))
cat(
  "\nNote: ratios vs raw lm slightly below 1 are normal (finite pwt, prior on tau). ",
  "ING + n_prior+p often deviates more from the pure (1-pwt) variance story than NG.\n",
  sep = ""
)

cat("\nSection 1: congruence M = U^{-T} vcov(lmb) U^{-1}  (V_lm = U'U from chol):\n")
cat(congruence_line(vcov(fit_p1_ng), "NG", V_lm, "vcov(lm)"), "\n", sep = "")
cat(congruence_line(vcov(fit_p1_ing), "ING (n_prior+p)", V_lm, "vcov(lm)"), "\n", sep = "")
cat(
  "  (If V_post = k * V_lm, all eigenvalues equal k; vs scaled target k = (1-pwt), compare to NG.)\n",
  sep = ""
)
cat(congruence_line(vcov(fit_p1_ng), "NG", V_lm_scaled_s1, "scaled vcov(lm)"), "\n", sep = "")
cat(
  congruence_line(vcov(fit_p1_ing), "ING (n_prior+p)", V_lm_scaled_s1, "scaled vcov(lm)"),
  "\n",
  sep = ""
)

## ----- Section 2: n-counting style variance -----
## Prior_Setup twice: n_prior_sec2; ING default shape_df, NG n_prior-p.
## Gamma rates for lmb: ING E[sigma^2] prior = rate/(shape-1) with d_np = RSS/(n-p);
## NG uses default rate when shape <= 1.

cat("\n")
cat("######################################################################\n")
cat("# Section 2: priors aimed at n-style (not n - p) counting          #\n")
cat("#   Data: each plant weight repeated ", rep_sec2, " times (n = ", n_eff_s2, ")   #\n", sep = "")
cat("#   n_prior = ", n_prior_sec2, "; ING: default shape_df; NG: n_prior-p           #\n", sep = "")
cat("#   ING: Gamma rate = (shape - 1) * d_np; NG: shape < 1 => default Prior_Setup rate. #\n")
cat("######################################################################\n")

p2_ing <- Prior_Setup(
  weight_s2 ~ group_s2,
  n_prior = n_prior_sec2,
  shape_df = "n_prior",
  intercept_source = "full_model",
  effects_source = "full_model"
)

p2_ng <- Prior_Setup(
  weight_s2 ~ group_s2,
  n_prior = n_prior_sec2,
  shape_df = "n_prior-p",
  intercept_source = "full_model",
  effects_source = "full_model"
)

## Prior on precision tau ~ Gamma(shape, rate) (glmbayes parameterization).
## sigma^2 = 1/tau has prior mean E[sigma^2] = rate / (shape - 1) for shape > 1.
## Classical unbiased residual variance: d_np = RSS/(n-p) = summary(lm)$sigma^2.
## Choose rate so E[sigma^2]_prior = d_np  =>  rate = (shape - 1) * d_np.
## Different shape_df => different shape => different rate (Prior_Setup $rate ignored here).
d_np_s2 <- as.numeric(summary(fit_lm_s2)$sigma^2)
shape_p2_ing <- p2_ing$shape
shape_p2_ng <- p2_ng$shape
rate_p2_ing <- (shape_p2_ing - 1) * d_np_s2
rate_p2_ng <- if (shape_p2_ng > 1) {
  (shape_p2_ng - 1) * d_np_s2
} else {
  p2_ng$rate
}

cat("\n======== Section 2: Gamma rate calibration (ING); NG default if shape <= 1 =========\n")
cat("d_np = summary(lm_s2)$sigma^2 (RSS/(n-p)) =", format(d_np_s2, digits = 8), "\n")
cat("ING: shape =", format(shape_p2_ing, digits = 6), "  rate = (shape-1)*d_np =", format(rate_p2_ing, digits = 8), "\n")
cat("  (Prior_Setup rate was", format(p2_ing$rate, digits = 8), ")\n")
if (shape_p2_ng > 1) {
  cat("NG:  shape =", format(shape_p2_ng, digits = 6), "  rate = (shape-1)*d_np =", format(rate_p2_ng, digits = 8), "\n")
} else {
  cat(
    "NG:  shape =", format(shape_p2_ng, digits = 6),
    "  (<= 1: no E[sigma^2] calibration)  rate = Prior_Setup default =", format(rate_p2_ng, digits = 8), "\n"
  )
}
cat("  (Prior_Setup rate was", format(p2_ng$rate, digits = 8), ")\n")

fit_p2_ing <- lmb(
  weight_s2 ~ group_s2,
  pfamily = dIndependent_Normal_Gamma(
    p2_ing$mu,
    p2_ing$Sigma,
    shape = shape_p2_ing,
    rate = rate_p2_ing
  ),
  n = n_mc
)

fit_p2_ng <- lmb(
  weight_s2 ~ group_s2,
  pfamily = dNormal_Gamma(
    p2_ng$mu,
    p2_ng$Sigma / p2_ng$dispersion,
    shape = shape_p2_ng,
    rate = rate_p2_ng
  ),
  n = n_mc
)

cat("\n======== Section 2: location (lm vs prior means, replicated n) =========\n")
print(stats::coef(fit_lm_s2))
print(cbind(lm = coef(fit_lm_s2), ing_mu = p2_ing$mu[, 1], ng_mu = p2_ng$mu[, 1]))

cat("\n======== Section 2: reference vcov(lm) on replicated data =========\n")
print(V_lm_s2)

cat("\n======== Section 2: posterior vcov(lmb), dIndependent_Normal_Gamma (default shape_df) =========\n")
print(vcov(fit_p2_ing))

cat("\n======== Section 2: posterior vcov(lmb), dNormal_Gamma (n_prior-p) =========\n")
print(vcov(fit_p2_ng))

cat("\n======== Section 2: ratios vcov(lmb) / vcov(lm) on replicated data =========\n")
R2_ing <- vcov(fit_p2_ing) / V_lm_s2
R2_ng <- vcov(fit_p2_ng) / V_lm_s2
cat("dIndependent_Normal_Gamma (default shape_df, calibrated rate):\n")
print(R2_ing)
cat("dNormal_Gamma (n_prior-p; Prior_Setup default rate, shape = 0.5):\n")
print(R2_ng)
cat(sprintf(
  "Mean diagonal ratio vs lm:  ING (default) %5.3f;  NG (n_prior-p) %5.3f\n",
  mean(diag(R2_ing)),
  mean(diag(R2_ng))
))
R2_ing_s <- vcov(fit_p2_ing) / V_lm_scaled_s2
R2_ng_s <- vcov(fit_p2_ng) / V_lm_scaled_s2
cat("\nSection 2: vcov(lmb) / scaled vcov(lm_s2)  (pwt = n_prior/(n_prior+n_effective)):\n")
cat("dIndependent_Normal_Gamma (default n_prior):\n")
print(R2_ing_s)
cat("dNormal_Gamma (n_prior-p; default rate):\n")
print(R2_ng_s)
cat(sprintf(
  "Mean diagonal vs scaled:  ING %5.3f;  NG (n_prior-p) %5.3f\n",
  mean(diag(R2_ing_s)),
  mean(diag(R2_ng_s))
))

cat("\nSection 2: congruence M = U^{-T} vcov(lmb) U^{-1}  (base = vcov(lm_s2)):\n")
cat(congruence_line(vcov(fit_p2_ing), "ING (default)", V_lm_s2, "vcov(lm_s2)"), "\n", sep = "")
cat(congruence_line(vcov(fit_p2_ng), "NG (n_prior-p)", V_lm_s2, "vcov(lm_s2)"), "\n", sep = "")
cat(congruence_line(vcov(fit_p2_ing), "ING (default)", V_lm_scaled_s2, "scaled vcov(lm_s2)"), "\n", sep = "")
cat(congruence_line(vcov(fit_p2_ng), "NG (n_prior-p)", V_lm_scaled_s2, "scaled vcov(lm_s2)"), "\n", sep = "")

cat("\nSee ?Prior_Setup, argument shape_df.\n")

invisible(NULL)
