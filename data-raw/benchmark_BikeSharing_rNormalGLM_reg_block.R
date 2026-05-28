# Benchmark and cross-validate Block 2 of the BikeSharing Gibbs sampler
# (Chapter 18 / demo("Ex_09_BikeSharingPoisson")) against rNormalGLM_reg_block().
#
# Block 2 in the vignette updates observation-level theta[i] with
#   rglmb(1, y_train[i], matrix(1,1,1), poisson(), dNormal(mu = mu_all[i], Sigma = sigma_theta_sq))
# in a loop over n_train. This script replaces that loop with one call to
# rNormalGLM_reg_block() and compares modes, draws (with matched seeds), and time.
#
# Run from package root:
#   Rscript data-raw/benchmark_BikeSharing_rNormalGLM_reg_block.R
# Optional package root:
#   Rscript data-raw/benchmark_BikeSharing_rNormalGLM_reg_block.R "C:/path/to/glmbayes"
#
# Does not modify inst/extdata/BikeSharing_ch14_gibbs.rds (vignette precompute).

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) {
  normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
} else {
  getwd()
}
owd <- setwd(root)
on.exit(setwd(owd), add = TRUE)

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload, e.g. install.packages('pkgload')")
}

pkgload::load_all(export_all = FALSE)

fmt_hms <- function(secs) {
  secs <- as.numeric(secs)
  if (!is.finite(secs) || secs < 0) secs <- 0
  h <- floor(secs / 3600)
  rem <- secs - h * 3600
  m <- floor(rem / 60)
  s <- rem - m * 60
  sprintf("%d h %d min %.2f s", h, m, s)
}

## --- Same data / design as demo Ex_09 and Chapter-18 -------------------------
data("BikeSharing", package = "glmbayes")

cont_vars <- c(
  "temp", "atemp", "hum", "windspeed",
  "hr_sin", "hr_cos", "mon_sin", "mon_cos"
)
BikeSharing_c <- BikeSharing
BikeSharing_c[cont_vars] <- scale(BikeSharing[cont_vars], center = TRUE, scale = FALSE)

form2 <- cnt ~ part_of_day + quarter + holiday + workingday + weathersit +
  hr_sin + hr_cos + mon_sin + mon_cos

pct_train <- 0.01
set.seed(42)
n <- nrow(BikeSharing_c)
idx_train <- sample(n, size = round(pct_train * n))

Bike_train <- BikeSharing_c[idx_train, ]
X_train <- model.matrix(form2, data = Bike_train)
y_train <- Bike_train$cnt
n_train <- length(y_train)
p <- ncol(X_train)

theta <- log(y_train + 0.5)
data_pop <- data.frame(theta = theta, Bike_train)
form_pop <- theta ~ part_of_day + quarter + holiday + workingday + weathersit +
  hr_sin + hr_cos + mon_sin + mon_cos
ps_pop <- Prior_Setup(form_pop, family = gaussian(), data = data_pop)

x_one <- matrix(1, n_train, 1)
colnames(x_one) <- "(Intercept)"

## --- One population (Block 1) update -> mu_all, sigma_theta_sq --------------
set.seed(123)
out_pop <- rglmb(
  1L, theta, X_train, family = gaussian(),
  pfamily = dNormal_Gamma(
    ps_pop$mu, Sigma_0 = ps_pop$Sigma_0,
    ps_pop$shape, ps_pop$rate
  ),
  use_parallel = FALSE,
  verbose = FALSE
)
beta <- as.vector(out_pop$coefficients[1, ])
sigma_theta_sq <- out_pop$dispersion[1]
mu_all <- as.vector(X_train %*% beta)

message("BikeSharing train n = ", n_train, ", p = ", p)
message("sigma_theta_sq = ", signif(sigma_theta_sq, 4))

## --- Helpers: Block 2 (theta | beta, sigma) ---------------------------------
block2_theta_rglmb_loop <- function(
    mu_all,
    sigma_theta_sq,
    y_train,
    seed = NULL,
    Gridtype = 2L,
    use_parallel = FALSE,
    use_opencl = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  theta <- numeric(n_train)
  for (i in seq_len(n_train)) {
    theta[i] <- rglmb(
      1L,
      y = y_train[i],
      x = matrix(1, 1, 1),
      family = poisson(),
      pfamily = dNormal(mu = mu_all[i], Sigma = sigma_theta_sq),
      Gridtype = Gridtype,
      use_parallel = use_parallel,
      use_opencl = use_opencl,
      verbose = FALSE
    )$coefficients[1, 1]
  }
  theta
}

block2_theta_rNormal_reg_loop <- function(
    mu_all,
    sigma_theta_sq,
    y_train,
    seed = NULL,
    Gridtype = 2L,
    use_parallel = FALSE,
    use_opencl = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  theta <- numeric(n_train)
  fam <- poisson()
  for (i in seq_len(n_train)) {
    pl <- list(
      mu = mu_all[i],
      Sigma = matrix(sigma_theta_sq, 1, 1),
      dispersion = 1,
      ddef = FALSE
    )
    theta[i] <- rNormal_reg(
      1L,
      y = y_train[i],
      x = matrix(1, 1, 1),
      prior_list = pl,
      family = fam,
      Gridtype = Gridtype,
      n_envopt = 1L,
      use_parallel = use_parallel,
      use_opencl = use_opencl,
      verbose = FALSE,
      progbar = FALSE
    )$coefficients[1, 1]
  }
  theta
}

block2_theta_reg_block <- function(
    mu_all,
    sigma_theta_sq,
    y_train,
    seed = NULL,
    Gridtype = 2L,
    use_parallel = FALSE,
    use_opencl = FALSE) {
  prior_lists <- lapply(mu_all, function(m) {
    list(
      mu = m,
      Sigma = matrix(sigma_theta_sq, 1, 1),
      dispersion = 1,
      ddef = FALSE
    )
  })
  if (!is.null(seed)) set.seed(seed)
  out <- rNormalGLM_reg_block(
    n = 1L,
    y = y_train,
    x = x_one,
    block = seq_len(n_train),
    prior_lists = prior_lists,
    family = poisson(),
    Gridtype = Gridtype,
    n_envopt = 1L,
    use_parallel = use_parallel,
    use_opencl = use_opencl,
    verbose = FALSE,
    progbar = FALSE
  )
  as.vector(out$coefficients[, 1])
}

block2_coef_mode_reg_block <- function(
    mu_all,
    sigma_theta_sq,
    y_train,
    Gridtype = 2L,
    use_parallel = FALSE,
    use_opencl = FALSE) {
  prior_lists <- lapply(mu_all, function(m) {
    list(
      mu = m,
      Sigma = matrix(sigma_theta_sq, 1, 1),
      dispersion = 1,
      ddef = FALSE
    )
  })
  out <- rNormalGLM_reg_block(
    n = 1L,
    y = y_train,
    x = x_one,
    block = seq_len(n_train),
    prior_lists = prior_lists,
    family = poisson(),
    Gridtype = Gridtype,
    n_envopt = 1L,
    use_parallel = use_parallel,
    use_opencl = use_opencl,
    verbose = FALSE,
    progbar = FALSE
  )
  as.vector(out$coef.mode[, 1])
}

## --- 1) Posterior mode (coef.mode): deterministic check ---------------------
message("\n=== coef.mode (optim): rNormal_reg loop vs rNormalGLM_reg_block ===")
mode_rnr <- {
  fam <- poisson()
  modes <- numeric(n_train)
  for (i in seq_len(n_train)) {
    pl <- list(
      mu = mu_all[i],
      Sigma = matrix(sigma_theta_sq, 1, 1),
      dispersion = 1,
      ddef = FALSE
    )
    cm <- rNormal_reg(
      1L, y_train[i], matrix(1, 1, 1), pl,
      family = fam, Gridtype = 2L, n_envopt = 1L,
      use_parallel = FALSE, use_opencl = FALSE, verbose = FALSE
    )$coef.mode
    modes[i] <- if (is.matrix(cm)) cm[1, 1] else as.numeric(cm)[1]
  }
  modes
}
mode_blk <- block2_coef_mode_reg_block(
  mu_all, sigma_theta_sq, y_train,
  use_parallel = FALSE, use_opencl = FALSE
)
mode_max_abs_diff <- max(abs(mode_rnr - mode_blk))
message("max |mode_loop - mode_block| = ", signif(mode_max_abs_diff, 6))
if (!isTRUE(all.equal(mode_rnr, mode_blk, tolerance = 1e-5, check.names = FALSE))) {
  warning("coef.mode differs between per-obs rNormal_reg and rNormalGLM_reg_block")
}

## --- 2) Stochastic draws: same seed -----------------------------------------
message("\n=== coefficients (one draw): seed-matched comparison ===")
seed_draw <- 2026L
theta_rglmb <- block2_theta_rglmb_loop(
  mu_all, sigma_theta_sq, y_train, seed = seed_draw,
  use_parallel = FALSE, use_opencl = FALSE
)
theta_rnr <- block2_theta_rNormal_reg_loop(
  mu_all, sigma_theta_sq, y_train, seed = seed_draw + 1L,
  use_parallel = FALSE, use_opencl = FALSE
)
theta_blk <- block2_theta_reg_block(
  mu_all, sigma_theta_sq, y_train, seed = seed_draw,
  use_parallel = FALSE, use_opencl = FALSE
)

diff_rglmb_blk <- max(abs(theta_rglmb - theta_blk))
diff_rnr_blk <- max(abs(theta_rnr - theta_blk))
message("max |theta_rglmb_loop - theta_block|     = ", signif(diff_rglmb_blk, 6))
message("max |theta_rNormal_reg_loop - theta_block| = ", signif(diff_rnr_blk, 6))

if (isTRUE(all.equal(theta_rglmb, theta_blk, tolerance = 1e-8))) {
  message("OK: rglmb loop matches rNormalGLM_reg_block under same seed.")
} else if (isTRUE(all.equal(theta_rglmb, theta_blk, tolerance = 1e-4))) {
  message("OK (tolerance 1e-4): rglmb loop ~ block (minor numerical difference).")
} else {
  message("Note: draws differ — RNG path may differ between loop and fused block; compare distributions below.")
}

## --- 3) Distribution check (replicates) -------------------------------------
message("\n=== distribution check (", 20L, " replicates, different seeds) ===")
n_rep <- 20L
mat_rglmb <- matrix(NA_real_, n_rep, n_train)
mat_blk <- matrix(NA_real_, n_rep, n_train)
for (r in seq_len(n_rep)) {
  s <- 1000L + r
  mat_rglmb[r, ] <- block2_theta_rglmb_loop(
    mu_all, sigma_theta_sq, y_train, seed = s,
    use_parallel = FALSE, use_opencl = FALSE
  )
  mat_blk[r, ] <- block2_theta_reg_block(
    mu_all, sigma_theta_sq, y_train, seed = s,
    use_parallel = FALSE, use_opencl = FALSE
  )
}
mean_rglmb <- colMeans(mat_rglmb)
mean_blk <- colMeans(mat_blk)
sd_rglmb <- apply(mat_rglmb, 2, sd)
sd_blk <- apply(mat_blk, 2, sd)
message("mean abs diff of replicate means: ",
        signif(mean(abs(mean_rglmb - mean_blk)), 6))
message("mean abs diff of replicate SDs:   ",
        signif(mean(abs(sd_rglmb - sd_blk)), 6))

## --- 4) Timing: one Block-2 update ------------------------------------------
message("\n=== timing: one full Block-2 update (n_train = ", n_train, ") ===")
n_time <- 5L
time_rglmb <- numeric(n_time)
time_rnr <- numeric(n_time)
time_blk <- numeric(n_time)

for (t in seq_len(n_time)) {
  time_rglmb[t] <- system.time({
    block2_theta_rglmb_loop(
      mu_all, sigma_theta_sq, y_train, seed = NULL,
      use_parallel = FALSE, use_opencl = FALSE
    )
  })["elapsed"]
  time_rnr[t] <- system.time({
    block2_theta_rNormal_reg_loop(
      mu_all, sigma_theta_sq, y_train, seed = NULL,
      use_parallel = FALSE, use_opencl = FALSE
    )
  })["elapsed"]
  time_blk[t] <- system.time({
    block2_theta_reg_block(
      mu_all, sigma_theta_sq, y_train, seed = NULL,
      use_parallel = FALSE, use_opencl = FALSE
    )
  })["elapsed"]
}

summ_time <- function(x) c(mean = mean(x), median = median(x), min = min(x), max = max(x))
s_rglmb <- summ_time(time_rglmb)
s_rnr <- summ_time(time_rnr)
s_blk <- summ_time(time_blk)

message("rglmb loop (vignette style):")
print(round(s_rglmb, 3))
message("rNormal_reg loop (same math, no rglmb wrapper):")
print(round(s_rnr, 3))
message("rNormalGLM_reg_block (Phase 1 R implementation):")
print(round(s_blk, 3))
message("speedup block vs rglmb (mean times): ",
        signif(s_rglmb["mean"] / s_blk["mean"], 3), "x")
message("speedup block vs rNormal_reg loop (mean times): ",
        signif(s_rnr["mean"] / s_blk["mean"], 3), "x")

## --- Save summary (optional local artifact) ---------------------------------
out_path <- file.path(root, "data-raw", "BikeSharing_block_reg_benchmark.rds")
benchmark <- list(
  n_train = n_train,
  p = p,
  sigma_theta_sq = sigma_theta_sq,
  mode_max_abs_diff = mode_max_abs_diff,
  draw_diff_rglmb_block = diff_rglmb_blk,
  draw_diff_rnr_block = diff_rnr_blk,
  replicate_mean_diff = mean(abs(mean_rglmb - mean_blk)),
  timing_seconds = list(
    rglmb_loop = s_rglmb,
    rNormal_reg_loop = s_rnr,
    rNormalGLM_reg_block = s_blk
  ),
  timestamp = Sys.time()
)
saveRDS(benchmark, out_path)
message("\nWrote summary: ", out_path)

message("\nDone.")
