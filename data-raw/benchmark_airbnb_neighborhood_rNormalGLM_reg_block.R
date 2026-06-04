# Benchmark multivariate neighborhood random effects on bayesrules::airbnb.
#
# Model: reviews ~ rating_c + room_type with a separate l1-dimensional coefficient
# vector per neighborhood (k blocks, l2_b > 1). NOT listing-level scalar intercepts.
#
#   Block 1 (hyper mean): rglmb_update on stacked beta (k * l1), design
#       kronecker(I_l1, 1_k) — samples mu | beta, Sigma_fixed.
#   Block 2 (data): rNormalGLM_reg_block_update — samples beta | y, mu, Sigma_fixed.
#   Sigma is fixed at the Prior_Setup calibration (not updated each iteration).
#
#   Rscript data-raw/benchmark_airbnb_neighborhood_rNormalGLM_reg_block.R
#   Rscript data-raw/benchmark_airbnb_neighborhood_rNormalGLM_reg_block.R quick
#   Rscript data-raw/benchmark_airbnb_neighborhood_rNormalGLM_reg_block.R legacy
#   Rscript ... opencl_random          # Block 2: rNormalGLM_reg_block_update
#   Rscript ... opencl_fixed           # Block 1: rglmb_update (hyper); rarely helps
#   Rscript ... opencl_random opencl_fixed
#
# Default: n_burn = 200, n_sim = 1000; use_opencl_fixed/random = FALSE.
# Append `quick` for airbnb_small + short chain.
# Append `opencl_random` and/or `opencl_fixed` to turn OpenCL on for that block.
#
# See also: benchmark_airbnb_rNormalGLM_reg_block.R (listing-level Ch. 18 style).

args <- commandArgs(trailingOnly = TRUE)
run_legacy <- any(tolower(args) %in% c("legacy", "--legacy", "-l"))
run_quick <- any(tolower(args) %in% c("quick", "--quick", "-q"))
use_opencl_fixed <- FALSE
use_opencl_random <- FALSE
if (any(tolower(args) %in% c("opencl_fixed", "use_opencl_fixed", "ocl_fixed"))) {
  use_opencl_fixed <- TRUE
}
if (any(tolower(args) %in% c("opencl_random", "use_opencl_random", "ocl_random"))) {
  use_opencl_random <- TRUE
}
path_args <- args[!tolower(args) %in% c(
  "legacy", "--legacy", "-l",
  "quick", "--quick", "-q",
  "opencl_fixed", "use_opencl_fixed", "ocl_fixed",
  "opencl_random", "use_opencl_random", "ocl_random"
)]

root <- if (length(path_args) >= 1L) {
  normalizePath(path_args[[1]], winslash = "/", mustWork = TRUE)
} else {
  getwd()
}
owd <- setwd(root)
on.exit(setwd(owd), add = TRUE)

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install pkgload, e.g. install.packages('pkgload')")
}
if (!requireNamespace("bayesrules", quietly = TRUE)) {
  stop("Install bayesrules, e.g. install.packages('bayesrules')")
}

pkgload::load_all(export_all = FALSE)

rglmb_update <- function(y, x, family = gaussian(), pfamily,
                         offset = NULL, weights = 1, Gridtype = 2L,
                         n_envopt = NULL, use_parallel = TRUE,
                         use_opencl = FALSE, verbose = FALSE,
                         mu_from_x = TRUE) {
  if (missing(pfamily)) stop("'pfamily' is required.", call. = FALSE)
  x <- as.matrix(x)
  y <- as.vector(y)
  if (nrow(x) != length(y)) stop("nrow(x) must equal length(y).", call. = FALSE)
  out <- rglmb(n = 1L, y = y, x = x, family = family, pfamily = pfamily,
               offset = offset, weights = weights, Gridtype = as.integer(Gridtype),
               n_envopt = n_envopt, use_parallel = use_parallel,
               use_opencl = use_opencl, verbose = verbose)
  beta  <- as.vector(out$coefficients[1L, ])
  names(beta) <- colnames(out$coefficients)
  disp  <- out$dispersion[1L]
  mu_all <- if (isTRUE(mu_from_x)) as.vector(x %*% beta) else NULL
  list(beta = beta, dispersion = disp, sigma_theta_sq = disp,
       sigma_theta = sqrt(disp), mu_all = mu_all, rglmb = out)
}

summ_time <- function(x) {
  c(mean = mean(x), median = median(x), min = min(x), max = max(x))
}

fmt_hms <- function(secs) {
  secs <- as.numeric(secs)
  if (!is.finite(secs) || secs < 0) {
    secs <- 0
  }
  h <- floor(secs / 3600)
  rem <- secs - h * 3600
  m <- floor(rem / 60)
  s <- rem - m * 60
  sprintf("%d h %d min %.2f s", h, m, s)
}

gibbs_report_interval <- function(n) {
  max(1L, min(50L, as.integer(n %/% 10L)))
}

gibbs_message_first_iter_estimates <- function(tag, sec_per_iter, n_burn, n_sim) {
  message(tag, "After iteration 1 — estimated burn-in: ",
          fmt_hms(sec_per_iter * n_burn),
          "; estimated total simulation (burn-in + main): ",
          fmt_hms(sec_per_iter * (n_burn + n_sim)),
          " (", signif(sec_per_iter, 4), " s/iteration)")
}

n_burn <- if (run_quick) 5L else 200L
n_sim  <- if (run_quick) 10L else 1000L
n_gibbs <- n_burn + n_sim

## --- Data -------------------------------------------------------------------
if (run_quick) {
  data("airbnb_small", package = "bayesrules")
  airbnb_use <- airbnb_small
  message("Quick mode: using bayesrules::airbnb_small")
} else {
  data("airbnb", package = "bayesrules")
  airbnb_use <- airbnb
}

airbnb_dat <- airbnb_use
airbnb_dat$rating_c <- airbnb_dat$rating - mean(airbnb_dat$rating)
airbnb_dat$room_type <- factor(airbnb_dat$room_type)
airbnb_dat <- airbnb_dat[complete.cases(airbnb_dat[, c(
  "reviews", "rating_c", "room_type", "neighborhood"
)]), ]

form_x <- reviews ~ rating_c + room_type
X_full <- model.matrix(form_x, data = airbnb_dat)
y_full <- airbnb_dat$reviews
block_full <- factor(airbnb_dat$neighborhood)
l2 <- length(y_full)
l1 <- ncol(X_full)
k <- nlevels(block_full)
beta_names <- colnames(X_full)
neighborhood_names <- levels(block_full)

block_info <- glmbayes:::normalize_block(block_full, l2)
stopifnot(block_info$k == k)

# --- Identifiability preflight -----------------------------------------------
# Intercept-only hyper (X_nbhd = NULL): Level 2 requires at least one
# full-rank block, which is automatic for any neighbourhood with n >= l1.
id_check <- block_check_identifiability_xy(
  x          = X_full,
  block      = block_full,
  X_nbhd     = NULL,
  on_failure = "stop"
)
stopifnot(id_check$action == "proceed")
message("Identifiability check passed — proceeding with full model (", k, " neighborhoods).")

message("airbnb (neighborhood RE): n = ", l2, ", k = ", k, ", l1 = ", l1)
message("  predictors: ", paste(beta_names, collapse = ", "))
message("listings per neighborhood: min = ", min(block_info$l2_blocks),
        ", median = ", median(block_info$l2_blocks),
        ", max = ", max(block_info$l2_blocks))
message("Block Gibbs: n_burn = ", n_burn, ", n_sim = ", n_sim,
        " (", n_gibbs, " iterations)")
message("Block 2: multivariate coef vector per neighborhood (NOT x_one / listing RE)")
message("use_opencl_fixed (Block 1 hyper / rglmb_update): ", use_opencl_fixed)
message("Main Gibbs run: use_opencl_fixed = FALSE, use_opencl_random = FALSE")
if (requireNamespace("glmbayes", quietly = TRUE)) {
  message("OpenCL available (glmbayesCore): ", glmbayesCore::glmbayesCore_has_opencl())
}
if (!run_legacy) {
  message("Legacy per-neighborhood rNormal_reg loop skipped (append 'legacy').")
}

airbnb_dat$eta_proxy <- log(y_full + 1)
ps_glm <- Prior_Setup(
  eta_proxy ~ rating_c + room_type,
  family = gaussian(),
  data = airbnb_dat
)
prior_template <- list(
  mu = as.numeric(ps_glm$mu),
  Sigma = ps_glm$Sigma,
  dispersion = 1,
  ddef = FALSE
)

fam <- poisson()

## Hyper block design: vec(beta_mat) is k*l1 with columns aligned to coef names
X_hyper <- kronecker(diag(l1), matrix(1, k, 1))
colnames(X_hyper) <- beta_names



init_b2 <- rNormalGLM_reg_block_update(
  y = y_full,
  x = X_full,
  block = block_full,
  prior_list = prior_template,
  family = fam,
  Gridtype = 2L,
  n_envopt = 1L,
  use_parallel = FALSE,
  use_opencl = use_opencl_random,
  verbose = FALSE,
  progbar = FALSE
)
beta_mat <- init_b2$coefficients
hyper_mu <- as.numeric(prior_template$mu)
hyper_Sigma <- prior_template$Sigma

## Hyper block prior: reuse population calibration from reviews ~ X (l1 coefficients)
pfamily_hyper <- dNormal_Gamma(
  ps_glm$mu, Sigma_0 = ps_glm$Sigma_0,
  ps_glm$shape, ps_glm$rate
)

run_neighborhood_gibbs <- function(store = FALSE,
                                   label = "",
                                   use_opencl_fixed = FALSE,
                                   use_opencl_random = FALSE) {
  hyper_mu_out <- if (store) matrix(0, nrow = n_sim, ncol = l1) else NULL
  coef_mean_out <- if (store) {
    a <- array(0, dim = c(n_sim, k, l1))
    dimnames(a) <- list(
      iteration = NULL,
      neighborhood = neighborhood_names,
      coefficient = beta_names
    )
    a
  } else {
    NULL
  }
  beta_loc <- beta_mat
  hyper_mu_loc <- hyper_mu
  hyper_Sigma_loc <- hyper_Sigma
  tag <- if (nzchar(label)) paste0("[", label, "] ") else ""
  t_chain_start <- Sys.time()

  burn_time <- system.time({
    for (iter in seq_len(n_burn)) {
      b2 <- rNormalGLM_reg_block_update(
        y = y_full,
        x = X_full,
        block = block_full,
        prior_list = list(
          mu = hyper_mu_loc,
          Sigma = hyper_Sigma_loc,
          dispersion = 1,
          ddef = FALSE
        ),
        family = fam,
        Gridtype = 2L,
        n_envopt = 1L,
        use_parallel = FALSE,
        use_opencl = use_opencl_random,
        verbose = FALSE,
        progbar = FALSE
      )
      beta_loc <- b2$coefficients
      b1 <- rglmb_update(
        y = as.vector(beta_loc),
        x = X_hyper,
        family = gaussian(),
        pfamily = pfamily_hyper,
        use_parallel = FALSE,
        use_opencl = use_opencl_fixed,
        verbose = FALSE
      )
      hyper_mu_loc <- b1$beta
      if (iter == 1L) {
        sec_per_iter <- as.numeric(difftime(Sys.time(), t_chain_start, units = "secs"))
        gibbs_message_first_iter_estimates(tag, sec_per_iter, n_burn, n_sim)
      }
    }
  })

  t_burn_elapsed <- as.numeric(burn_time["elapsed"])
  sec_per_iter <- t_burn_elapsed / n_burn
  message(tag, "Burn-in complete: ", fmt_hms(t_burn_elapsed),
          " (", n_burn, " iterations, ",
          signif(sec_per_iter, 4), " s/iteration)")
  message(tag, "Estimated time remaining (main phase, ",
          n_sim, " iterations): ",
          fmt_hms(sec_per_iter * n_sim))

  report_every <- gibbs_report_interval(n_sim)
  t_main_start <- Sys.time()
  sim_time <- system.time({
    for (iter in seq_len(n_sim)) {
      b2 <- rNormalGLM_reg_block_update(
        y = y_full,
        x = X_full,
        block = block_full,
        prior_list = list(
          mu = hyper_mu_loc,
          Sigma = hyper_Sigma_loc,
          dispersion = 1,
          ddef = FALSE
        ),
        family = fam,
        Gridtype = 2L,
        n_envopt = 1L,
        use_parallel = FALSE,
        use_opencl = use_opencl_random,
        verbose = FALSE,
        progbar = FALSE
      )
      beta_loc <- b2$coefficients
      b1 <- rglmb_update(
        y = as.vector(beta_loc),
        x = X_hyper,
        family = gaussian(),
        pfamily = pfamily_hyper,
        use_parallel = FALSE,
        use_opencl = use_opencl_fixed,
        verbose = FALSE
      )
      hyper_mu_loc <- b1$beta
      if (store) {
        hyper_mu_out[iter, ] <- hyper_mu_loc
        coef_mean_out[iter, , ] <- beta_loc
      }
      if (iter == 1L) {
        sec_main <- as.numeric(difftime(Sys.time(), t_main_start, units = "secs"))
        message(tag, "After 1 main iteration — estimated main phase: ",
                fmt_hms(sec_main * n_sim),
                "; estimated total simulation: ",
                fmt_hms(t_burn_elapsed + sec_main * n_sim),
                " (", signif(sec_main, 4), " s/iteration in main)")
      }
      if (iter %% report_every == 0L || iter == n_sim) {
        elapsed_main <- as.numeric(difftime(Sys.time(), t_main_start, units = "secs"))
        remaining <- elapsed_main / iter * (n_sim - iter)
        message(tag, "Main: ", iter, "/", n_sim,
                " — elapsed ", fmt_hms(elapsed_main),
                ", ETA ", fmt_hms(remaining))
      }
    }
  })

  list(
    burn_time = burn_time,
    sim_time = sim_time,
    hyper_mu_out = hyper_mu_out,
    coef_block_out = coef_mean_out,
    beta_final = beta_loc,
    hyper_mu_final = hyper_mu_loc,
    hyper_Sigma_final = hyper_Sigma_loc,
    use_opencl_fixed = use_opencl_fixed,
    use_opencl_random = use_opencl_random
  )
}

# =============================================================================
# Neighborhood two-block Gibbs
# =============================================================================

message("\n========== Neighborhood multivariate Block Gibbs ==========")
message("Started: ", format(Sys.time(), usetz = TRUE))

set.seed(123)
gibbs_nbhd <- run_neighborhood_gibbs(
  store = TRUE,
  use_opencl_fixed = FALSE,
  use_opencl_random = FALSE
)

t_burn <- as.numeric(gibbs_nbhd$burn_time["elapsed"])
t_sim  <- as.numeric(gibbs_nbhd$sim_time["elapsed"])
t_total <- t_burn + t_sim

message("\n--- timing: neighborhood Block Gibbs ---")
message("burn-in (", n_burn, "): ", fmt_hms(t_burn), " (", signif(t_burn, 4), " s)")
message("main (", n_sim, "):     ", fmt_hms(t_sim), " (", signif(t_sim, 4), " s)")
message("TOTAL (", n_gibbs, "):   ", fmt_hms(t_total), " (", signif(t_total, 4), " s)")
message("mean seconds per iteration: ", signif(t_total / n_gibbs, 4))

benchmark <- list(
  dataset = if (run_quick) "bayesrules::airbnb_small" else "bayesrules::airbnb",
  model = "reviews ~ rating_c + room_type; random coef vector per neighborhood",
  n = l2,
  k = k,
  l1 = l1,
  l2_blocks_summary = summary(block_info$l2_blocks),
  quick_mode = run_quick,
  use_opencl_fixed = FALSE,
  use_opencl_random = FALSE,
  block1_method = "rglmb_update (hyper pool on vec(beta), design k*l1 x l1)",
  block2_method = "rNormalGLM_reg_block_update (k neighborhoods, l1 coefs each)",
  timing_seconds = list(
    burn_in = t_burn,
    main = t_sim,
    total = t_total,
    per_iteration_mean = t_total / n_gibbs
  ),
  hyper_mu_posterior_mean = colMeans(gibbs_nbhd$hyper_mu_out),
  beta_names = beta_names,
  timestamp = Sys.time()
)

if (requireNamespace("coda", quietly = TRUE) && n_sim >= 20L) {
  mcmc_hyper <- coda::mcmc(gibbs_nbhd$hyper_mu_out)
  colnames(mcmc_hyper) <- beta_names
  message("\n--- CODA summary (hyper mean mu across neighborhoods) ---")
  print(summary(mcmc_hyper))
  benchmark$effective_size_hyper_mu <- as.list(coda::effectiveSize(mcmc_hyper))
}

## Posterior mean neighborhood coefficients (rows = neighborhoods, cols = coefs)
coef_post_mean <- apply(gibbs_nbhd$coef_block_out, c(2L, 3L), mean)
dimnames(coef_post_mean) <- list(neighborhood_names, beta_names)

coef_table <- rbind(
  coef_post_mean,
  Average = colMeans(coef_post_mean)
)
dimnames(coef_table) <- list(
  c(neighborhood_names, "Average"),
  beta_names
)

message("\n--- Posterior mean coefficients by neighborhood ---")
print(as.data.frame(round(coef_table, 4)), digits = 4)

benchmark$neighborhood_names <- neighborhood_names
benchmark$coef_post_mean_by_neighborhood <- coef_post_mean
benchmark$coef_table_with_average_row <- coef_table
benchmark$coef_average_by_neighborhood <- colMeans(coef_post_mean)

# =============================================================================
# LEGACY: rNormal_reg loop per neighborhood (optional)
# =============================================================================

if (run_legacy) {
  message("\n========== LEGACY: rNormal_reg per neighborhood ==========")
  prior_block <- glmbayes:::normalize_prior_for_blocks(
    prior_list = list(
      mu = gibbs_nbhd$hyper_mu_final,
      Sigma = gibbs_nbhd$hyper_Sigma_final,
      dispersion = 1,
      ddef = FALSE
    ),
    prior_lists = NULL,
    block_info = block_info,
    l1 = l1
  )
  n_time <- if (run_quick) 2L else 3L
  time_blk <- numeric(n_time)
  time_rnr <- numeric(n_time)
  for (t in seq_len(n_time)) {
    time_blk[t] <- system.time({
      rNormalGLM_reg_block_update(
        y = y_full,
        x = X_full,
        block = block_full,
        prior_list = list(
          mu = gibbs_nbhd$hyper_mu_final,
          Sigma = gibbs_nbhd$hyper_Sigma_final,
          dispersion = 1,
          ddef = FALSE
        ),
        family = fam,
        Gridtype = 2L,
        n_envopt = 1L,
        use_parallel = FALSE,
        use_opencl = use_opencl_random,
        verbose = FALSE,
        progbar = FALSE
      )
    })["elapsed"]
    time_rnr[t] <- system.time({
      coef_rnr <- matrix(NA_real_, nrow = k, ncol = l1)
      for (b in seq_len(k)) {
        rows_b <- block_info$rows[[b]]
        out_b <- rNormal_reg(
          n = 1L,
          y = y_full[rows_b],
          x = X_full[rows_b, , drop = FALSE],
          prior_list = prior_block[[b]],
          family = fam,
          Gridtype = 2L,
          n_envopt = 1L,
          use_parallel = FALSE,
          use_opencl = use_opencl_random,
          verbose = FALSE,
          progbar = FALSE
        )
        cb <- out_b$coefficients
        coef_rnr[b, ] <- if (is.matrix(cb)) cb[1L, ] else as.numeric(cb)
      }
    })["elapsed"]
  }
  s_blk <- summ_time(time_blk)
  s_rnr <- summ_time(time_rnr)
  message("rNormalGLM_reg_block_update (one step):")
  print(round(s_blk, 3))
  message("rNormal_reg loop:")
  print(round(s_rnr, 3))
  benchmark$legacy <- list(
    timing_seconds = list(
      rNormalGLM_reg_block = s_blk,
      rNormal_reg_loop = s_rnr
    ),
    speedup_block_vs_loop = s_rnr["mean"] / s_blk["mean"]
  )
}

out_path <- file.path(root, "data-raw", "Airbnb_neighborhood_block_reg_benchmark.rds")
saveRDS(benchmark, out_path)
message("\nWrote summary: ", out_path)
message("Done.")
