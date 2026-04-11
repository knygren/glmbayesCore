# Precompute outputs for vignettes/Chapter-11.Rmd (Section 3.3):
#   (1) lmb with dIndependent_Normal_Gamma, n = 10000
#   (2) Two-block Gibbs (1000 burn-in, 10000 stored iterations), set.seed(180)
#
# Run from package root, after the package loads with your current code:
#   Rscript data-raw/make_Chapter11_Dobson_two_sampler_output.R
# Optional path:
#   Rscript data-raw/make_Chapter11_Dobson_two_sampler_output.R "C:/path/to/glmbayes"
#
# Writes: inst/extdata/Chapter11_Dobson_two_sampler.rds

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

## Dobson (1990) plant weight — same as vignette
ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
weight <- c(ctl, trt)

ps <- Prior_Setup(weight ~ group)
x <- ps$x
p <- ncol(x)
y <- ps$y
mu <- ps$mu
V <- ps$Sigma
shape <- ps$shape
rate <- ps$rate

message("Independent Normal-Gamma lmb (n = 10000) ...")
lmb_D9_v3 <- lmb(
  n = 10000L,
  weight ~ group,
  dIndependent_Normal_Gamma(
    ps$mu,
    ps$Sigma,
    shape = ps$shape + p / 2,
    rate = ps$rate,
    max_disp_perc = 0.99,
    disp_lower = NULL,
    disp_upper = NULL
  )
)

message("Two-block Gibbs: burn-in 1000, then 10000 iterations ...")
set.seed(180)
dispersion2 <- ps$dispersion

for (i in seq_len(1000L)) {
  out1 <- rlmb(
    n = 1L, y = y, x = x,
    pfamily = dNormal(mu = mu, Sigma = V, dispersion = dispersion2)
  )
  out2 <- rlmb(
    n = 1L, y = y, x = x,
    pfamily = dGamma(
      shape = shape, rate = rate,
      beta = out1$coefficients[1, ]
    )
  )
  dispersion2 <- out2$dispersion
}

n_sim_gibbs <- 10000L
beta_out <- matrix(0, nrow = n_sim_gibbs, ncol = 2L)
disp_out <- numeric(n_sim_gibbs)
coef_names <- colnames(out1$coefficients)

for (i in seq_len(n_sim_gibbs)) {
  out1 <- rlmb(
    n = 1L, y = y, x = x,
    pfamily = dNormal(mu = mu, Sigma = V, dispersion = dispersion2)
  )
  out2 <- rlmb(
    n = 1L, y = y, x = x,
    pfamily = dGamma(
      shape = shape, rate = rate,
      beta = out1$coefficients[1, ]
    )
  )
  dispersion2 <- out2$dispersion
  beta_out[i, ] <- out1$coefficients[1, seq_len(2L)]
  disp_out[i] <- out2$dispersion
}

colnames(beta_out) <- coef_names

out <- list(
  indep_norm_gamma = list(
    n_draws = nrow(lmb_D9_v3$coefficients),
    coefficients = lmb_D9_v3$coefficients,
    dispersion = lmb_D9_v3$dispersion,
    coef_colnames = colnames(lmb_D9_v3$coefficients)
  ),
  gibbs_two_block = list(
    seed = 180L,
    n_burn = 1000L,
    n_sim = n_sim_gibbs,
    beta_out = beta_out,
    disp_out = disp_out,
    coef_colnames = colnames(beta_out)
  ),
  prior_setup = list(
    mu = mu,
    Sigma = V,
    shape = shape,
    rate = rate,
    dispersion_ml = ps$dispersion
  ),
  data_digest = list(
    formula = "weight ~ group",
    n_obs = length(y)
  ),
  package_version = as.character(utils::packageVersion("glmbayes")),
  generated_at = Sys.time()
)

dir.create(file.path(root, "inst", "extdata"), recursive = TRUE, showWarnings = FALSE)
dest <- file.path(root, "inst", "extdata", "Chapter11_Dobson_two_sampler.rds")
saveRDS(out, dest, compress = "xz")
message("Wrote ", normalizePath(dest, winslash = "/"))
