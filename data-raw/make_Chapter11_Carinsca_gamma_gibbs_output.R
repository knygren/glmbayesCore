# Precompute two-block Gibbs output for vignettes/Chapter-11.Rmd (Section 4.2,
# Carinsca Gamma regression): 1000 burn-in iterations, then 1000 stored draws.
# Matches vignette chunk `Block_Gibbs_gamma_Regression` (kept in the vignette with eval = FALSE).
#
# Run from package root:
#   Rscript data-raw/make_Chapter11_Carinsca_gamma_gibbs_output.R
# Optional path:
#   Rscript data-raw/make_Chapter11_Carinsca_gamma_gibbs_output.R "C:/path/to/glmbayes"
#
# Writes: inst/extdata/Chapter11_Carinsca_gamma_gibbs.rds

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

data(carinsca)
carinsca$Merit <- ordered(carinsca$Merit)
carinsca$Class <- factor(carinsca$Class)
oldopt <- options(contrasts = c("contr.treatment", "contr.treatment"))
on.exit(options(oldopt), add = TRUE)
Claims <- carinsca$Claims
Merit <- carinsca$Merit
Class <- carinsca$Class
Cost <- carinsca$Cost
Claims_Adj <- Claims / 1000

glm.carinsca <- glm(
  Cost / Claims ~ Merit + Class,
  family = Gamma(link = "log"),
  weights = Claims_Adj,
  x = TRUE
)

ps <- Prior_Setup(
  Cost / Claims ~ Merit + Class,
  family = Gamma(link = "log"),
  weights = Claims_Adj
)
mu <- ps$mu
V <- ps$Sigma
shape <- ps$shape
rate <- ps$rate
x <- ps$x
y <- ps$y

dispersion2 <- gamma.dispersion(glm.carinsca)

seed <- 190L
set.seed(seed)

message("Carinsca Gamma two-block Gibbs: burn-in 1000 ...")
suppressWarnings(suppressMessages(
  for (i in seq_len(1000L)) {
    out1 <- rglmb(
      n = 1L, y = y, x = x,
      family = Gamma(link = "log"),
      pfamily = dNormal(mu = mu, Sigma = V, dispersion = dispersion2),
      weights = Claims_Adj
    )
    out2 <- rglmb(
      n = 1L, y = y, x = x,
      family = Gamma(link = "log"),
      pfamily = dGamma(
        shape = shape, rate = rate,
        beta = out1$coefficients[1, ]
      ),
      weights = Claims_Adj
    )
    dispersion2 <- out2$dispersion
  }
))

n_sim <- 1000L
beta_out <- matrix(0, nrow = n_sim, ncol = ncol(x))
disp_out <- numeric(n_sim)
iters_out <- numeric(n_sim)

message("Carinsca Gamma two-block Gibbs: storing ", n_sim, " draws ...")
suppressWarnings(suppressMessages(
  for (i in seq_len(n_sim)) {
    out1 <- rglmb(
      n = 1L, y = y, x = x,
      family = Gamma(link = "log"),
      pfamily = dNormal(mu = mu, Sigma = V, dispersion = dispersion2),
      weights = Claims_Adj
    )
    out2 <- rglmb(
      n = 1L, y = y, x = x,
      family = Gamma(link = "log"),
      pfamily = dGamma(
        shape = shape, rate = rate,
        beta = out1$coefficients[1, ]
      ),
      weights = Claims_Adj
    )
    dispersion2 <- out2$dispersion
    beta_out[i, ] <- out1$coefficients[1, seq_len(ncol(x))]
    disp_out[i] <- out2$dispersion
    iters_out[i] <- out2$iters
  }
))

colnames(beta_out) <- colnames(out1$coefficients)

out <- list(
  gibbs_gamma = list(
    seed = seed,
    n_burn = 1000L,
    n_sim = n_sim,
    beta_out = beta_out,
    disp_out = disp_out,
    iters_out = iters_out,
    coef_colnames = colnames(beta_out)
  ),
  prior_digest = list(
    formula = "Cost/Claims ~ Merit + Class",
    family = "Gamma(log)",
    n_obs = length(y),
    n_coef = ncol(x)
  ),
  package_version = as.character(utils::packageVersion("glmbayes")),
  generated_at = Sys.time()
)

dir.create(file.path(root, "inst", "extdata"), recursive = TRUE, showWarnings = FALSE)
dest <- file.path(root, "inst", "extdata", "Chapter11_Carinsca_gamma_gibbs.rds")
saveRDS(out, dest, compress = "xz")
message("Wrote ", normalizePath(dest, winslash = "/"))
