data("BikeSharing")

# Center continuous predictors
cont_vars <- c("temp", "atemp", "hum", "windspeed", "hr_sin", "hr_cos", "mon_sin", "mon_cos")
BikeSharing_c <- BikeSharing
BikeSharing_c[cont_vars] <- scale(BikeSharing[cont_vars], center = TRUE, scale = FALSE)


# Formula (all variable model)
form <- cnt ~ part_of_day + quarter + holiday + workingday + weathersit +
  hr_sin + hr_cos + mon_sin + mon_cos + temp + atemp + hum + windspeed


# Formula (Limited variable model)
form2 <- cnt ~ part_of_day + quarter + holiday + workingday + weathersit +
  hr_sin + hr_cos + mon_sin + mon_cos


pct_train  <- 0.01   # 2% of data (~350 obs) for fast demo; use 0.05+ for analysis

# 5% train / 95% test split
set.seed(42)
n <- nrow(BikeSharing_c)
idx_train <- sample(n, size = round(pct_train * n))  # ~869 obs
idx_test  <- setdiff(seq_len(n), idx_train)

Bike_train <- BikeSharing_c[idx_train, ]
Bike_test  <- BikeSharing_c[idx_test, ]

# Design matrices and response
X_train <- model.matrix(form2, data = Bike_train)
X_test  <- model.matrix(form2, data = Bike_test)
y_train <- Bike_train$cnt
y_test  <- Bike_test$cnt
n_train <- length(y_train)
n_test  <- length(y_test)
p       <- ncol(X_train)

# Initial theta and prior for population block
theta <- log(y_train + 0.5)
data_pop <- data.frame(theta = theta, Bike_train)
form_pop <- theta ~ part_of_day + quarter + holiday + workingday + weathersit +
  hr_sin + hr_cos + mon_sin + mon_cos
ps_pop <- Prior_Setup(form_pop, family = gaussian(), data = data_pop)


# Demo-friendly settings (increase for your own analysis)
n_burn     <- 200    # 200 burn-in for demo; 
n_sim      <- 1000    # 1000 draws for analysis



beta_out   <- matrix(0, nrow = n_sim, ncol = p)
sigma_out  <- numeric(n_sim)
theta_out  <- matrix(0, nrow = n_sim, ncol = n_train)

set.seed(123)

burn_time <- system.time({
  for (k in seq_len(n_burn)) {
    out_pop <- rglmb(1, theta, X_train, family = gaussian(),
                     pfamily = dNormal_Gamma(ps_pop$mu, ps_pop$Sigma / ps_pop$dispersion,
                                             ps_pop$shape, ps_pop$rate))
    beta   <- as.vector(out_pop$coefficients[1, ])
    sigma_theta_sq <- out_pop$dispersion[1]
    mu_all <- as.vector(X_train %*% beta)
    
    for (i in seq_len(n_train)) {
      theta[i] <- rglmb(1, y_train[i], matrix(1, 1, 1), family = poisson(),
                        pfamily = dNormal(mu = mu_all[i], Sigma = sigma_theta_sq))$coefficients[1, 1]
    }
  }
})

cat("Burn-in time:\n")
print(burn_time)

sim_time <- system.time({
  for (k in seq_len(n_sim)) {
    out_pop <- rglmb(1, theta, X_train, family = gaussian(),
                     pfamily = dNormal_Gamma(ps_pop$mu, ps_pop$Sigma / ps_pop$dispersion,
                                             ps_pop$shape, ps_pop$rate))
    beta   <- as.vector(out_pop$coefficients[1, ])
    sigma_theta_sq <- out_pop$dispersion[1]
    mu_all <- as.vector(X_train %*% beta)
    
    for (i in seq_len(n_train)) {
      theta[i] <- rglmb(1, y_train[i], matrix(1, 1, 1), family = poisson(),
                        pfamily = dNormal(mu = mu_all[i], Sigma = sigma_theta_sq))$coefficients[1, 1]
    }
    beta_out[k, ]  <- beta
    sigma_out[k]   <- sqrt(sigma_theta_sq)
    theta_out[k, ] <- theta
  }
})

cat("\nMain simulation time:\n")
print(sim_time)



# After the Gibbs loop, bind beta and sigma (exclude theta_out)
beta_names <- colnames(X_train)
mcmc_main <- coda::mcmc(cbind(
  beta_out,
  sigma_theta = sigma_out
))
colnames(mcmc_main) <- c(beta_names, "sigma_theta")

# CODA summary
cat("\nCODA summary (main coefficients + sigma):\n")
print(summary(mcmc_main))

cat("\nEffective sample size:\n")
print(coda::effectiveSize(mcmc_main))

# Optional: Geweke diagnostic, autocorrelation
cat("\nGeweke z-scores (first 5):\n")
print(head(coda::geweke.diag(mcmc_main)$z, 5))
cat("\nLag-1 own autocorrelations:\n")
ac1 <- coda::autocorr(mcmc_main, lag = 1)


# ac1 is array (nvar, nvar, 1); extract correlation matrix and then diagonal
ac1_mat <- drop(ac1)
own_ac1 <- diag(ac1_mat)
names(own_ac1) <- colnames(mcmc_main)
print(round(own_ac1, 4))

#plot(mcmc_main)


# Posterior means
beta_mean  <- colMeans(beta_out)
sigma_mean <- mean(sigma_out)

# Out-of-sample predictions on test set
# Option A: point prediction (conditional mean ignoring random effect)
y_pred_cond <- exp(X_test %*% beta_mean)

# Option B: posterior predictive samples (includes random-effect uncertainty)
n_pred <- 1000
y_pred_samples <- matrix(0, nrow = n_pred, ncol = n_test)
for (s in seq_len(n_pred)) {
  idx_s <- sample(n_sim, 1)
  beta_s  <- beta_out[idx_s, ]
  sigma_s <- sigma_out[idx_s]
  theta_test <- rnorm(n_test, mean = X_test %*% beta_s, sd = sigma_s)
  y_pred_samples[s, ] <- rpois(n_test, lambda = exp(theta_test))
}
y_pred_mean <- colMeans(y_pred_samples)
y_pred_sd   <- apply(y_pred_samples, 2, sd)

# Evaluate (e.g., MAE, RMSE)
mae  <- mean(abs(y_test - y_pred_mean))
rmse <- sqrt(mean((y_test - y_pred_mean)^2))


# Option A: conditional mean (no random effect)
y_pred_cond <- exp(X_test %*% beta_mean)
mae_cond  <- mean(abs(y_test - y_pred_cond))
rmse_cond <- sqrt(mean((y_test - y_pred_cond)^2))

# Option B: posterior predictive mean
mae_pp   <- mean(abs(y_test - y_pred_mean))
rmse_pp  <- sqrt(mean((y_test - y_pred_mean)^2))

cat("Option A (conditional): MAE =", round(mae_cond, 2), " RMSE =", round(rmse_cond, 2), "\n")
cat("Option B (post. pred.): MAE =", round(mae_pp, 2), " RMSE =", round(rmse_pp, 2), "\n")
