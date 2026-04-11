library(MASS)

# Original dataset
data("Boston")

# Identify predictors (exclude response)
predictors <- setdiff(names(Boston), "medv")

# Mean-center predictors (subtract column means, keep original units)
Boston_centered <- Boston
Boston_centered[predictors] <- scale(Boston[predictors], center = TRUE, scale = FALSE)

# Quick check: column means should be ~0
colMeans(Boston_centered[predictors])

# --------------------------------
# Explicit regression formula with all predictors
# --------------------------------
form <- medv ~ 
  crim + zn + 
  indus + chas + nox + age +  dis + rad + tax + ptratio + black + lstat+ rm 

# --------------------------------
# Classical linear model
# --------------------------------
lm.boston <- lm(form, data = Boston_centered, x = TRUE, y = TRUE)
summary(lm.boston)

readline("\n--- End of Classical Model block. Press <Enter> to continue... ---")

# --------------------------------
# Conjugate Normal prior (fixed dispersion)
# --------------------------------
ps <- Prior_Setup(form, gaussian(), data = Boston_centered)
p <- ncol(ps$x)

lmb.boston <- lmb(
  form,
  data     = Boston_centered,
  pfamily  = dNormal(mu = ps$mu,
                     Sigma = ps$Sigma,
                     dispersion = ps$dispersion)
)


summary(lmb.boston)

readline("\n--- End of Normal Prior sampler block. Press <Enter> to continue... ---")


# --------------------------------
# Conjugate Normal-Gamma prior
# --------------------------------
lmb.boston_v2 <- lmb(
  form,
  data     = Boston_centered,
  pfamily  = dNormal_Gamma(mu     = ps$mu,
                           Sigma_0 = ps$Sigma_0,
                           shape  = ps$shape,
                           rate   = ps$rate)
)


summary(lmb.boston_v2)

readline("\n--- End of Normal-Gamma Prior sampler block. Press <Enter> to continue... ---")

# --------------------------------
# Independent Normal-Gamma prior (non-conjugate, uses envelope)
# --------------------------------
lmb.boston_v3 <- glmb(n = 1000,
                      form,
                      data       = Boston_centered,
                      family     = gaussian(),
                      pfamily    = dIndependent_Normal_Gamma(ps$mu, ps$Sigma,
                                                             shape = ps$shape + p / 2,
                                                             rate  = ps$rate
                      ),
                      use_parallel = TRUE,
                      use_opencl = TRUE
)


summary(lmb.boston_v3)


