############################### Cleveland OpenCL example (rglmb) ####################

data("Cleveland")
head(Cleveland)
summary(Cleveland)

## OpenCL-accelerated Bayesian logistic regression via rglmb.
## Runs only when this build was compiled with OpenCL support.
\donttest{
if (glmbayesCore_has_opencl()) {
  form <- hd ~ age + sex + cp + trestbps + chol +
    fbs + restecg + thalach + exang + oldpeak + slope + ca + thal

  ps <- Prior_Setup(
    form,
    family = binomial(logit),
    data = Cleveland
  )

  ## Prior_Setup keeps the factor response; rglmb needs numeric 0/1
  ## (same coding as glm: second factor level is success).
  y <- ps$y
  if (is.factor(y)) {
    y <- as.numeric(y) - 1L
  } else {
    y <- as.numeric(y)
  }

  fit <- rglmb(
    n = 1000,
    y = y,
    x = as.matrix(ps$x),
    family = binomial(link = "logit"),
    pfamily = dNormal(mu = ps$mu, Sigma = ps$Sigma),
    weights = if (!is.null(ps$weights)) ps$weights else rep(1, length(y)),
    Gridtype = 2,
    use_parallel = TRUE,
    use_opencl = TRUE,
    verbose = FALSE
  )
  summary(fit)
}
}
###############################################################################
## End of Cleveland OpenCL example
###############################################################################
