############################### Boston_centered OpenCL example (rlmb) ####################

data("Boston_centered")
head(Boston_centered)
summary(Boston_centered)

## Predictors are mean-centered (column means ~0)
predictors <- setdiff(names(Boston_centered), "medv")
colMeans(Boston_centered[predictors])

form <- medv ~
  crim + zn +
  indus + chas + nox + age + dis + rad + tax + ptratio + black + lstat + rm

## Independent Normal-Gamma (OpenCL path when available), via rlmb
\donttest{
if (glmbayesCore_has_opencl()) {
  ps <- Prior_Setup(form, gaussian(), data = Boston_centered)

  fit <- rlmb(
    n = 1000L,
    y = ps$y,
    x = as.matrix(ps$x),
    pfamily = dIndependent_Normal_Gamma(
      ps$mu,
      ps$Sigma,
      shape = ps$shape_ING,
      rate = ps$rate
    ),
    use_parallel = TRUE,
    use_opencl = TRUE,
    verbose = FALSE
  )
  summary(fit)
}
}
###############################################################################
## End of Boston_centered OpenCL example
###############################################################################
