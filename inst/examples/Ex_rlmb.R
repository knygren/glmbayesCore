## Main Example based on Dobson Plant Weight Data 
## Use demo(Ex_07_Schools) for a longer/more complex model

## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
## Page 9: Plant Weight Data.

ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
group  <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
weight <- c(ctl, trt)

ps <- Prior_Setup(weight ~ group)
p  <- ncol(ps$x)
x  <- ps$x
mu <- ps$mu
V  <- ps$Sigma
y <- ps$y
shape    <- ps$shape
rate     <- ps$rate


## Two-Block Gibbs sampler for Plant Weight regression model
set.seed(180)

n_burnin=1000
n_samples=1000

## Initilize dispersion to ML estimate
dispersion2 <- ps$dispersion

## Run  burn-in iterations
for (i in 1:n_burnin) {
  ## Update block for regression coefficients
  out1 <- rlmb( n = 1, y = y, x = x,
    pfamily = dNormal(mu = mu, Sigma = V, dispersion = dispersion2) )
  
  ## Update block for dispersion
  out2 <- rlmb(n = 1, y = y, x = x,
  pfamily = dGamma(shape = shape, rate = rate, beta = out1$coefficients[1, ]))
  dispersion2 <- out2$dispersion
}

## Create Objects to store outputs
beta_out <- matrix(0, nrow = n_samples, ncol = 2)
disp_out <- rep(0, n_samples)

for (i in 1:n_samples) {
  ## Update block for regression coefficients
  out1 <- rlmb( n = 1, y = y, x = x,
                pfamily = dNormal(mu = mu, Sigma = V, dispersion = dispersion2) )

  ## Update block for dispersion
  out2 <- rlmb(n = 1, y = y, x = x,
               pfamily = dGamma(shape = shape, rate = rate, beta = out1$coefficients[1, ]))
  dispersion2 <- out2$dispersion
  
  ## Store output
  
  beta_out[i, 1:2] <- out1$coefficients[1, 1:2]
  disp_out[i]      <- out2$dispersion
}

mcmc_two_block <- coda::mcmc(cbind(    beta1 = beta_out[, 1],beta2 = beta_out[, 2],
                                       dispersion = disp_out  ))
  
## Review output
cat("\nCODA summary (Two-block Gibbs):\n")
print(summary(mcmc_two_block))
cat("\nEffective sample size (dispersion):\n")
print(coda::effectiveSize(mcmc_two_block)["dispersion"])


## Same model using the lmb function

lmb.D9 <- lmb(n= 4000,weight  ~ group,
  pfamily = dIndependent_Normal_Gamma(ps$mu, ps$Sigma, shape = ps$shape + p / 2, rate  = ps$rate))

## lmb summary
summary(lmb.D9)


## rlmb with dGamma prior (dispersion-only; coefficients fixed)
out_rlmb_dGamma <- rlmb(n = 100, y = y, x = x,
  pfamily = dGamma(shape = shape, rate = rate, beta = ps$coefficients),
  weights = rep(1, length(y)))
summary(out_rlmb_dGamma)

