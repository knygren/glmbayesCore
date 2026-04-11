## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
## Page 9: Plant Weight Data.
ctl <- c(4.17, 5.58, 5.18, 6.11, 4.50, 4.61, 5.17, 4.53, 5.33, 5.14)
trt <- c(4.81, 4.17, 4.41, 3.59, 5.87, 3.83, 6.03, 4.89, 4.32, 4.69)
group <- gl(2, 10, 20, labels = c("Ctl", "Trt"))
weight <- c(ctl, trt)

ps <- Prior_Setup(weight ~ group, gaussian())
p <- ncol(ps$x)
# mu <- ps$mu
# V  <- ps$Sigma
# mu[1, 1] <- mean(weight)

## Classical model
lm.D9 <- lm(weight ~ group, x = TRUE, y = TRUE)
summary(lm.D9)
vcov(lm.D9)
s <- summary(lm.D9)
disp_classical <- s$sigma^2
cat("Classical lm dispersion (sigma^2 = RSS/(n-p)):", disp_classical, "\n")

## Conjugate Normal Prior (fixed dispersion)
lmb.D9 <- lmb(n=1000,
  weight ~ group,
  pfamily = dNormal(mu = ps$mu, ps$Sigma, dispersion = ps$dispersion)
)
summary(lmb.D9)
vcov(lmb.D9)

## Conjugate Normal_Gamma Prior (second argument is Sigma_0 from Prior_Setup)
lmb.D9_v2 <- lmb(n=1000,
  weight ~ group,
  pfamily = dNormal_Gamma(
    ps$mu,
    Sigma_0 = ps$Sigma_0,
    shape = ps$shape,
    rate  = ps$rate
  )
)
summary(lmb.D9_v2)
vcov(lmb.D9_v2)

## Independent_Normal_Gamma_Prior (same mu, Sigma, rate as dNormal_Gamma; shape + p/2)

lmb.D9_v3 <- lmb(n=1000,
  weight ~ group,
  dIndependent_Normal_Gamma(
    ps$mu,
    ps$Sigma,
    shape = ps$shape + p / 2,
    rate  = ps$rate
  )
)
summary(lmb.D9_v3)

vcov(lm.D9)
0.99*vcov(lm.D9)
vcov(lmb.D9)
vcov(lmb.D9_v2)
vcov(lmb.D9_v3)

disp_classical
mean(lmb.D9$dispersion)
mean(lmb.D9_v2$dispersion)
mean(lmb.D9_v3$dispersion)


## anova 
anova(lmb.D9)

## lmb with dGamma prior (dispersion-only; coefficients fixed)
out_lmb_dGamma <- lmb(n=1000,
  weight ~ group,
  pfamily = dGamma(shape = ps$shape, rate = ps$rate, beta = ps$coefficients))
summary(out_lmb_dGamma)

