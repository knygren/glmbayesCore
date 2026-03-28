
## ----dobson-------------------------------------------------------------------
## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
outcome <- gl(3, 1, 9)
treatment <- gl(3, 3)
print(d.AD <- data.frame(treatment, outcome, counts))

## ----Prior_and_Calls,results = "hide"-----------------------------------------
ps <- Prior_Setup(counts ~ outcome + treatment, family = poisson())
## Call to glmb
glmb.D93 <- glmb(
  n = 1000,
  counts ~ outcome + treatment,
  family = poisson(),
  pfamily = dNormal(mu = ps$mu, Sigma = ps$Sigma)
)

## ----glmb nobs----------------------------------------------------------------
## nobs for the glmb function
nobs(glmb.D93)


