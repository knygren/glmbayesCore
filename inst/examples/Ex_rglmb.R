data(menarche2)

summary(menarche2)
plot(Menarche/Total ~ Age, data=menarche2)
design_df <- data.frame(
  Age = menarche2$Age,
  Age2 = menarche2$Age - 13,
  Proportion = menarche2$Menarche / menarche2$Total,
  Total = menarche2$Total
)
x <- model.matrix(~ Age2, data = design_df)
y <- design_df$Proportion
wt <- design_df$Total

# Extract coefficient names from design matrix
coef_names <- colnames(x)

# Set up prior mean with names
mu <- matrix(0, nrow = length(coef_names), ncol = 1)
mu[2, 1] <- (log(0.9 / 0.1) - log(0.5 / 0.5)) / 3
rownames(mu) <- coef_names

# Set up prior covariance matrix with named rows and columns
V1<-1*diag(as.numeric(2.0))

# 2 standard deviations for prior estimate at age 13 between 0.1 and 0.9
## Specifies uncertainty around the point estimates
V1[1, 1] <- ((log(0.9 / 0.1) - log(0.5 / 0.5)) / 2)^2
V1[2, 2] <- (3 * mu[2, 1] / 2)^2 # Allows slope to be up to 3 times as large as point estimate 
rownames(V1) <- coef_names
colnames(V1) <- coef_names

out<-rglmb(n = 1000, y=y, x=x, pfamily=dNormal(mu=mu,Sigma=V1), weights = wt, 
           family = binomial(logit)) 
summary(out)

# Add mean(out$iters to rglmb summary function)
mean(out$iters)
