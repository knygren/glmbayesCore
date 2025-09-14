## Dobson (1990) Page 93: Randomized Controlled Trial :
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))
## Set up Prior for Poisson Model
ps=Prior_Setup(counts ~ outcome + treatment, family=poisson())
ps

## Normal prior for glmb
glmb.D93<-glmb(counts ~ outcome + treatment, family = poisson(),pfamily=dNormal(mu=ps$mu,Sigma=ps$Sigma))


## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
## Page 9: Plant Weight Data.
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)

## Set up prior for gaussian model
ps2=Prior_Setup(weight ~ group,family=gaussian())
ps2

## Conjugate Normal Prior (fixed dispersion)
lmb.D9=lmb(weight ~ group,pfamily=dNormal(mu=ps2$mu,ps2$Sigma,dispersion = ps2$dispersion))

## Conjugate Normal_Gamma Prior 
lmb.D9_v2=lmb(weight ~ group,pfamily=dNormal_Gamma(ps2$mu,ps2$Sigma/ps2$dispersion,shape=ps2$shape,rate=ps2$rate))

## Independent_Normal_Gamma_Prior 
lmb.D9_v3=lmb(weight ~ group,dIndependent_Normal_Gamma(ps2$mu,ps2$Sigma,shape=ps2$shape,rate=ps2$rate))
