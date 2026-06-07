## Conditionally independent block GLM draw (Poisson, 3 outcome groups)

set.seed(42)

## Dobson (1990) p. 93 RCT
counts <- c(18, 17, 15, 20, 10, 20, 25, 13, 12)
outcome <- gl(3, 1, 9)
treatment <- gl(3, 3)
d.AD <- data.frame(outcome, treatment, counts)

ps <- glmbayesCore::Prior_Setup(counts ~ treatment, family = poisson(), data = d.AD)
y <- ps$y
x <- ps$x
block <- outcome

out <- block_rNormalGLM(
  n = 1L,
  y = y,
  x = x,
  block = block,
  prior_list = list(mu = ps$mu, Sigma = ps$Sigma),
  family = poisson(),
  use_parallel = FALSE
)

out$coefficients
out$coef.mode
