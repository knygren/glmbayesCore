## Prior_SetupGroup — independent Prior_Setup() per Species block (iris)

data("iris", package = "datasets")

ps_block <- Prior_SetupGroup(
  Sepal.Length ~ Sepal.Width + Petal.Length,
  group = "Species",
  data = iris,
  family = gaussian()
)

print(ps_block)
names(ps_block)
ps_block$setosa$mu

## Per-block pfamilies via pfamily_list(ps_block)
pf <- pfamily_list(ps_block)
names(pf)
pf$setosa$prior_list$mu
