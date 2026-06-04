# Legacy / optional code (not part of the installed package)

Files here are excluded from `R CMD build` via `.Rbuildignore`.
Copy or source manually when experimenting; they are not exported from **glmbayes**.

## `pp_check.glmb.R`

Former S3 method `pp_check.glmb()` (wrapper around **bayesplot** for `glmb` fits).
Removed from `R/` so **glmbayes** does not depend on **bayesplot** (~41 packages on R-universe). **bayesplot** is not in **Imports** or **Suggests**.

To use locally:

1. `install.packages("bayesplot")`
2. `source("legacy_code/pp_check.glmb.R")` from the package root (after `library(glmbayes)`)
3. Call `bayesplot::pp_check(fit, ...)` or register the method by sourcing the file

Vignette chunks that previously called this method are commented out; see
`vignettes/Chapter-05.Rmd`, `Chapter-06.Rmd`, and `Chapter-12.Rmd` (in **glmbayes**).
