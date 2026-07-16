# glmbayesCore 0.1.0

## Initial CRAN submission

* First CRAN release of **glmbayesCore**. This package provides the compiled
  sampling engine and related “core” R interfaces (prior families, matrix-level
  samplers such as `rglmb` / `rlmb`, envelope construction, and optional OpenCL
  acceleration) that previously lived inside **glmbayes** and are now packaged
  separately for reuse.

* End users of Bayesian GLM/LM modelling should continue to use **glmbayes**
  for the formula interface and S3 methods. **glmbayesCore** is the developer
  / backend layer those packages build on.

* Future work includes pointing **glmbayes** (and other downstream packages)
  at this Core package as their shared sampling engine, so that the iid
  envelope stack is maintained in one place.
