# glmbayes

![GitHub release (latest by date)](https://img.shields.io/github/v/release/knygren/glmbayes?label=version)
![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)
![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/knygren/glmbayes/R-CMD-check.yaml?label=R%20CMD%20Check)

glmbayes provides independent and identically distributed (iid) samples for Bayesian Generalized Linear Models (GLMs).
Its primary interface, glmb(), serves as a Bayesian analogue to R's glm() function, supporting Gaussian, Poisson,
Binomial, and Gamma families under log-concave likelihoods. Sampling for most models is performed using accept-reject
methods based on likelihood subgradients (Nygren and Nygren, 2006). For Gaussian models, the package also includes
lmb(), a Bayesian counterpart to R's lm().

The package includes a rich set of supporting tools for prior specification, model diagnostics, and method functions
that mirror those for lm() and glm(). Most functions are extensively documented, and a comprehensive set of vignettes
are available to guide users through the package's capabilities.

The package is currently available on R-Universe, with plans for a future CRAN submission. For recent updates and planned enhancements, see
https://github.com/knygren/glmbayes/blob/main/NEWS.md

## Installation

To install the current development version (excluding OpenCL functionality):

install.packages("glmbayes",
                 repos = c("https://knygren.r-universe.dev",
                           "https://cloud.r-project.org"))

To install a version suitable for large models with GPU acceleration, follow the instructions from 

**Chapter 12 — Large Models: GPU Acceleration using OpenCL**
https://knygren.r-universe.dev/articles/glmbayes/Chapter-12.html

## Minimal Working Example

    library(glmbayes)

    # Dobson (1990), p. 93: Randomized Controlled Trial
    counts <- c(18,17,15,20,10,20,25,13,12)
    outcome <- gl(3,1,9)
    treatment <- gl(3,3)
    print(d.AD <- data.frame(treatment, outcome, counts))

    ## Classical glm
    glm.D93 <- glm(counts ~ outcome + treatment,
                   family = poisson())

    ## Bayesian glmb
    # Step 1: Set up prior
    ps <- Prior_Setup(counts ~ outcome + treatment, family = poisson())
    mu <- ps$mu
    V  <- ps$Sigma

    # Step 2: Fit using glmb
    glmb.D93 <- glmb(counts ~ outcome + treatment,
                     family = poisson(),
                     pfamily = dNormal(mu = mu, Sigma = V))

    summary(glmb.D93)

## Supported families, links, and pfamilies

As with the glm() function, models are specified by providing a symbolic description of the linear predictor
(using a formula) and a description of the error distribution (using a family and a link function). In addition, glmb()
also requires a prior specification, provided through a pfamily object. The available combinations of likelihood families,
link functions, and compatible pfamilies are:

| family                     | Available Link Functions           | Compatible pfamilies                                      |
|----------------------------|------------------------------------|------------------------------------------------------------|
| Gaussian                   | identity                           | dNormal, dGamma, dNormal_Gamma, dIndependent_Normal_Gamma |
| Poisson / Quasi-Poisson    | log                                | dNormal                                                    |
| Binomial / Quasi-Binomial  | logit, probit, cloglog             | dNormal                                                    |
| Gamma                      | log                                | dNormal, dGamma                                            |

More specifically, each pfamily constructor requires parameters associated with the prior. The available pfamilies and their usage signatures are:

- dNormal(mu, Sigma, dispersion = NULL)
- dGamma(shape, rate, beta, disp_lower = NULL, disp_upper = NULL)
- dNormal_Gamma(mu, Sigma, shape, rate)
- dIndependent_Normal_Gamma(
    mu,
    Sigma,
    shape,
    rate,
    max_disp_perc = 0.99,
    disp_lower = NULL,
    disp_upper = NULL
  )

To facilitate prior specification, the package provides a Prior_Setup() function, which extracts the needed prior parameters based on the same
symbolic model description and family specification used by glm(). By default, Prior_Setup() returns a reasonable prior specification
(described elsewhere), and optional arguments allow users to request alternative prior structures.

All supported models feature log-concave likelihoods, enabling efficient iid sampling via enveloping functions and
subgradient-based accept-reject algorithms (for models where other standard iid sampling algorithms are unavailable).

## Examples and Demos

Use example() and demo() to explore built-in examples and demos for supported families and links:

    ## Bayesian linear regression
    example("lmb")

    ## Bayesian Generalized Linear Models
    example("glmb")

    ## Two-block Gibbs sampler compared to iid-sampling for linear model
    example("rlmb")

    ## Larger Binomial Model with GPU Acceleration (requires OpenCL)
    example("Cleveland")

    ## Rubin/Gelman's Hierarchical 8 schools Model
    demo("Ex_07_Schools")

## Methodology

For generalized linear models where well known sampling methods are unavailable, sampling follows the
framework from Nygren and Nygren (2006), using likelihood subgradients to construct enveloping functions for
the posterior distribution. When the posterior is approximately normal, the expected number of draws per
acceptance is bounded as per that paper and as discussed in our vignettes.

## GPU Acceleration Using OpenCL

The implemented algorithms tend to have acceptable performance on CPUs up to around 10-14 dimensions.
For larger models, the envelope construction is embarrassingly parallel. To accelerate envelope construction
in such cases, the package provides optional GPU acceleration using OpenCL. This requires that users have
GPU enabled machines and an OpenCL installation. These features are discussed in more detail in two of
our vignettes.

## Vignettes

The glmbayes package includes a comprehensive set of vignettes organized into five major parts.
These vignettes guide users from introductory material through applied modeling, advanced topics,
and the underlying simulation methods that support the package.

### Part 1: An Introduction
Overview of the package, its design philosophy, and the basic workflow for
fitting Bayesian linear and generalized linear models. It introduces the core functions, model
objects, and the structure of the modeling interface.

- **Chapter 00 — Introduction**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-00.html

- **Chapter 01 — Getting Started with glmbayes**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-01.html

### Part 2: Estimating Bayesian Linear Models
These chapters focus on Bayesian linear regression using the Gaussian family. Topics include
model fitting, prior construction, posterior summaries, predictions, and deviance residuals.
This part establishes the foundation for understanding the Bayesian GLM framework used throughout
the package.

- **Chapter 02 — Estimating Bayesian Linear Models**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-02.html

- **Chapter 03 — Tailoring Priors - Leveraging the Prior_Setup Function**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-03.html

- **Chapter 04 — Reviewing Model Predictions, Deviance Residuals and Model Statistics**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-04.html

### Part 3: Generalized Linear Models
This part presents Bayesian GLMs across the major likelihood families, including binomial,
quasi-binomial, Poisson, quasi-Poisson, and Gamma models. It covers model specification,
link functions, log-concavity, diagnostics, and interpretation of posterior results.

- **Chapter 05 — Foundations of GLMs - Families, Links, and Log-Concave Likelihoods**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-05.html

- **Chapter 06 — Estimating Bayesian Generalized Linear Models**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-06.html

- **Chapter 07 — Models for the Binomial Family**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-07.html

- **Chapter 08 — Models for the Poisson Family**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-08.html

- **Chapter 09 — Models for the Gamma Family**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-09.html

### Part 4: Advanced Topics
These chapters explore more complex modeling scenarios and computational strategies, such as
informative priors, two-block Gibbs sampling, models with
unknown dispersion parameters, and large-scale model fitting using GPU acceleration 
using OpenCL.

- **Chapter 10 — Informative Priors: Centering and priors with differentiate weights**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-10.html

- **Chapter 11 — Estimating Models with unknown dispersion parameters**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-11.html

- **Chapter 12 — Large Models: GPU Acceleration using OpenCL**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-12.html

- **Chapter 13 — Block Gibbs Sampling Using rlmb and rglmb**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-13.html

### Part 5: Simulation Methods and Technical Implementation
This part documents the mathematical and algorithmic foundations of the package. Topics include
estimation procedures, likelihood subgradient densities, envelope construction, accept-reject
sampling, and technical reports on sampler design including implementation aspects for GPU acceleration using
OpenCL.

- **Chapter A01 — A detailed overview of the glmbayes package**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-A01.html

- **Chapter A02 — Overview of Estimation Procedures**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-A02.html

- **Chapter A03 — Methods Available in glmbayes**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-A03.html

- **Chapter A04 — Directional Tail Diagnostics for Prior-Posterior Disagreement**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-A04.html

- **Chapter A05 — Simulation Methods - Likelihood Subgradient Densities**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-A05.html

- **Chapter A06 — Accept-Reject Sampling for Dispersion in Gamma Regression**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-A06.html

- **Chapter A07 — Accept-Reject Sampling for gaussian Regression models with independent normal-gamma priors**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-A07.html

- **Chapter A08 — Overview of Envelope Related Functions**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-A08.html

- **Chapter A09 — Parallel Sampling Implementation using RcppParallel**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-A09.html

- **Chapter A10 — Accelerated EnvelopeBuild Implementation using OpenCL**  
https://knygren.r-universe.dev/articles/glmbayes/Chapter-A10.html


Together, these vignettes form a comprehensive reference that supports users at all levels, from first-time Bayesian GLM users to researchers interested in the mathematical and computational details behind the samplers.

## Feature Highlights

- S3 interface mirroring the structure of base glm()
- Accept-reject sampling for log-concave likelihoods
- Samplers for both fixed and variable dispersion
- Extensive vignettes to guide users through the package's capabilities
- Modular prior setup and checking tools

## Limitations

- Non-log-concave likelihoods are not currently supported
- Dispersion estimation requires a second sampler (rglmbdisp())

## Future Plans

- Full CRAN submission 
- Make Gamma regression for dGamma prior consistent with gamma.dispersion() function
- Update proofs
   - Vignette A6 to reflect similar approach as gamma.dispersion
   - Vignette A7 to reflect latest update to algorithm (with optimized bounds)
- Draft and Modify additional Vignettes
   - Chapter-10: Informative Priors: Centering and Priors with differentiated prior weights
   - Appendix Section: Restructure and Improve
- Improve References for Vignettes and Select Functions