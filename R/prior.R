#' Setup Prior Objects
#' 
#' Helper function to facilitate the Setup of Prior Distributions for glm models.
#' @name Prior_Setup
#' @param na.action how \code{NAs} are treated. The default is first, any \code{\link{na.action}} attribute of 
#' data, second a \code{na.action} setting of \link{options}, and third \code{na.fail} if that is unset. 
#' The \code{factory-fresh} default is \code{na.omit}. Another possible value is \code{NULL}.
#' @param family a description of the error distribution and link function to be used in the model.
#' @param pwt Weight on the prior relative to the likelihood function at the maximum likelihood 
#' estimate. If supplied, this value is used directly (scalar or one value per coefficient). 
#' If \code{n_prior} is provided and \code{pwt} is still a **scalar** and \code{sd} was **not**
#' supplied, \code{pwt} is set to \code{n_prior / (n_prior + n_effective)}. If \code{length(pwt) > 1}
#' (including from \code{sd}) or \code{sd} was supplied, \code{n_prior} does **not** overwrite
#' \code{pwt}; it is used only as a scalar for Gamma / S_marg steps. If \code{sd} is provided,
#' \code{pwt} is computed from the prior standard deviations. If none of these are supplied,
#' \code{pwt} defaults to \code{pwt_default_low} for models with fewer than 14 coefficients, and 
#' \code{pwt_default_high} otherwise.
#' @param pwt_default_low Default prior weight used when \code{pwt} is not supplied and the model 
#' dimension is below 14. Defaults to 0.01.
#' @param pwt_default_high Default prior weight used when \code{pwt} is not supplied and the model 
#' dimension is 14 or greater. Defaults to 0.05.
#' @param n_prior Optional scalar effective prior sample size (on the \code{n_effective} scale).
#' If provided with scalar \code{pwt} and without \code{sd}, \code{pwt} is recomputed from
#' \code{n_prior}. With vector \code{pwt} or with \code{sd}, \code{pwt} is left unchanged and
#' \code{n_prior} is used for the Gamma prior on precision and related Gaussian calibration only.
#' If missing and \code{pwt} is scalar, \code{n_prior = (pwt/(1-pwt))*n_effective}.
#' @param sd Optional vector argument with the prior standard deviations for the coefficients
#' @param shape_df How the Gamma **shape** on residual precision is built from the scalar
#'   \eqn{n_{\mathrm{prior}}} and the number of coefficients \eqn{p=\texttt{ncol}(x)}.
#'   Here **df** means the **numerator** (effective prior count) before dividing by 2 for the
#'   Gamma shape---not residual degrees of freedom and not a \code{data.frame}.
#'   \describe{
#'     \item{\code{"n_prior"}}{(default) \eqn{\texttt{shape} = n_{\mathrm{prior}}/2}. With
#'       \code{\link{dIndependent_Normal_Gamma}}, the weak-prior limit of the log-target is
#'       close to the **joint** Gaussian log-likelihood in \eqn{(\beta,\tau)} (full-sample
#'       counting for dispersion). With \code{\link{dNormal_Gamma}}, the conjugate conditional
#'       \eqn{\pi(\beta\mid\tau)} adds an extra \eqn{(p/2)\log\tau} term, so the **same**
#'       \code{shape} implies **different** limiting dispersion behavior across these
#'       \code{pfamily}s.}
#'     \item{\code{"n_prior+p"}}{\eqn{\texttt{shape} = (n_{\mathrm{prior}}+p)/2}; \code{rate} is still
#'       \eqn{(n_{\mathrm{prior}}/2)\cdot\texttt{dispersion}} (same as \code{"n_prior"}), i.e. the usual
#'       \code{shape <- shape + p/2} with \code{rate} fixed. Intended for
#'       \code{\link{dIndependent_Normal_Gamma}} when you want marginal coefficient uncertainty
#'       in the **\eqn{n_{\mathrm{effective}}-p}** / \code{summary(\link[stats]{lm})} spirit
#'       as \code{pwt} approaches 0 (weak-prior limit). The weak-prior log-target **no longer** coincides with the raw
#'       log-likelihood because a \eqn{(p/2)\log\tau} prior contribution remains.}
#'     \item{\code{"n_prior-p"}}{\eqn{\texttt{shape} = (n_{\mathrm{prior}}-p)/2}; \code{rate} is still
#'       \eqn{(n_{\mathrm{prior}}/2)\cdot\texttt{dispersion}}. Intended for
#'       \code{\link{dNormal_Gamma}} when you want to **offset** the \eqn{(p/2)\log\tau} from
#'       \eqn{\pi(\beta\mid\tau)} so the total prior \eqn{\log\tau} term aligns with an
#'       independent Gamma\eqn{(n_{\mathrm{prior}}/2,\cdot)} specification---pushing weak-prior
#'       marginal dispersion toward **full-sample** (\eqn{n_{\mathrm{effective}}}) counting.
#'       **Requires** \eqn{n_{\mathrm{prior}} > p} so \eqn{\texttt{shape} > 0} (proper Gamma
#'       prior); the posterior may remain proper in some cases when this fails, but the prior
#'       is then invalid for the compound sampler.}
#'   }
#'   Ignored when \code{shape} and \code{rate} are not computed (non-scalar \code{n_prior}, or
#'   no dispersion for the family).
#' @param disp_type For Gaussian models, how \code{dispersion} and \code{rate} are set.
#'   \code{"Post_mean"} (default) chooses them so \code{dispersion} is **consistent with the
#'   posterior distribution** for dispersion under the prior and data (see Details).
#'   \code{"OLS_mean"} uses \eqn{\mathrm{RSS}_w/(n_{\mathrm{effective}}-2)}.
#' @param intercept_source Specifies the method through which the prior mean for the intercept term is set. Options are based on the null intercept only model (null_model) or full_models. The default is the null model which is safer if variables are not centered. 
#' @param effects_source Specifies the method through which the prior means for the effects terms are set. Options are null_effects (prior means set to zero) or full_model (effect means set to match maximum likelihood estimates).  
#' @param mu Optional vector argument with the prior means for the coefficients
#' @param x An object of class \code{"PriorSetup"}
#' @inheritParams stats::glm
#' @inheritParams stats::model.frame
#' @details
#' `Prior_Setup()` initializes a structured set of prior parameters for generalized linear models (GLMs), supporting both Gaussian and non-Gaussian families.
#'  It is designed to provide the full set of inputs required for multiple prior specifications (referred to as "pfamilies"), including conjugate normal priors, 
#'  Normal-Gamma priors, and independent normal-gamma priors.
#'
#' The function returns a list containing:
#' * `mu`: the prior mean vector
#' * `Sigma`: the prior variance-covariance matrix
#' * `dispersion`: the estimated dispersion (for Gaussian models)
#' * `shape`: the shape parameter for Normal-Gamma priors (if applicable)
#' * `rate`: the rate parameter for Normal-Gamma priors (if applicable)
#' * `coefficients`: named numeric vector, same length and names as columns of `x`; a point
#'   estimate for the regression coefficients (currently `coef()` of the internal full-model
#'   GLM fit, i.e. the MLE under that model). Intended mainly as `beta` in \code{\link{dGamma}()}
#'   and related `prior_list` entries when updating dispersion with fixed coefficients.
#'   This is distinct from `mu`, which holds prior means.
#' * `model`: the model frame used to construct the design matrix
#' * `x`: the model matrix used
#' * `y`: the response used
#' * `call`: the matched call to `Prior_Setup()`
#' * `PriorSettings`: a list of metadata including:
#'   - `pwt`: prior weight (scalar or vector)
#'   - `n_prior`: effective prior sample size
#'   - `n_likelihood`: effective likelihood sample size
#'   - `intercept_source`: method used to set the prior mean for the intercept
#'   - `effects_source`: method used to set the prior mean for the effects
#'
#' ### Connection to \code{pfamily} constructors and \code{prior_list}
#'
#' Call \code{Prior_Setup(formula, family = \dots)} (with \code{data},
#' \code{weights}, etc. as needed) so \code{mu}, \code{Sigma}, and (for
#' Gaussian families) \code{dispersion}, \code{shape}, and \code{rate} are
#' on the scale implied by the likelihood and your \code{pwt} / \code{n_prior}.
#' The recommended mapping from the returned list into a \code{\link{pfamily}}
#' or into \code{prior_list} for \code{\link{simfuncs}} follows the worked
#' patterns in \code{example("Prior_Setup")} (\code{inst/examples/Ex_Prior_Setup.R}).
#' Distinct priors use **different** combinations of \code{Sigma} versus
#' \code{Sigma/dispersion}; mixing these up is a common error.
#'
#' #### \code{\link{dNormal}()}
#' * **Typical use (e.g. Poisson \code{\link{glmb}}):** \code{dNormal(mu = ps$mu, Sigma = ps$Sigma)}.
#' * **Gaussian with fixed dispersion (\code{\link{lmb}} / \code{\link{rlmb}}):**
#'   pass \code{dispersion = ps2$dispersion} as well:
#'   \code{dNormal(mu = ps2$mu, Sigma = ps2$Sigma, dispersion = ps2$dispersion)}.
#' * Matrix-input fits use the same \code{pfamily} with \code{y} and \code{x}
#'   taken from \code{ps$y} and \code{as.matrix(ps$x)} (see examples).
#'
#' #### \code{\link{dNormal_Gamma}()} (conjugate Normal--Gamma, Gaussian)
#' The second argument is the prior covariance on the **precision-weighted**
#' coefficient scale: use **\code{ps2$Sigma / ps2$dispersion}**, not \code{ps2$Sigma}
#' alone. Recommended call:
#' \code{dNormal_Gamma(ps2$mu, ps2$Sigma / ps2$dispersion, shape = ps2$shape, rate = ps2$rate)}.
#' For \code{\link{rNormalGamma_reg}}, build
#' \code{prior_list = list(mu = ps2$mu, Sigma = ps2$Sigma / ps2$dispersion, shape = ps2$shape, rate = ps2$rate)}
#' (again **divided** \code{Sigma}).
#'
#' #### \code{\link{dIndependent_Normal_Gamma}()} (Gaussian, non-conjugate joint \eqn{(\beta,\phi)})
#' Here the second argument is the **full** prior covariance on \eqn{\beta}:
#' \code{dIndependent_Normal_Gamma(ps2$mu, ps2$Sigma, shape = ps2$shape, rate = ps2$rate)}.
#' Same pattern for \code{\link{rglmb}} / \code{\link{rlmb}} with that \code{pfamily}.
#' For \code{\link{rindepNormalGamma_reg}}, the template \code{prior_list} also
#' uses **undivided** \code{Sigma}, plus \code{dispersion}, \code{Precision = solve(Sigma)},
#' and envelope controls such as \code{max_disp_perc} (see examples).
#'
#' #### \code{\link{dGamma}()} (Gamma on precision / dispersion with fixed \eqn{\beta})
#' \code{Prior_Setup()} supplies \code{shape}, \code{rate}, and \code{coefficients}
#' (full-model MLE by default). Typical use:
#' \code{dGamma(shape = ps2$shape, rate = ps2$rate, beta = ps2$coefficients)}.
#' For \code{\link{rGamma_reg}}, pass
#' \code{prior_list = list(beta = ps2$coefficients, shape = ps2$shape, rate = ps2$rate)}.
#'
#' For a runnable comparison of \code{shape_df} settings with \code{\link{lmb}}, see the package
#' demo \code{demo(Ex_10_Prior_Setup_shape_df, package = "glmbayes")}.
#'
#' #### References and further reading
#' Zellner-style scaling of \code{Sigma} from the likelihood
#' \insertCite{zellner1986gprior}{glmbayes}; conjugate Gaussian theory
#' \insertCite{Raiffa1961}{glmbayes}; envelope sampling for \code{dNormal()}
#' on non-Gaussian families and for \code{dIndependent_Normal_Gamma()}
#' \insertCite{Nygren2006}{glmbayes}; Normal--Gamma / GLM background
#' \insertCite{Gelman2013,Dobson1990,McCullagh1989}{glmbayes}; prior tailoring
#' \insertCite{glmbayesChapter03}{glmbayes}.
#'
#' ### Inputs to the function
#'
#' The inputs to `Prior_Setup()` fall into three conceptual categories:
#'
#' **1. Model specification**
#' * `formula`: defines the structure of the GLM, including response and predictors.
#' * `family`: specifies the error distribution and link function (e.g., `gaussian()`, `binomial()`).
#' * `data`: optional data frame used to evaluate the formula and extract variables.
#'
#' **2. Prior variance-covariance specification**
#' * `pwt`: prior weight relative to the likelihood. If scalar, used to compute Zellner's g-prior.
#' * `n_prior`: optional scalar effective prior sample size. Replaces scalar `pwt` only when
#'   `pwt` is scalar and `sd` is not used; otherwise supplies precision-prior / calibration only.
#' * `sd`: optional vector of prior standard deviations. If provided, used to compute `pwt`.
#' * `shape_df`: maps scalar `n_prior` and `p = ncol(x)` to the Gamma **shape** on residual precision
#'   (see argument description and Details).
#'
#' **3. Prior mean specification**
#' * `intercept_source`: method for setting the prior mean of the intercept (`"null_model"` or `"full_model"`).
#' * `effects_source`: method for setting the prior mean of the effects (`"null_effects"` or `"full_model"`).
#' * `mu`: optional user-specified prior mean vector. Overrides other centering logic if provided.
#'
#' ### Mathematical structure and interpretation
#'
#' When `pwt` is a scalar, the prior covariance is scaled from the likelihood covariance using a Zellner g-prior:
#' \deqn{\Sigma_{Prior} = \frac{1 - pwt}{pwt} \cdot V_{MLE}}
#' where \eqn{V_{MLE}} is the variance-covariance matrix of the maximum likelihood estimator.
#'
#' If the likelihood covariance is not full rank, the function aborts with an error, as a g-prior cannot be constructed.
#'  While Bayesian models can still be estimated in such cases, users should proceed with caution.
#'  
#' The corresponding prior precision is:
#' \deqn{P_{Prior} = \frac{pwt}{1 - pwt} \cdot P_{MLE}}
#' and the posterior mean (in the Gaussian case) simplifies to:
#' \deqn{\mu_{Post} = (1 - pwt) \cdot \widehat{\beta}_{MLE} + pwt \cdot \mu_{Prior}}
#'
#' Default prior centering reflects classical reference structures:
#' \itemize{
#'   \item The intercept prior is centered at the null model estimate (`intercept_source = "null_model"`), consistent with how R\eqn{^2}, F-tests, and t-statistics are defined.
#'   \item The effect priors are centered at zero (`effects_source = "null_effects"`), supporting intuitive posterior summaries and direct interpretation of tail probabilities.
#' }
#' 
#' These defaults ensure that posterior means remain close to classical estimates when `pwt` is small, 
#' while allowing additional flexible prior specifications through `mu`, `sd`, and `n_prior`. 
#' 
#' If `n_prior` is not provided, it is derived from `pwt` and `n_effective` (stored as
#' \code{PriorSettings$n_effective}; for Gaussian with weights, \eqn{n_{\mathrm{effective}}=\sum w_i}):
#' \deqn{n_{\mathrm{prior}} = \frac{pwt}{1 - pwt} \cdot n_{\mathrm{effective}}.}
#' The field \code{n_likelihood} in \code{PriorSettings} is set equal to \code{n_effective} for compatibility.
#'
#' If `n_prior` is provided and `pwt` is still scalar and `sd` was not used, `pwt` is set to:
#' \deqn{pwt = \frac{n_{\mathrm{prior}}}{n_{\mathrm{prior}} + n_{\mathrm{effective}}}}
#' (With vector `pwt` or with `sd`, `pwt` is not overwritten.)
#'
#' When applicable, `Prior_Setup()` computes the shape and rate parameters for a Gamma prior on the residual precision (inverse variance), used in
#' compound prior families such as `dNormal_Gamma()`, `dIndependent_Normal_Gamma()`, and `dGamma()`. Let \eqn{p=\texttt{ncol}(x)}. The **shape** is
#' \deqn{\text{shape} = \frac{n_{\mathrm{shape}}}{2}}
#' where \eqn{n_{\mathrm{shape}}} is set by `shape_df`:
#' \describe{
#'   \item{\code{shape_df = "n_prior"}}{\eqn{n_{\mathrm{shape}} = n_{\mathrm{prior}}}}
#'   \item{\code{shape_df = "n_prior+p"}}{\eqn{n_{\mathrm{shape}} = n_{\mathrm{prior}} + p}}
#'   \item{\code{shape_df = "n_prior-p"}}{\eqn{n_{\mathrm{shape}} = n_{\mathrm{prior}} - p}
#'     (requires \eqn{n_{\mathrm{prior}} > p})}
#' }
#' The **rate** is tied to the **nominal** prior sample size \eqn{n_{\mathrm{prior}}}, not to the
#' adjusted **shape** when \code{shape_df} adds or subtracts \eqn{p}:
#' \deqn{\text{rate} = \frac{n_{\mathrm{prior}}}{2} \cdot \text{dispersion}}
#' So for \code{shape_df = "n_prior+p"}, \eqn{\texttt{shape} = (n_{\mathrm{prior}}+p)/2} but
#' \eqn{\texttt{rate}} is the same as for \code{shape_df = "n_prior"} with the same \eqn{n_{\mathrm{prior}}}
#' and \code{dispersion}. This matches the common recipe for \code{\link{dIndependent_Normal_Gamma}}:
#' increase Gamma **shape** by \eqn{p/2} while leaving **rate** at the \code{"n_prior"} value.
#' With \code{disp_type = "Post_mean"}, the joint search updates \code{dispersion} and then sets
#' \eqn{\texttt{rate} = (n_{\mathrm{prior}}/2)\cdot\texttt{dispersion}} at the optimum.
#'    
#' ### Gaussian dispersion and \code{disp_type}
#'
#' For \code{family = gaussian()}, the function always computes a **classical** weighted
#' residual-variance ratio from the internal \code{\link[stats]{glm.fit}} (not \code{\link[stats]{glm}}):
#' with \eqn{\mathrm{RSS}_w = \sum_i w_i r_i^2} and \eqn{n_{\mathrm{effective}} = \sum_i w_i},
#' \deqn{d_{\mathrm{OLS}} = \frac{\mathrm{RSS}_w}{n_{\mathrm{effective}} - 2}}
#' requiring \eqn{n_{\mathrm{effective}} > 2}. This is **not**
#' \code{summary(glm.fit(...))$dispersion}. That value is the **returned** \code{dispersion}
#' when \code{disp_type = "OLS_mean"}; when \code{disp_type = "Post_mean"}, it is used only as an
#' internal starting value.
#'
#' If \code{n_prior} is available so a Gamma prior on precision is defined (\code{shape},
#' \code{rate}), the returned \code{dispersion} and \code{rate} are set according to
#' \code{disp_type}:
#' * \code{disp_type = "Post_mean"} (default): choose \code{dispersion} and \code{rate} so the
#'   prior's dispersion input is **consistent with the posterior distribution** for residual
#'   dispersion (given the rest of the prior and the data).
#' * \code{disp_type = "OLS_mean"}: return \eqn{d_{\mathrm{OLS}}} as \code{dispersion} and
#'   \eqn{\mathrm{rate} = (n_{\mathrm{prior}}/2)\cdot \texttt{dispersion}} (for any \code{shape_df}).
#'
#' The posterior shape and rate for residual precision (Gaussian sampling fragment) are:
#' \deqn{\text{shape}_{\mathrm{post}} = \text{shape} + \frac{n_{\mathrm{effective}}}{2}
#'       = \frac{n_{\mathrm{shape}} + n_{\mathrm{effective}}}{2}}
#' \deqn{\text{rate}_{\mathrm{post}} = \text{rate} + \frac{\mathrm{RSS}_w}{2}.}
#' With the default \code{shape_df = "n_prior"}, \eqn{n_{\mathrm{shape}}=n_{\mathrm{prior}}} and
#' \eqn{\text{rate} = \texttt{dispersion}\,n_{\mathrm{prior}}/2}, so
#' \eqn{\text{rate}_{\mathrm{post}} = \texttt{dispersion}\,n_{\mathrm{prior}}/2 + \mathrm{RSS}_w/2}
#' as before.
#'
#' This structure allows the prior to contribute pseudo-observations to the residual precision estimate, enabling adaptive shrinkage and hierarchical 
#' regularization---especially valuable in small-sample or high-dimensional settings.
#'
#' @return A list with items related to the prior.
#' \item{mu}{A prior mean vector}
#' \item{Sigma}{A prior variance-covariance matrix}
#' \item{dispersion}{Empirical bayes estimate for the dispersion (gaussian model only)}
#' \item{shape}{Derived prior shape parameter (gaussian model only). Defaults to n_prior/2 where n_prior is derived from pwt if not provided}
#' \item{rate}{Derived prior rate parameter (gaussian model only). Defaults to (n_prior*dispersion)/2 where n_prior is derived from pwt if not provided}
#' \item{coefficients}{Named numeric vector of prior-implied posterior-mean coefficients.
#'   For \code{gaussian()} this uses the closed-form Zellner blend
#'   \eqn{(1-\texttt{pwt})\hat\beta+\texttt{pwt}\mu} (scalar or per-coefficient \code{pwt})
#'   when available; otherwise it falls back to the internal full-model GLM coefficients
#'   (same convention as \code{\link[stats]{coef}}).}
#' \item{model}{The model frame from \code{object} if it exists}
#' \item{x}{The design matrix from \code{object} if it exists}
#' \item{PriorSettings}{A list containing prior configuration details}
#' @family prior
#' @seealso
#' \code{\link{pfamily}} for prior-family objects and the constructors
#' \code{\link{dNormal}}, \code{\link{dNormal_Gamma}}, \code{\link{dGamma}},
#' and \code{\link{dIndependent_Normal_Gamma}}.
#'
#' \code{\link{glmb}}, \code{\link{lmb}} for formula-based fits with a
#' \code{pfamily} built from \code{Prior_Setup()} output; \code{\link{rglmb}},
#' \code{\link{rlmb}} for matrix-based sampling that consumes the same prior
#' structure; \code{\link{simfuncs}} for functions that take a \code{prior_list}
#' assembled from those components (including \code{\link{rindepNormalGamma_reg}}
#' for \code{\link{dIndependent_Normal_Gamma}()}).
#'
#' \insertCite{glmbayesChapter03}{glmbayes} for prior tailoring and examples.
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @example inst/examples/Ex_Prior_Setup.R
#' @export
#' @rdname Prior_Setup
#' @order 1

## Note arguments outside of first two are currently not used

Prior_Setup <- function(
    formula,
    family      = gaussian(),
    data=NULL,
    weights=NULL,
    subset=NULL,
    na.action   = na.fail,
    offset=NULL,
    contrasts   = NULL,
    pwt         = NULL,
    pwt_default_low = 0.01,      # new: low-d default
    pwt_default_high = 0.05,     # new: high-d default
    n_prior     = NULL,
    sd          = NULL,
    shape_df    = c("n_prior", "n_prior+p", "n_prior-p"),
    disp_type   = c("Post_mean", "OLS_mean"),
    intercept_source = c("null_model", "full_model"),
    effects_source   = c("null_effects",  "full_model"),
    mu          = NULL,  ...
  ) 
  
  {

  call <- match.call()  
  intercept_source <- match.arg(intercept_source)
  effects_source <- match.arg(effects_source)
  shape_df <- match.arg(shape_df)
  disp_type <- match.arg(disp_type)
  
  
  
  #mf<-model.frame(formula,data,subset=subset,na.action=na.action,
  #                drop.unused.levels=drop.unused.levels,xlev=xlev)
  
  
  if (is.character(family)) 
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) 
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  
  if (missing(data))   data <- environment(formula)
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
               "etastart", "mustart", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  
##  mf<-model.frame(formula,data)
  
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  
  if (family$family == "binomial"
      && is.numeric(Y) && is.vector(Y)
      && all(Y >= 0 & Y <= 1)    # use <= instead of <=
      && is.null(weights)) {
    warning(
      "You supplied a proportion response (0 <= y <= 1) to a binomial family\n",
      "without `weights`. Each case will be treated as a single trial (n=1).\n",
      "If you meant to model counts, either use `cbind(success, failure)`\n",
      "or supply `weights =` the number of trials."
    )
  }
  
  
##  X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts) else matrix(, NROW(Y), 0L)
  X <- if (!is.empty.model(mt)) model.matrix(mt, mf, ...) else matrix(, NROW(Y), 0L)
  
  ## --- WEIGHT HANDLING -------------------------------------------------------
  
  # Extract raw weights from the model frame (may be NULL)
  raw_wt <- model.weights(mf)
  
  # Number of observations (always correct)
  n_obs <- nrow(X)
  
  # Case 1: User supplied scalar weight (e.g., weights = 4)
  # model.frame() cannot accept scalar weights, so we expand them *after* mf is built
  if (!is.null(raw_wt) && length(raw_wt) == 1L) {
    weights <- rep(raw_wt, n_obs)
  }
  
  # Case 2: User supplied a full-length weight vector
  else if (!is.null(raw_wt)) {
    
    if (!is.numeric(raw_wt))
      stop("'weights' must be numeric")
    
    if (any(raw_wt < 0))
      stop("negative weights not allowed")
    
    if (length(raw_wt) != n_obs)
      stop("weights must be either a scalar or have length equal to number of observations")
    
    weights <- raw_wt
  }
  
  # Case 3: No weights supplied -> default depends on family
  else {
    
    if (family$family %in% c("gaussian", "Gamma")) {
      # Gaussian/Gamma treat weights as replication weights -> default = 1
      weights <- rep(1, n_obs)
    } else {
      # Binomial/Poisson: GLM semantics require weights = NULL
      weights <- NULL
    }
  }
  
  # Compute effective sample size
  if (is.null(weights)) {
    n_effective <- n_obs
  } else {
    n_effective <- sum(weights)
  }
  
  
  #######################################33
  
  
  
  
  
  
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(Y)) 
      stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                    length(offset), NROW(Y)), domain = NA)
  }
  
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")
  
  
  x<-X  
##  x<-model.matrix(formula,mf)
  
  nvar=ncol(x)
  
  if (is.null(pwt)) {
    pwt <- if (nvar < 14) pwt_default_low else pwt_default_high
    ## n_prior (later) recomputes pwt; avoid implying the default applies.
    if (is.null(n_prior)) {
      message("Using default pwt = ", pwt,
              " (", if (nvar < 14) "low-d" else "high-d", " default).")
    }
  }
  
  ## Make sure the *columns* of x are named correctly:
  
  
  ## validate pwt  
  if (!is.numeric(pwt) || any(is.na(pwt))) {  
    stop("pwt must be numeric and non NA, either length 1 or length ", nvar)  
  }  
  if (! (length(pwt) %in% c(1, nvar)) ) {  
    stop("pwt must have length 1 or length(coef) = ", nvar,  
         "; you supplied length ", length(pwt))  
  }  
  if (any(pwt <= 0 | pwt >= 1)) {  
    stop("All elements of pwt must lie strictly between 0 and 1; you supplied:",  
         paste0(round(pwt, 3), collapse = ", "))  
  }  
  
  
  
  var_names <- colnames(x)
  colnames(x) <- var_names

  mu_internal <- matrix(0, nrow = nvar, ncol = 1, dimnames = list(var_names, "mu"))
 
  
  

  glm_full <- glm.fit(
    x       = X,
    y       = Y,
    weights = weights,
    offset  = offset,
    family  = family
    ,control = glm.control(...)
  )
  

  glm_full$call      <- call
  glm_full$formula   <- formula
  glm_full$terms     <- mt
  glm_full$data      <- mf
  glm_full$offset    <- offset
  glm_full$contrasts <- attr(X, "contrasts")
  glm_full$xlevels   <- .getXlevels(mt, mf)
  class(glm_full)    <- c("glm", "lm")
  

  V0 <- vcov(glm_full)
  
  glm_summary=summary(glm_full)
  
  ##n_likelihood <- glm_summary$df.residual + glm_summary$df[1]  # residual df + model rank
  n_likelihood <- n_effective
  
  # If sd is provided, use it to compute pwt
if (!is.null(sd)) {
  if (!is.numeric(sd) || any(is.na(sd))) {
    stop("sd must be a numeric vector with no missing values.")
  }
  if (length(sd) != nvar) {
    stop("Length of sd must match number of coefficients (", nvar, ").")
  }

  # Compute pwt from sd and V0
  V0_diag <- diag(V0)
  if (any(V0_diag <= 0)) {
    stop("Diagonal entries of V0 must be positive to compute pwt from sd.")
  }

  pwt <- V0_diag / (V0_diag + sd^2)
  message("Computed pwt from user-specified prior standard deviations (sd).")
}
  
    
  ## n_prior may imply pwt only when pwt is still a single scalar not from `sd`.
  ## If length(pwt) > 1 (vector pwt) or `sd` was supplied, do not overwrite pwt;
  ## n_prior is then used only downstream as a scalar for Gamma / S_marg remap.
  if (!is.null(n_prior)) {
    if (!is.numeric(n_prior) || length(n_prior) != 1 || n_prior <= 0) {
      stop("n_prior must be a single positive numeric value")
    }
    if (length(pwt) == 1L && is.null(sd)) {
##    pwt <- n_prior / (n_prior + n_likelihood)
      pwt <- n_prior / (n_prior + n_effective)
      message("Computed pwt = ", round(pwt, 4),
              " from n_prior = ", n_prior,
              " and n_effective = ", n_effective)
    }
  }
  
  # Compute n_prior if not supplied and pwt is scalar
  if (is.null(n_prior) && length(pwt) == 1L) {
    n_prior <- (pwt/(1-pwt)) * n_effective
  #  message("Computed n_prior = ", round(n_prior, 4),
  #          " from pwt = ", round(pwt, 4),
  #          " and n_likelihood = ", n_likelihood)
  }
  if (identical(family$family, "gaussian") && is.null(n_prior)) {
    stop(
      "For Gaussian models, a scalar effective prior sample size `n_prior` is required. ",
      "Use scalar `pwt` (it implies `n_prior`), or supply `n_prior` explicitly. ",
      "Per-coefficient `sd` implies vector `pwt`; in that case you must pass `n_prior`.",
      call. = FALSE
    )
  }

  ## --- CONDITIONAL DISPERSION (Gaussian): explicit ratio from glm.fit object -----
  ## Uses stats::glm.fit (not glm()). Default: dispersion = RSS_w / (n_effective - 2).
  ## # Old MLE-style ratio (retained for reference, not used):
  ## # dispersion <- rss_weighted / n_effective
  ## With rate = dispersion * shape (shape = n_prior/2), posterior summaries of tau = 1/d
  ## depend on pwt unless additional structure holds (see Details).
  rss_weighted_stored <- NA_real_
  dispersion_classical <- NA_real_
  if (family$family == "gaussian") {
    res <- residuals(glm_full, type = "response")
    w   <- glm_full$prior.weights
    rss_weighted <- sum(w * res^2)
    if (!is.finite(rss_weighted) || rss_weighted <= 0) {
      stop("Weighted RSS must be strictly positive for Gaussian dispersion priors.")
    }
    if (!is.finite(n_effective) || n_effective <= 0) {
      stop("n_effective must be strictly positive to compute Gaussian dispersion.")
    }
    if (n_effective <= 2) {
      stop(
        "Gaussian dispersion requires n_effective > 2 ",
        "(denominator n_effective - 2)."
      )
    }
    dispersion <- rss_weighted / (n_effective - 2)
    if (!is.finite(dispersion) || dispersion <= 0) {
      stop("Computed Gaussian dispersion must be strictly positive.")
    }
    rss_weighted_stored <- rss_weighted
    dispersion_classical <- dispersion
    
  } else if (family$family == "Gamma") {
    
    # MASS::gamma.dispersion() already returns the correct quasi-likelihood
    # dispersion estimate for Gamma GLMs.
    dispersion <- MASS::gamma.dispersion(glm_full)
    
  } else {
    
    dispersion <- NULL
  }
  

    if (!is.matrix(V0) || nrow(V0) != ncol(V0)) {
    stop("vcov(glm_full) (V0) must be a square matrix.")
  }
  if (anyNA(V0)) {
    stop("vcov(glm_full) (V0) contains missing values.")
  }
  
  # 2. symmetry (up to numerical tolerance)
  if (!isSymmetric(V0, tol = sqrt(.Machine$double.eps))) {
    stop("vcov(glm_full) (V0) is not symmetric.")
  }
  
  # 3. positive-definiteness via Cholesky
  pd_try <- try(chol(V0), silent = TRUE)
  if (inherits(pd_try, "try-error")) {
    stop(
      "Variance-covariance matrix V0 is not positive-definite.\n",
      "This usually means the classical GLM is rank-deficient."
    )
  }

  if (var_names[1] == "(Intercept)") {
    # build 1-column design matrix for intercept only
    X0 <- matrix(1, nrow = NROW(Y), ncol = 1,
                 dimnames = list(NULL, "(Intercept)"))
    
    # fit intercept-only model via glm.fit()
    fit0 <- glm.fit(
      x       = X0,
      y       = Y,
      weights = weights,
      offset  = offset,
      family  = family,
      control = glm.control(...)
    )
    
    # pick the intercept from null or full model
    chosen_int <- switch(
      intercept_source,
      null_model = fit0$coefficients[1],
      full_model = glm_full$coefficients[1]
    )
    
    mu_internal[1, 1] <- chosen_int
  }
  

  # 5) effects prior means
  if (nvar > 1) {
    effect_names <- var_names[-1]
    if (effects_source == "full_model") {
      coefs <- coef(glm_full)[effect_names]
      mu_internal[effect_names, 1] <- coefs
    }
    # else null_effects leaves mu[...] as zero
  }
  

  # Validate user-supplied mu if provided
  if (!is.null(mu)) {
    if (!is.numeric(mu)) {
      stop("mu must be numeric.")
    }
    if (is.vector(mu)) {
      if (length(mu) != nvar) {
        stop("Length of mu vector must match number of coefficients (", nvar, ").")
      }
      mu <- matrix(mu, ncol = 1, dimnames = list(var_names, "mu"))
    } else if (is.matrix(mu)) {
      if (!all(dim(mu) == c(nvar, 1))) {
        stop("mu matrix must have dimensions [", nvar, ", 1].")
      }
      rownames(mu) <- var_names
      colnames(mu) <- "mu"
    } else {
      stop("mu must be either a numeric vector or a matrix.")
    }
    message("Using user-specified prior mean vector (mu).")
  } else {
    mu <- mu_internal
  }
  
    
  
#  Sigma=as.matrix(diag(nvar))
 
#  Sigma=(1-pwt)/pwt*V0
   
  ## build prior covariance  
  if (length(pwt) == 1L) {  
    ## full matrix prior  
    Sigma <- ((1 - pwt) / pwt) * V0  
  }
  else {  
    scale_vec <- sqrt((1 - pwt) / pwt)
    scale_mat <- outer(scale_vec, scale_vec)
    Sigma <- V0 * scale_mat
  }  

  rownames(mu)=var_names
  colnames(mu)=c("mu")
  rownames(Sigma)=var_names
  colnames(Sigma)=var_names

  ## --- d_P (diagnostic only; not returned or used downstream) -----------------
  ## S = y'Wy + mu' P0 mu - mu_n' Pn mu_n with P0 = Sigma^{-1}, PL = X' W X,
  ## Pn = P0 + PL, mu_n = Pn^{-1}(P0 mu + X'Wy).  d_P = S / (n_effective - 2).
  d_P <- NA_real_
  if (identical(family$family, "gaussian")) {
    if (is.finite(n_effective) && n_effective > 2) {
      y_num <- as.numeric(Y)
      if (length(y_num) == n_obs) {
        w_vec <- if (is.null(weights)) rep(1, n_obs) else as.numeric(weights)
        P0 <- tryCatch(solve(Sigma), error = function(e) NULL)
        if (!is.null(P0)) {
          XtW <- sweep(X, 1, w_vec, `*`)
          PL <- crossprod(XtW, X)
          Pn <- P0 + PL
          mu_mat <- matrix(as.numeric(mu), ncol = 1L)
          rhs <- P0 %*% mu_mat + crossprod(X, w_vec * y_num)
          mu_n <- solve(Pn, rhs)
          ytWy <- sum(w_vec * y_num^2)
          S_quad <- ytWy +
            as.numeric(crossprod(mu_mat, P0 %*% mu_mat)) -
            as.numeric(crossprod(mu_n, Pn %*% mu_n))
          d_P <- S_quad / (n_effective - 2)
        }
      }
    }
  }

  ## Gamma on precision: shape = n_shape_num/2, rate = dispersion * (n_prior/2).
  ## n_shape_num from shape_df: n_prior, n_prior+p, or n_prior-p (latter needs n_prior > p).
  ## Rate uses n_prior/2 always so n_prior+p matches "shape += p/2, rate unchanged" vs n_prior.
  ## In `disp_type="Post_mean"`, `Sigma` is rescaled below; `dispersion` stays d_OLS and `rate` stays (n_prior/2)*dispersion.
  dispersion_for_shape_rate <- dispersion
  if (!is.null(n_prior) && length(n_prior) == 1L && !is.null(dispersion_for_shape_rate)) {
    ## n_prior is interpreted as effective prior sample size, on the same scale as sum(weights).
    p_coef <- nvar
    n_shape_num <- switch(
      shape_df,
      "n_prior"   = n_prior,
      "n_prior+p" = n_prior + p_coef,
      "n_prior-p" = {
        if (!is.finite(n_prior) || !is.finite(p_coef) || n_prior <= p_coef) {
          stop(
            "shape_df = \"n_prior-p\" requires n_prior > p (number of coefficients). ",
            "Got n_prior = ", n_prior, " and p = ", p_coef, "."
          )
        }
        n_prior - p_coef
      }
    )
    shape <- n_shape_num / 2
    if (!is.finite(shape) || shape <= 0) {
      stop("Computed shape must be strictly positive.")
    }
    rate <- dispersion_for_shape_rate * (n_prior / 2)
    if (!is.finite(rate) || rate <= 0) {
      stop("Computed rate must be strictly positive.")
    }
  } else {
    shape <- NULL
    rate <- NULL
  }

  ## --- S_marg_new / S_marg_sigma0_vcov before Post_mean / Nelder (full Sigma_pre_nm: scalar or vector pwt)
  ## Sigma_0 = Sigma_pre_nm / d with d = d_OLS or d_vcov = summary(glm)$dispersion (cancels vcov scale in V0).
  ## S_marg keeps the same value for downstream b_0 / remap logic (identical to S_marg_new here).
  S_marg <- NA_real_
  S_marg_new <- NA_real_
  S_marg_sigma0_vcov <- NA_real_
  S_marg_post_mean <- NA_real_
  S_marg_scalar_zellner <- NA_real_
  Sigma_pre_nm <- Sigma
  if (identical(family$family, "gaussian") &&
      is.finite(rss_weighted_stored) && rss_weighted_stored > 0 &&
      is.finite(dispersion_classical) && dispersion_classical > 0) {
    bh_sm <- coef(glm_full)
    if (length(bh_sm) == nvar && all(is.finite(bh_sm))) {
      w_sm <- if (is.null(weights)) rep(1, n_obs) else as.numeric(weights)
      XtW_sm <- sweep(X, 1, w_sm, `*`)
      Gm_sm <- crossprod(XtW_sm, X)
      dlt_sm <- matrix(bh_sm, ncol = 1L) - matrix(as.numeric(mu), ncol = 1L)
      Sigma_0_pre_nm <- Sigma_pre_nm / dispersion_classical
      Ginv_sm <- tryCatch(
        solve(Gm_sm),
        error = function(e) {
          stop(
            "Prior_Setup: cannot invert weighted Gram matrix X'WX (S_marg_new). ",
            conditionMessage(e),
            call. = FALSE
          )
        }
      )
      M_sm <- Sigma_0_pre_nm + Ginv_sm
      Mi_sm <- tryCatch(
        solve(M_sm),
        error = function(e) {
          stop(
            "Prior_Setup: cannot invert Sigma_0 + (X'WX)^{-1} (S_marg_new). ",
            conditionMessage(e),
            call. = FALSE
          )
        }
      )
      quad_sm <- as.numeric(crossprod(dlt_sm, Mi_sm %*% dlt_sm))
      if (!is.finite(quad_sm) || quad_sm < 0) {
        stop(
          "Prior_Setup: S_marg quadratic form is not finite or nonnegative.",
          call. = FALSE
        )
      }
      S_marg_new <- rss_weighted_stored + quad_sm
      S_marg <- S_marg_new
      d_vcov <- glm_summary$dispersion
      if (is.finite(d_vcov) && d_vcov > 0) {
        Sigma_0_vcov <- Sigma_pre_nm / d_vcov
        M_vc <- Sigma_0_vcov + Ginv_sm
        Mi_vc <- tryCatch(solve(M_vc), error = function(e) NULL)
        if (!is.null(Mi_vc)) {
          quad_vc <- as.numeric(crossprod(dlt_sm, Mi_vc %*% dlt_sm))
          if (is.finite(quad_vc) && quad_vc >= 0) {
            S_marg_sigma0_vcov <- rss_weighted_stored + quad_vc
          }
        }
      }
      if (length(pwt) == 1L && is.finite(pwt)) {
        quad_s_sm <- as.numeric(pwt * crossprod(dlt_sm, Gm_sm %*% dlt_sm))
        if (is.finite(quad_s_sm) && quad_s_sm >= 0) {
          S_marg_scalar_zellner <- rss_weighted_stored + quad_s_sm
        }
      }
    }
  }

  ## Nelder–Mead Post_mean joint optimization: entire block below is commented out (not deleted).
#   ## --- Gaussian: Post_mean - joint fixed point for d and prior scale lambda ----------
#   ## Sigma_g = ((1-pwt)/pwt)*V0.  Search (d, lambda) with returned Sigma = Sigma_g*lambda so that, using the
#   ## same fragment as dNormal_Gamma (rNormal_reg.wfit with P = d*Sigma^{-1}, full w):
#   ##   (i) posterior mean of beta equals (1-pwt)*beta_MLE + pwt*mu (per coefficient, where finite);
#   ##   (ii) d = E[dispersion|y] = (shape*d + S/2)/(a_post-1), a_post = shape + n_eff/2.
#   ## The Nelder–Mead optimization is retained (not replaced by a shortcut) so results stay available
#   ## for validation against direct or closed-form recipes (temporary \code{message} when run).
#   if (identical(family$family, "gaussian") &&
#       identical(disp_type, "Post_mean") &&
#       !is.null(shape) && is.finite(shape) && shape > 0 &&
#       is.finite(n_effective) && n_effective > 0) {
#     Sigma_g <- Sigma
#     sigma_mle_sq <- glm_summary$dispersion
#     mu_fp <- as.numeric(mu)
#     w_fp <- if (is.null(weights)) rep(1, n_obs) else as.numeric(weights)
#     off_fp <- offset
#     a_post_fp <- shape + n_effective / 2
# 
#     if (a_post_fp > 1 && is.finite(dispersion) && dispersion > 0 &&
#         is.finite(sigma_mle_sq) && sigma_mle_sq > 0) {
#       d_fp <- dispersion
#       lam_start <- d_fp^2 / sigma_mle_sq
#       if (!is.finite(lam_start) || lam_start <= 0) {
#         lam_start <- 1
#       }
# 
#       ng_wfit_dlam <- function(d, lam) {
#         if (!is.finite(d) || d <= 0 || !is.finite(lam) || lam <= 0) {
#           return(NULL)
#         }
#         Sig <- Sigma_g * lam
#         R <- tryCatch(chol(Sig), error = function(e) NULL)
#         if (is.null(R)) {
#           return(NULL)
#         }
#         Pinv <- chol2inv(R)
#         Pinv <- 0.5 * (Pinv + t(Pinv))
#         P_ng <- d * Pinv
#         tryCatch(
#           rNormal_reg.wfit(
#             x = X, y = Y, P = P_ng, mu = mu_fp, w = w_fp,
#             offset = off_fp
#           ),
#           error = function(e) NULL
#         )
#       }
# 
#       ## dNormal(mu, Sigma, dispersion=d): prior precision inv(Sigma); likelihood uses w/d.
#       post_coef_dNormal_dlam <- function(d, lam) {
#         if (!is.finite(d) || d <= 0 || !is.finite(lam) || lam <= 0) {
#           return(rep(NA_real_, length(mu_fp)))
#         }
#         Sig <- Sigma_g * lam
#         R <- tryCatch(chol(Sig), error = function(e) NULL)
#         if (is.null(R)) {
#           return(rep(NA_real_, length(mu_fp)))
#         }
#         Pinv <- chol2inv(R)
#         Pinv <- 0.5 * (Pinv + t(Pinv))
#         fit <- rNormal_reg.wfit(
#           x = X, y = Y, P = Pinv, mu = mu_fp, w = w_fp / d,
#           offset = off_fp
#         )
#         as.numeric(fit$Btilde)
#       }
# 
#       cmle <- coef(glm_full)
#       mle_fp <- vapply(
#         var_names,
#         function(nm) {
#           if (!is.null(names(cmle)) && nm %in% names(cmle)) {
#             v <- unname(cmle[nm])
#             if (length(v) == 1L && is.finite(v)) v else NA_real_
#           } else {
#             NA_real_
#           }
#         },
#         NA_real_
#       )
#       blend_fp <- (1 - pwt) * mle_fp + pwt * mu_fp
# 
#       Sigma_start <- Sigma_g * lam_start
#       fit_ng0 <- ng_wfit_dlam(d_fp, lam_start)
#       coef_ng0 <- if (!is.null(fit_ng0)) {
#         as.numeric(fit_ng0$Btilde)
#       } else {
#         rep(NA_real_, length(mu_fp))
#       }
#       coef_n0 <- post_coef_dNormal_dlam(d_fp, lam_start)
# 
#       ## Post_mean console diagnostics (commented out; delete when no longer needed)
#       # cat(
#       #   "\nPrior_Setup Post_mean diagnostics (before Nelder-Mead):\n",
#       #   "  prior mean (mu):\n",
#       #   sep = ""
#       # )
#       # print(stats::setNames(mu_fp, var_names))
#       # cat("  MLE (weighted glm_full coefficients):\n")
#       # print(stats::setNames(mle_fp, var_names))
#       # if (length(pwt) == 1L) {
#       #   cat(
#       #     "  (1-pwt)*MLE + pwt*mu  with pwt = ", format(pwt, digits = 10), ":\n",
#       #     sep = ""
#       #   )
#       # } else {
#       #   cat("  (1-pwt)*MLE + pwt*mu  (vector pwt per coefficient):\n")
#       # }
#       # print(stats::setNames(blend_fp, var_names))
#       # cat(
#       #   "  starting dispersion d0 = ", format(d_fp, digits = 10),
#       #   "  lambda0 = ", format(lam_start, digits = 10), "\n",
#       #   "  starting Sigma = Sigma_g * lambda0 (rounded):\n",
#       #   sep = ""
#       # )
#       # print(round(Sigma_start, 6))
#       # cat("  posterior mean coefficients - dNormal_Gamma @ (d0, lambda0):\n")
#       # print(stats::setNames(coef_ng0, var_names))
#       # cat("  posterior mean coefficients - dNormal @ (d0, lambda0):\n")
#       # print(stats::setNames(coef_n0, var_names))
# 
#       Q_joint_fp <- function(uv) {
#         d <- exp(uv[1L])
#         lam <- exp(uv[2L])
#         fit <- ng_wfit_dlam(d, lam)
#         if (is.null(fit)) {
#           return(.Machine$double.xmax^0.25)
#         }
#         post <- as.numeric(fit$Btilde)
#         if (length(post) != length(blend_fp)) {
#           return(.Machine$double.xmax^0.25)
#         }
#         idx <- is.finite(blend_fp) & is.finite(post)
#         if (!any(idx)) {
#           return(.Machine$double.xmax^0.25)
#         }
#         diff_b <- post[idx] - blend_fp[idx]
#         w_beta <- 1 / max(mean(blend_fp[idx]^2, na.rm = TRUE), 1e-12)
#         term_beta <- w_beta * mean(diff_b^2)
#         rate_m <- (n_prior / 2) * d
#         S_m <- as.numeric(fit$S)
#         if (length(S_m) != 1L || !is.finite(S_m)) {
#           return(.Machine$double.xmax^0.25)
#         }
#         b_m <- rate_m + 0.5 * S_m
#         m_disp <- b_m / (a_post_fp - 1)
#         if (!is.finite(m_disp) || m_disp <= 0) {
#           return(.Machine$double.xmax^0.25)
#         }
#         term_disp <- (log(m_disp) - log(d))^2
#         term_disp + term_beta
#       }
# 
#       opt_fp <- tryCatch(
#         suppressWarnings(
#           stats::optim(
#             par = c(log(d_fp), log(lam_start)),
#             fn = Q_joint_fp,
#             method = "Nelder-Mead",
#             control = list(maxit = 3000L, reltol = 1e-10)
#           )
#         ),
#         error = function(e) NULL
#       )
#       ok <- !is.null(opt_fp) && length(opt_fp$par) == 2L &&
#         is.finite(opt_fp$value) && opt_fp$value < .Machine$double.xmax^0.2
#       nm_dispersion <- NA_real_
#       nm_lambda <- NA_real_
#       post_mean_return_applied <- FALSE
#       if (ok) {
#         u_opt <- as.numeric(opt_fp$par)
#         d_nm <- exp(u_opt[1L])
#         lam_nm <- exp(u_opt[2L])
#         nm_dispersion <- d_nm
#         nm_lambda <- lam_nm
#         if (is.finite(d_nm) && d_nm > 0 && is.finite(lam_nm) && lam_nm > 0) {
#           dispersion <- d_nm
#           rate <- (n_prior / 2) * d_nm
#           Sigma <- Sigma_g * lam_nm
#           fit_ng1 <- ng_wfit_dlam(d_nm, lam_nm)
#           coef_ng1 <- if (!is.null(fit_ng1)) {
#             as.numeric(fit_ng1$Btilde)
#           } else {
#             rep(NA_real_, length(mu_fp))
#           }
#           coef_n1 <- post_coef_dNormal_dlam(d_nm, lam_nm)
#           post_mean_return_applied <- TRUE
#           ## Post_mean optim result diagnostics (commented out; delete when no longer needed)
#           # cat(
#           #   "  final dispersion d* = ", format(d_nm, digits = 10),
#           #   "  lambda* = ", format(lam_nm, digits = 10),
#           #   "  (optim convergence = ", opt_fp$convergence, ")\n",
#           #   "  joint objective Q = ", format(opt_fp$value, digits = 10), "\n",
#           #   "  final Sigma = Sigma_g * lambda* (rounded):\n",
#           #   sep = ""
#           # )
#           # print(round(Sigma, 6))
#           # cat("  posterior mean coefficients - dNormal_Gamma @ (d*, lambda*):\n")
#           # print(stats::setNames(coef_ng1, var_names))
#           # cat("  posterior mean coefficients - dNormal @ (d*, lambda*):\n")
#           # print(stats::setNames(coef_n1, var_names))
#           # cat("\n")
#         }
#       } else {
#         ## Nelder-Mead failure message (commented out; delete when no longer needed)
#         # cat(
#         #   "  Nelder-Mead did not return a usable optimum; Sigma and dispersion unchanged.\n\n"
#         # )
#       }
# 
#       ## Temporary Nelder–Mead diagnostics (not returned; delete when no longer needed).
#       post_mean_diag <- list(
#         ok = ok,
#         return_update_applied = post_mean_return_applied,
#         starter_dispersion = d_fp,
#         starter_lambda = lam_start,
#         dispersion = nm_dispersion,
#         lambda = nm_lambda,
#         convergence = if (!is.null(opt_fp)) opt_fp$convergence else NA_integer_,
#         objective_Q = if (!is.null(opt_fp)) opt_fp$value else NA_real_,
#         par_log = if (!is.null(opt_fp)) as.numeric(opt_fp$par) else rep(NA_real_, 2L)
#       )
#       nm <- post_mean_diag
#       message(
#         "Prior_Setup temporary Post_mean / Nelder diagnostics:\n",
#         paste0(
#           c(
#             paste0("  ok = ", nm$ok),
#             paste0("  return_update_applied = ", nm$return_update_applied),
#             paste0("  starter_dispersion = ", format(nm$starter_dispersion, digits = 10)),
#             paste0("  starter_lambda     = ", format(nm$starter_lambda, digits = 10)),
#             paste0("  dispersion (opt)   = ", format(nm$dispersion, digits = 10)),
#             paste0("  lambda (opt)       = ", format(nm$lambda, digits = 10)),
#             paste0("  objective_Q        = ", format(nm$objective_Q, digits = 10)),
#             paste0("  convergence        = ", nm$convergence)
#           ),
#           collapse = "\n"
#         )
#       )
#     }
#   }

  ## S_marg_post_mean: same RSS + quadratic as S_marg_new, but Sigma_0 = Sigma/dispersion **after**
  ## Post_mean (Nelder--Mead may have updated Sigma and dispersion).
  if (identical(family$family, "gaussian") &&
      is.finite(rss_weighted_stored) && rss_weighted_stored > 0 &&
      !is.null(dispersion) && is.finite(dispersion) && dispersion > 0) {
    bh_pm <- coef(glm_full)
    if (length(bh_pm) == nvar && all(is.finite(bh_pm))) {
      w_pm <- if (is.null(weights)) rep(1, n_obs) else as.numeric(weights)
      XtW_pm <- sweep(X, 1, w_pm, `*`)
      Gm_pm <- crossprod(XtW_pm, X)
      dlt_pm <- matrix(bh_pm, ncol = 1L) - matrix(as.numeric(mu), ncol = 1L)
      Sigma_0_pm <- Sigma / dispersion
      Ginv_pm <- tryCatch(
        solve(Gm_pm),
        error = function(e) {
          stop(
            "Prior_Setup: cannot invert weighted Gram matrix X'WX (S_marg_post_mean). ",
            conditionMessage(e),
            call. = FALSE
          )
        }
      )
      M_pm <- Sigma_0_pm + Ginv_pm
      Mi_pm <- tryCatch(
        solve(M_pm),
        error = function(e) {
          stop(
            "Prior_Setup: cannot invert Sigma/dispersion + (X'WX)^{-1} (S_marg_post_mean). ",
            conditionMessage(e),
            call. = FALSE
          )
        }
      )
      quad_pm <- as.numeric(crossprod(dlt_pm, Mi_pm %*% dlt_pm))
      if (!is.finite(quad_pm) || quad_pm < 0) {
        stop(
          "Prior_Setup: S_marg_post_mean quadratic form is not finite or nonnegative.",
          call. = FALSE
        )
      }
      S_marg_post_mean <- rss_weighted_stored + quad_pm
    }
  }

  ## Temporary S_marg printed diagnostics: entire block below commented out (not deleted).
#   ## Temporary S_marg diagnostics (Gaussian only; not returned; delete when no longer needed).
#   if (identical(family$family, "gaussian")) {
#     lines_sm <- character(0)
#     if (is.finite(S_marg_new)) {
#       lines_sm <- c(
#         lines_sm,
#         paste0("  S_marg_new (Sigma_pre/d_OLS, pre-Post_mean): ", format(S_marg_new, digits = 10))
#       )
#     }
#     if (is.finite(S_marg_sigma0_vcov)) {
#       lines_sm <- c(
#         lines_sm,
#         paste0(
#           "  S_marg_sigma0_vcov (Sigma_pre/summary(glm)$dispersion): ",
#           format(S_marg_sigma0_vcov, digits = 10)
#         )
#       )
#     }
#     if (is.finite(S_marg_new) && is.finite(S_marg_sigma0_vcov)) {
#       dvc <- abs(S_marg_new - S_marg_sigma0_vcov)
#       rvc <- if (abs(S_marg_new) > 0) dvc / abs(S_marg_new) else dvc
#       lines_sm <- c(
#         lines_sm,
#         paste0(
#           "  |S_marg_new - S_marg_sigma0_vcov| = ", format(dvc, digits = 10),
#           ", rel. diff = ", format(rvc, digits = 10)
#         )
#       )
#     }
#     if (is.finite(S_marg_post_mean)) {
#       lines_sm <- c(
#         lines_sm,
#         paste0(
#           "  S_marg_post_mean (Sigma/dispersion after Post_mean): ",
#           format(S_marg_post_mean, digits = 10)
#         )
#       )
#     }
#     if (is.finite(S_marg_new) && is.finite(S_marg_post_mean)) {
#       dpm <- abs(S_marg_new - S_marg_post_mean)
#       rpm <- if (abs(S_marg_new) > 0) dpm / abs(S_marg_new) else dpm
#       lines_sm <- c(
#         lines_sm,
#         paste0(
#           "  |S_marg_new - S_marg_post_mean| = ", format(dpm, digits = 10),
#           ", rel. diff = ", format(rpm, digits = 10)
#         )
#       )
#     }
#     if (length(lines_sm) > 0L) {
#       message("Prior_Setup temporary S_marg diagnostics:\n", paste(lines_sm, collapse = "\n"))
#     }
#   }

  ## Gaussian with shape but iteration did not run (e.g. chol failed): fall back rate.
  if (identical(family$family, "gaussian") && !is.null(shape) && is.null(rate) &&
      !is.null(dispersion_for_shape_rate) && !is.null(n_prior) && length(n_prior) == 1L) {
    rate <- dispersion_for_shape_rate * (n_prior / 2)
    if (!is.finite(rate) || rate <= 0) {
      stop("Computed rate must be strictly positive.")
    }
  }

  ## --- Marginal sum of squares S_marg (Gaussian): S_marg matches S_marg_new; fallback below --------
  ## Fallback below only if the pre-Post_mean computation did not yield a finite value.
  ## General: RSS + (beta_hat - mu)^T (Sigma_0 + G^{-1})^{-1} (beta_hat - mu),
  ##   Sigma_0 = Sigma/d_OLS with d_OLS = RSS_w/(n_effective-2) (= dispersion_classical).
  ## Scalar-pwt Zellner check: RSS + pwt * (beta_hat - mu)^T G (beta_hat - mu).
  dispersion_marginal <- NA_real_
  dispersion_nelder_mead <- NA_real_
  b_0_rate_dispersion <- NA_real_
  b_0_S_marg_formula <- NA_real_
  E_phi_sigma2_special <- NA_real_
  if (identical(family$family, "gaussian") &&
      is.finite(rss_weighted_stored) && rss_weighted_stored > 0 &&
      !is.null(dispersion) && is.finite(dispersion) && dispersion > 0) {
    dispersion_nelder_mead <- dispersion
    if (!is.null(n_prior) && length(n_prior) == 1L &&
        is.finite(n_prior) && n_prior > 0 &&
        is.finite(dispersion_nelder_mead)) {
      b_0_rate_dispersion <- (n_prior / 2) * dispersion_nelder_mead
    }
    bh <- coef(glm_full)
    if (!is.null(shape) && is.finite(shape) && shape > 0 &&
        (length(bh) != nvar || !all(is.finite(bh)))) {
      stop(
        "Prior_Setup: Gaussian Normal-Gamma prior requires finite GLM coefficients for every coefficient.",
        call. = FALSE
      )
    }
    if (length(bh) == nvar && all(is.finite(bh))) {
      w_g <- if (is.null(weights)) rep(1, n_obs) else as.numeric(weights)
      XtW <- sweep(X, 1, w_g, `*`)
      Gm <- crossprod(XtW, X)
      dlt <- matrix(bh, ncol = 1L) - matrix(as.numeric(mu), ncol = 1L)
      Ginv <- tryCatch(
        solve(Gm),
        error = function(e) {
          stop(
            "Prior_Setup: cannot invert weighted Gram matrix X'WX. ",
            conditionMessage(e),
            call. = FALSE
          )
        }
      )
      if (!is.finite(S_marg) &&
          is.finite(dispersion_classical) && dispersion_classical > 0) {
        Sigma_0 <- Sigma_pre_nm / dispersion_classical
        M <- Sigma_0 + Ginv
        Mi <- tryCatch(
          solve(M),
          error = function(e) {
            stop(
              "Prior_Setup: cannot invert Sigma_0 + (X'WX)^{-1} (S_marg fallback). ",
              conditionMessage(e),
              call. = FALSE
            )
          }
        )
        quad <- as.numeric(crossprod(dlt, Mi %*% dlt))
        if (!is.finite(quad) || quad < 0) {
          stop(
            "Prior_Setup: S_marg fallback quadratic form is not finite or nonnegative.",
            call. = FALSE
          )
        }
        S_marg <- rss_weighted_stored + quad
      }
      if (!is.finite(S_marg_scalar_zellner) && length(pwt) == 1L && is.finite(pwt)) {
        quad_s <- as.numeric(pwt * crossprod(dlt, Gm %*% dlt))
        if (is.finite(quad_s) && quad_s >= 0) {
          S_marg_scalar_zellner <- rss_weighted_stored + quad_s
        }
      }
      ## Same denominator as classical Gaussian dispersion (RSS_w/(n_effective-2)).
      if (is.finite(S_marg) && is.finite(n_effective) && n_effective > 2) {
        dispersion_marginal <- S_marg / (n_effective - 2)
      }
      if (!is.null(n_prior) && length(n_prior) == 1L &&
          is.finite(n_prior) && n_prior > 0 &&
          is.finite(n_effective) && n_effective > 0 &&
          is.finite(S_marg)) {
        b_0_S_marg_formula <- 0.5 * (n_prior / n_effective) * S_marg
        den_phi <- n_prior + n_effective - 2
        if (is.finite(den_phi) && den_phi > 0) {
          E_phi_sigma2_special <-
            S_marg * (n_effective + n_prior) / n_effective / den_phi
        }
      }
      ## S_marg validation message (commented out; delete when no longer needed)
      # if (length(pwt) == 1L && is.finite(S_marg) && is.finite(S_marg_scalar_zellner)) {
      #   diff_sm <- abs(S_marg - S_marg_scalar_zellner)
      #   rel_sm <- if (abs(S_marg) > 0) diff_sm / abs(S_marg) else diff_sm
      #   message(
      #     "Prior_Setup S_marg validation (scalar pwt): general = ",
      #     format(S_marg, digits = 10),
      #     ", Zellner scalar = ",
      #     format(S_marg_scalar_zellner, digits = 10),
      #     ", |diff| = ",
      #     format(diff_sm, digits = 10),
      #     ", rel_err = ",
      #     format(rel_sm, digits = 10),
      #     "."
      #   )
      # }
      ## Compare three dispersion notions and Zellner Sigmas (commented out; delete when no longer needed)
      # if (length(pwt) == 1L && is.finite(pwt) && pwt > 0 && pwt < 1 &&
      #     !is.null(Ginv)) {
      #   cat(
      #     "\nPrior_Setup dispersion comparison (Gaussian, scalar pwt):\n",
      #     "  classical RSS_w/(n_effective-2)     = ",
      #     format(dispersion_classical, digits = 10),
      #     "\n  marginal S_marg/(n_effective-2)   = ",
      #     format(dispersion_marginal, digits = 10),
      #     "\n  Nelder-Mead / returned (pre-marg) = ",
      #     format(dispersion_nelder_mead, digits = 10),
      #     "\n",
      #     sep = ""
      #   )
      #   gfac <- (1 - pwt) / pwt
      #   Sigma_classical_disp <- gfac * dispersion_classical * Ginv
      #   Sigma_marginal_disp <- gfac * dispersion_marginal * Ginv
      #   dimnames(Sigma_classical_disp) <- list(var_names, var_names)
      #   dimnames(Sigma_marginal_disp) <- list(var_names, var_names)
      #   cat("  Sigma from classical d ((1-pwt)/pwt * d * G^{-1}), rounded:\n")
      #   print(round(Sigma_classical_disp, 6))
      #   cat("  Sigma from marginal d ((1-pwt)/pwt * (S_marg/(n_w-2)) * G^{-1}), rounded:\n")
      #   print(round(Sigma_marginal_disp, 6))
      #   cat("  Sigma from Nelder-Mead / Post_mean (current matrix), rounded:\n")
      #   print(round(Sigma, 6))
      #   if (is.finite(b_0_rate_dispersion) || is.finite(b_0_S_marg_formula)) {
      #     cat(
      #       "  Gamma prior rate b_0:\n",
      #       "    (n_prior/2)*d_Nelder (standard)     = ",
      #       format(b_0_rate_dispersion, digits = 10),
      #       "\n",
      #       "    (1/2)(n_prior/n_w)*S_marg (special) = ",
      #       format(b_0_S_marg_formula, digits = 10),
      #       "\n",
      #       sep = ""
      #     )
      #   }
      #   if (is.finite(E_phi_sigma2_special)) {
      #     cat(
      #       "  E[sigma^2]=E[tau^-1] (special b_0; b_n/(a_n-1)) = ",
      #       format(E_phi_sigma2_special, digits = 10),
      #       "\n",
      #       "    = S_marg*(n_w+n_prior)/n_w/(n_prior+n_w-2)\n",
      #       sep = ""
      #     )
      #     ## Cov(beta|y) = (b_n/(a_n-1)) V_n with V_n = (n_w/(n_prior+n_w)) G^{-1}
      #     ## (Zellner scalar pwt); algebraically (S_marg/(n_prior+n_w-2)) G^{-1}.
      #     v_n_fac <- n_effective / (n_prior + n_effective)
      #     Cov_beta_y_special <- E_phi_sigma2_special * v_n_fac * Ginv
      #     dimnames(Cov_beta_y_special) <- list(var_names, var_names)
      #     cat(
      #       "  Cov(beta|y) (special b_0; (b_n/(a_n-1))*V_n), rounded:\n",
      #       "    V_n = (n_w/(n_prior+n_w))*G^{-1}  [same as (S_marg/(n_prior+n_w-2))*G^{-1}]\n",
      #       sep = ""
      #     )
      #     print(round(Cov_beta_y_special, 6))
      #   }
      #   cat("\n")
      # } else if (identical(family$family, "gaussian") && is.finite(dispersion_marginal)) {
      #   cat(
      #     "\nPrior_Setup dispersion comparison (Gaussian, vector pwt):\n",
      #     "  classical RSS_w/(n_effective-2) = ",
      #     format(dispersion_classical, digits = 10),
      #     "\n  marginal S_marg/(n_effective-2) = ",
      #     format(dispersion_marginal, digits = 10),
      #     "\n  Nelder-Mead / current d       = ",
      #     format(dispersion_nelder_mead, digits = 10),
      #     "\n  (Zellner Sigma triple skipped: pwt not scalar)\n",
      #     sep = ""
      #   )
      #   if (is.finite(b_0_rate_dispersion) || is.finite(b_0_S_marg_formula)) {
      #     cat(
      #       "  Gamma prior rate b_0:\n",
      #       "    (n_prior/2)*d_Nelder (standard)     = ",
      #       format(b_0_rate_dispersion, digits = 10),
      #       "\n",
      #       "    (1/2)(n_prior/n_w)*S_marg (special) = ",
      #       format(b_0_S_marg_formula, digits = 10),
      #       "\n",
      #       sep = ""
      #     )
      #   }
      #   if (is.finite(E_phi_sigma2_special)) {
      #     cat(
      #       "  E[sigma^2]=E[tau^-1] (special b_0; b_n/(a_n-1)) = ",
      #       format(E_phi_sigma2_special, digits = 10),
      #       "\n",
      #       "    = S_marg*(n_w+n_prior)/n_w/(n_prior+n_w-2)\n",
      #       sep = ""
      #     )
      #   }
      #   cat("\n\n")
      # }
      ## Remap returned prior to the special S_marg path:
      ## - dispersion: E[sigma^2|y] from b_n/(a_n-1)
      ## - Sigma: PRIOR covariance (not posterior), Zellner form
      ##          ((1-pwt)/pwt) * dispersion * G^{-1} = (n_w/n_prior) * dispersion * G^{-1}
      ## - rate: b_0 = (1/2)(n_prior/n_w) S_marg
      if (!is.null(shape) && is.finite(shape) && shape > 0) {
        if (!is.finite(S_marg)) {
          stop(
            "Prior_Setup: S_marg is not finite; cannot remap prior to the marginal sum-of-squares calibration.",
            call. = FALSE
          )
        }
        den_phi_remap <- n_prior + n_effective - 2
        if (!is.finite(den_phi_remap) || den_phi_remap <= 0) {
          stop(
            "Prior_Setup: require n_prior + n_effective > 2 for S_marg remap.",
            call. = FALSE
          )
        }
        if (!is.finite(E_phi_sigma2_special) || E_phi_sigma2_special <= 0) {
          stop(
            "Prior_Setup: E[sigma^2|y] for the S_marg path is missing or not positive.",
            call. = FALSE
          )
        }
        if (!is.finite(b_0_S_marg_formula) || b_0_S_marg_formula <= 0) {
          stop(
            "Prior_Setup: Gamma rate term b_0 from S_marg is missing or not positive.",
            call. = FALSE
          )
        }
        dispersion <- E_phi_sigma2_special
        Sigma <- (n_effective / n_prior) * dispersion * Ginv
        rownames(Sigma) <- var_names
        colnames(Sigma) <- var_names
        rate <- b_0_S_marg_formula
      }
    }
  }

  coefficients_mle <- coef(glm_full)
  coefficients <- coefficients_mle
  ## Default returned coefficients: closed-form posterior mean blend when available.
  if (identical(family$family, "gaussian") &&
      length(coefficients_mle) == nvar &&
      !is.null(mu) && length(mu) == nvar &&
      all(is.finite(as.numeric(mu)))) {
    mle_fp <- vapply(
      var_names,
      function(nm) {
        if (!is.null(names(coefficients_mle)) && nm %in% names(coefficients_mle)) {
          v <- unname(coefficients_mle[nm])
          if (length(v) == 1L && is.finite(v)) v else NA_real_
        } else {
          NA_real_
        }
      },
      NA_real_
    )
    mu_fp <- as.numeric(mu)
    if (length(pwt) == 1L && is.finite(pwt)) {
      coefficients <- (1 - pwt) * mle_fp + pwt * mu_fp
      names(coefficients) <- var_names
    } else if (length(pwt) == nvar && all(is.finite(pwt))) {
      coefficients <- (1 - pwt) * mle_fp + pwt * mu_fp
      names(coefficients) <- var_names
    }
  }
  
  prior_list <- list(
    mu = mu,
    Sigma = Sigma,
    dispersion = dispersion,
    shape = shape,
    rate = rate,
    coefficients = coefficients,
    model = mf,
    x = x,
    y = Y,
    call = call,
    PriorSettings = list(
      pwt = pwt,
      n_prior = n_prior,
      intercept_source = intercept_source,
      effects_source = effects_source,
      ## For now retain n_likelihood for backward compatibility
      n_likelihood = n_likelihood,
      n_effective = n_effective
    )
  )
  
  class(prior_list) <- "PriorSetup"
  return(prior_list)
  
}



#' @export
#' @method print PriorSetup
#' @rdname Prior_Setup

print.PriorSetup <- function(x, ...) {

  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    
  settings <- x$PriorSettings
  
  if (!is.null(settings$pwt) && length(settings$pwt) == 1L) {
    g <- (1 - settings$pwt) / settings$pwt
    cat("Setting up a Zellner g-type prior: \n")
    cat("  pwt =", round(settings$pwt, 4), "\n")
    cat("  g   = (1 - pwt)/pwt =", round(g, 4), "\n\n")
  }
  
  if (!is.null(settings$n_prior) && !is.null(settings$n_likelihood)) {
    if (!is.null(settings$pwt) && length(settings$pwt) > 1L) {
      cat("Note: n_prior was provided to Prior_Setup (scalar prior sample size for precision / calibration; pwt is per coefficient):\n")
    } else {
      n_eff_print <- settings$n_effective
      if (is.null(n_eff_print)) {
        n_eff_print <- settings$n_likelihood
      }
      pwt_s <- settings$pwt
      n_prior_implies_pwt <- FALSE
      if (length(pwt_s) == 1L && is.finite(pwt_s) && pwt_s > 0 && pwt_s < 1 &&
          is.finite(n_eff_print) && n_eff_print > 0 &&
          is.finite(settings$n_prior)) {
        n_prior_from_pwt <- (pwt_s / (1 - pwt_s)) * n_eff_print
        np_stored <- as.numeric(settings$n_prior)
        scale_np <- max(abs(np_stored), abs(n_prior_from_pwt), 1e-12)
        n_prior_implies_pwt <- abs(np_stored - n_prior_from_pwt) <= 1e-10 * scale_np
      }
      if (n_prior_implies_pwt) {
        cat("Note: n_prior was computed as (pwt / (1 - pwt)) * n_likelihood: \n")
      } else {
        cat("Note: n_prior was provided to Prior_Setup (scalar prior sample size for precision / calibration):\n")
      }
    }
    cat("  n_prior      =", round(settings$n_prior, 4), "\n")
    cat("  n_likelihood =", round(settings$n_likelihood, 4), "\n\n")
  }

  if (!is.null(settings$pwt) && length(settings$pwt) > 1L) {
    cat("Note: Differential prior weights (pwt) were specified per coefficient.\n\n")
  }
  
  
  cat("Prior Setup Summary\n")
  cat("====================\n\n")
  
  # Check for Zellner g-prior structure
  Sigma <- x$Sigma
  mu <- x$mu
  var_names <- rownames(mu)
  nvar <- length(var_names)
  
  # Extract diagonal SDs
  prior_sd <- sqrt(diag(Sigma))
  
  # Always compute prior correlation matrix
  prior_cor <- cov2cor(Sigma)
  
  # Extract pwt vector for display
  if (!is.null(settings$pwt)) {
    if (length(settings$pwt) == 1L) {
      pwt_vec <- rep(settings$pwt, nvar)
    } else if (length(settings$pwt) == nvar) {
      pwt_vec <- settings$pwt
    } else {
      warning("Length of pwt does not match number of variables; skipping pwt column.")
      pwt_vec <- rep(NA_real_, nvar)
    }
  } else {
    pwt_vec <- rep(NA_real_, nvar)
  }
  

  # Compute 95% intervals
  z <- qnorm(0.975)
  lower <- mu[, 1] - z * prior_sd
  upper <- mu[, 1] + z * prior_sd
  

  # Build output table

  out <- data.frame(
    Prior.Mean = round(mu[, 1], 6),
    Prior.SD   = round(prior_sd, 6),
    CI.Lower   = round(lower, 6),
    CI.Upper   = round(upper, 6),
    pwt        = round(pwt_vec, 6)
  )
  
  # Print table
  cat("Prior Estimates with 95% Confidence Intervals\n")
  invisible(print(out))
  
  if (nvar <= 10) {
    cat("\nPrior Correlation Matrix\n")
    invisible(print(round(prior_cor, 4)))
  }
  
  # Optional: print dispersion
  if (!is.null(x$dispersion)) {
    cat("\nConditional Dispersion (Gaussian family): ", round(x$dispersion, 4), "\n\n")
  }

  if (!is.null(x$shape) && !is.null(x$rate)) {
    cat("Gamma Prior on Residual Precision:\n")
    sdf <- x$PriorSettings$shape_df
    if (!is.null(sdf)) {
      cat("  shape_df =", dQuote(sdf, FALSE), "\n")
    }
    for (nm in c("shape", "rate")) {
      z <- x[[nm]]
      zch <- if (length(z) != 1L || !is.numeric(z)) {
        as.character(z)
      } else if (!is.finite(z)) {
        as.character(z)
      } else {
        az <- abs(z)
        if (az == 0) {
          "0"
        } else if (az < 1e-4) {
          format(z, scientific = TRUE, digits = 6L)
        } else {
          as.character(round(z, 4L))
        }
      }
      cat("  ", nm, " = ", zch, "\n", sep = "")
    }
    cat(
      "  Expected precision (inverse variance) =",
      format(signif(x$shape / x$rate, 6), scientific = TRUE, digits = 7),
      "which implies 1/Expected precision  =",
      format(signif(x$rate / x$shape, 6), scientific = TRUE, digits = 7),
      "\n\n"
    )
    cat("  Applicable to gaussian models with compound pfamilies (e.g., dNormal_Gamma, dIndependent_Normal_Gamma),\n")
    cat("  as well as for Gamma regression, quasipoisson, and quasibinomial models.\n\n")
  }
  
  if (is.null(x$shape) && is.null(x$rate) && !is.null(x$dispersion)) {
    cat("Note: Gaussian family detected, but shape/rate parameters were not computed.\n")
    cat("This may occur if n_prior is not scalar.\n\n")
  }
  
  invisible(x)
}



#' Checks for Prior-data conflicts
#'
#' Checks if the credible intervals for the prior overlap with the implied confidence intervals
#' from the classical model (obtained via \code{\link[stats]{glm}}). The approach relates to
#' prior-data conflict checks \insertCite{EvansMoshonov2006}{glmbayes}.
#'
#' @param level the confidence level at which the Prior-data conflict should be checked.
#' @inheritParams glmb
#' @return A vector where each item provided the ratio of the absolue value for the difference between the 
#' prior and maximum likelihood estimate divided by the length of the sum of half of the two intervals 
#' (where normality is assumed)
#' @seealso \code{\link{Prior_Setup}}, \code{\link{glmb}}; see \insertCite{glmbayesChapter03}{glmbayes} for prior tailoring.
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @family prior
#' @example inst/examples/Ex_Prior_Check.R
#' @export
#' @rdname Prior_Check
#' @order 1

Prior_Check<-function(formula,family,pfamily,level=0.95,data=NULL, weights, subset,na.action, 
                      start = NULL, etastart, mustart, offset ,control = list(...) , model = TRUE, 
                      method = "glm.fit",x = FALSE, y = TRUE, contrasts = NULL, ...){
  
  pf=pfamily
  prior_list=pfamily$prior_list
  
  ## For now, the below is really only correct for the dNormal pfamily
  
  mu=prior_list$mu
  Sigma=prior_list$Sigma
  
  
  object=glm(formula=formula,family=family,data=data)
  
  Like_est=object$coefficients
  Like_std=summary(object)$coefficients[,2]
  
  if(is.null(mu)){
    print("No Prior mean vector provided. Variables with needed Priors are:")
    print(names(Like_est))
    names(Like_est)    
    
  }
  
  
  if(level<=0.5) stop("level must be greater than 0.5")
  
  Sigma=as.matrix(Sigma)
  Prior_std=sqrt(diag(Sigma))
  
  print("Variables in the Model Are:")
  print(names(Like_est))
  std_dev_sum=qnorm(level)*(Prior_std+Like_std)
  
  abs_ratio=matrix(rep(0,length(Like_est),nrow=length(Like_est),ncol=1))
  abs_diff=abs(mu-Like_est)
  abs_ratio[1:length(Like_est),]=abs_diff/std_dev_sum
  
  rownames(abs_ratio)=names(Like_est)
  colnames(abs_ratio)=c("abs_ratio")
  max_abs_ratio=max(abs_ratio)
  
  if(max_abs_ratio>1) {
    print("At least one of the maximum likelihood estimates appears to be inconsistent with the prior")
  }
  
  else{
    print("The maximum likelihood estimates for all coefficients appear to be roughly consistent with the prior.")
  }
  return(abs_ratio)
  
}




