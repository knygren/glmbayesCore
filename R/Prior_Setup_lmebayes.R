#' Prior setup for the two-block Gibbs lmebayes sampler
#'
#' Calibrates priors for the level-2 fixed effects (\code{fixef}) of a
#' hierarchical mixed model using the reference \code{lmer}/\code{glmer} fits
#' on \strong{all} groups (from \code{\link{model_setup}}).  Per-group design
#' rank (\code{re_rank}) is a diagnostic check only and does not subset the
#' data.  For \code{family = binomial()}, \code{Prior_Setup_lmebayes()} also
#' fits a classical \code{glm} within each algebraically full-rank group and
#' records \code{design$re_estimable} (whether a finite MLE exists); this does
#' not subset the reference \code{glmer} fit.  Random-effect variances are
#' treated as fixed at their mixed-model estimates.  The returned object
#' two-block Gibbs sampler:
#'
#' \strong{Block 1} (per-group, independent):
#' \deqn{p(\mathbf{b}_j \mid \mathbf{y}, \mathrm{fixef}, \sigma^2, \Sigma_b)
#'       = \mathcal{N}(\boldsymbol{\mu}_{b,j}^*, \boldsymbol{\Sigma}_{b,j}^*)}
#' \deqn{\boldsymbol{\Sigma}_{b,j}^{*-1}
#'       = \mathbf{Z}_j'\mathbf{Z}_j / \sigma^2
#'         + \mathrm{diag}(1/\tau^2_k)}
#' when \code{family = gaussian()}.  For non-Gaussian families there is no
#' observation-level dispersion; Block~1 uses \code{dNormal} with
#' \code{ddef = TRUE} (see \code{\link[glmbayesCore]{dNormal}}).
#'
#' \strong{Block 2} (per-RE coefficient \eqn{k}, independent):
#' \deqn{p(\mathrm{fixef}_k \mid \mathbf{b}_k, \tau^2_k)
#'       = \mathcal{N}(\boldsymbol{\mu}_{\mathrm{fixef},k}^*,
#'                     \boldsymbol{\Sigma}_{\mathrm{fixef},k}^*)}
#' \deqn{\boldsymbol{\Sigma}_{\mathrm{fixef},k}^{*-1}
#'       = \mathbf{X}_k'\mathbf{X}_k / \tau^2_k
#'         + \boldsymbol{\Sigma}_{\mathrm{fixef},k}^{-1}}
#'
#' @param formula Mixed-model formula passed to \code{\link{model_setup}},
#'   whose reference \code{lmer}/\code{glmer} fits (all groups) supply the
#'   calibration quantities.
#' @param data Data frame containing all variables in \code{formula}.
#' @param family Model \code{\link[stats]{family}}.  Default \code{gaussian()}.
#'   Non-Gaussian families use \code{\link[lme4]{glmer}} for calibration;
#'   \code{dispersion_ranef} is omitted (analogous to
#'   \code{\link[glmbayesCore]{Prior_Setup}} for flat GLMs).
#' @param pwt Prior weight(s) in \eqn{(0, 1)}.  Either a \strong{scalar}
#'   (applied to every random-effect component and every Block~2 predictor),
#'   or a \strong{list with one element per random-effect component} (named
#'   with the RE coefficient names in any order, or unnamed positional).
#'   Each list element is a scalar (recycled over that component's Block~2
#'   predictors) or a vector of length \eqn{p_k} (optionally named with the
#'   predictor column names of \code{X_hyper[[k]]}, reordered to match).
#'   The prior covariance for each \code{fixef_k} block is scaled relative to
#'   \code{vcov(fit_ref)} following the \code{\link[glmbayesCore]{Prior_Setup}}
#'   convention: \eqn{(1-\mathrm{pwt})/\mathrm{pwt}} for a scalar, and
#'   elementwise
#'   \eqn{\sqrt{(1-\mathrm{pwt}_i)/\mathrm{pwt}_i}\,
#'        \sqrt{(1-\mathrm{pwt}_j)/\mathrm{pwt}_j}} for vectors.
#' @param pwt_dispersion Optional \emph{relative} prior weight(s) in
#'   \eqn{(0, 1)} for the Block~2 dispersion (precision) prior, decoupled from
#'   \code{pwt}.  A scalar, or a list / numeric vector with one value per
#'   random-effect component (named or positional).  Converted internally to
#'   an effective prior sample size \eqn{n_k = J\,w_k/(1-w_k)} where \eqn{J}
#'   is the number of groups.  At most one of \code{pwt_dispersion} and
#'   \code{n_prior_dispersion} may be supplied; when neither is, the value
#'   is derived from \code{pwt} (the mean across a component's predictors),
#'   keeping the Block~2 \eqn{\tau^2_k} dispersion prior consistent with the
#'   coefficient prior strength.  Weak values carry no computational penalty for
#'   \code{dIndependent_Normal_Gamma} sampling: the \eqn{\tau^2} truncation
#'   window comes from limiting-posterior quantiles independent of the
#'   prior strength (see \code{ing_prior} below).
#' @param pwt_measurement Optional relative prior weight(s) in \eqn{(0, 1)} for
#'   the Block~1 observation \eqn{\sigma^2} Gamma prior (Gaussian models only),
#'   decoupled from \code{pwt} (Block~2 fixef) and from \code{pwt_dispersion}
#'   (Block~2 \eqn{\tau^2_k}).  Either a \strong{scalar} converted to
#'   \eqn{n_{\mathrm{prior}} = w/(1-w)\times n} on the total observation count
#'   \eqn{n} (pooled \code{ing_prior_measurement}), or a \strong{named or
#'   positional vector of length \eqn{J}} (\eqn{J =} \code{nlevels(groups)})
#'   with one weight per group level for per-group \code{dGamma_list()}
#'   calibration (\eqn{n_{\mathrm{prior},j} = w_j/(1-w_j)\times n_j}).  When
#'   \code{pwt_measurement} is a vector, the pooled
#'   \code{ing_prior_measurement} continues to use the default
#'   \code{w = 0.01} on total \eqn{n} (unless \code{n_prior_measurement} is
#'   supplied explicitly).  At most one of \code{pwt_measurement} and
#'   \code{n_prior_measurement} may be supplied; when neither is, scalar
#'   \code{pwt_measurement = 0.01} is used for both pooled and per-group paths.
#' @param n_prior_measurement Optional positive scalar: absolute effective
#'   prior sample size for the Block~1 \eqn{\sigma^2} prior (observation units,
#'   not groups).  See \code{pwt_measurement}.
#' @param dispformula One-sided formula selecting the Block~1 measurement-
#'   dispersion structure: \code{~1} (default, pooled) or \code{~<group_name>}
#'   (matching the random-effects grouping factor exactly, requesting
#'   per-group dispersion).  Gaussian models only compute the per-group
#'   \code{ing_prior_measurement_group} calibration (a within-group
#'   regression fit for every group level, used only by
#'   \code{\link{dGamma_list}()}) when \code{dispformula} is the group
#'   formula; \code{dispformula = ~1} skips it entirely, so \code{pwt_measurement}
#'   / \code{n_prior_measurement} must then be scalar (pooled).  Mirrors the
#'   \code{dispformula} argument on \code{lmerb()}/\code{glmerb()} in
#'   \pkg{lmebayes}, which gates the analogous choice of sampler route; the
#'   two are independent arguments that must be kept consistent by the caller.
#' @param max_disp_perc Scalar in \eqn{(0.5, 1)}, default \code{0.99}.
#'   Tail probability used to compute the \eqn{\sigma^2} (Block~1) and
#'   \eqn{\tau^2_k} (Block~2) truncation windows stored in
#'   \code{ing_prior_measurement} and \code{ing_prior} respectively.
#'   The window is the central \eqn{2 \times \mathrm{max\_disp\_perc} - 1}
#'   mass interval: \code{disp_lower} = \eqn{(1-\mathrm{max\_disp\_perc})}
#'   quantile and \code{disp_upper} = \eqn{\mathrm{max\_disp\_perc}} quantile
#'   of the relevant Gamma (precision) distribution, inverted to the
#'   dispersion scale.  Tighter values (e.g.\ \code{0.95}) shrink the
#'   truncation window and typically improve Block~1 acceptance rates at the
#'   cost of slightly less envelope coverage; looser values (e.g.\ \code{0.999})
#'   widen the window.  Passed to \code{\link{dGamma}()} as \code{max_disp_perc}.
#' @param n_prior_dispersion Optional \emph{absolute} effective prior sample
#'   size(s) (in group units) for the Block~2 \eqn{\tau^2_k} dispersion prior.
#'   A positive scalar, or a list / numeric vector with one value per
#'   random-effect component (named or positional).  See \code{pwt_dispersion}.
#' @param intercept_source Character string controlling the prior mean for the
#'   global intercept hyperparameter \code{(Intercept)::(Intercept)} only.
#'   One of \code{"null_model"} (default) or \code{"full_model"}.  When
#'   \code{"null_model"}, the prior mean is taken from a random-intercept-only
#'   reference fit \code{y ~ 1 + (1 | group)} that omits all fixed-effect
#'   predictors (analogous to \code{\link[glmbayesCore]{Prior_Setup}} with
#'   \code{intercept_source = "null_model"}).  When \code{"full_model"}, the
#'   full-model MLE intercept is used.
#' @param effects_source Character string controlling the prior mean for all
#'   other Block~2 hyperparameters (including population-mean slopes stored as
#'   \code{(Intercept)} columns in non-intercept RE components, and any
#'   non-intercept columns in \code{X_hyper}).  One of \code{"null_effects"}
#'   (default) or \code{"full_model"}.  When \code{"null_effects"}, prior
#'   means are set to zero.  When \code{"full_model"}, full-model MLE values
#'   are used.
#'
#' @return Object of class \code{"lmebayes_prior_setup"} with fields:
#'   \describe{
#'     \item{\code{formula}}{Model formula.}
#'     \item{\code{family}}{Family object.}
#'     \item{\code{pwt}}{Prior weight(s) used: the scalar as supplied, or the
#'       canonical named list of per-predictor weight vectors when a list was
#'       supplied.}
#'     \item{\code{pwt_dispersion}}{Named per-component vector of relative
#'       dispersion prior weights (always present; consistent with
#'       \code{n_prior_dispersion} via \eqn{w_k = n_k/(n_k + J)}).}
#'     \item{\code{n_prior_dispersion}}{Named per-component vector of
#'       effective prior sample sizes for the Block~2 dispersion prior
#'       (always present; used by
#'       \code{\link[=pfamily_list.lmebayes_prior_setup]{pfamily_list}()} to
#'       calibrate \code{dIndependent_Normal_Gamma} components).}
#'     \item{\code{pwt_measurement}}{Gaussian models only: scalar or length-\eqn{J}
#'       vector of relative prior weights for Block~1 \eqn{\sigma^2}.}
#'     \item{\code{n_prior_measurement}}{Gaussian models only: scalar effective
#'       prior sample size for pooled Block~1 \eqn{\sigma^2} on the observation
#'       scale.}
#'     \item{\code{pwt_measurement_group}}{Gaussian models only: named length-\eqn{J}
#'       vector of resolved per-group measurement prior weights (see
#'       \code{dGamma_list()}).}
#'     \item{\code{n_prior_measurement_group}}{Gaussian models only: named
#'       length-\eqn{J} vector of per-group \eqn{n_{\mathrm{prior},j}}.}
#'     \item{\code{block_formula}}{Within-group Block~1 formula: response
#'       regressed on the random-coefficient structure only (columns of
#'       \code{design$re_coef_names}); level-2 hyper covariates are excluded.
#'       Used by \code{dGamma_list()}.}
#'     \item{\code{sd_tau}}{Named vector \code{sqrt(vcov_re)} from the reference
#'       fit; shared population RE standard deviations for per-group calibration.}
#'     \item{\code{data}}{Data frame passed to \code{Prior_Setup_lmebayes()}
#'       (reference for \code{dGamma_list()} diagnostics).}
#'     \item{\code{max_disp_perc}}{The \code{max_disp_perc} value used for both
#'       the \eqn{\sigma^2} and \eqn{\tau^2_k} truncation windows.}
#'     \item{\code{design}}{Full \code{\link{model_setup}} object (all groups).}
#'     \item{\code{mer_fit}}{Reference \code{lmer}/\code{glmer} fit on all
#'       groups (the full-formula fit from \code{\link{model_setup}}), always
#'       present regardless of \code{dispformula} (backs
#'       \code{dispersion_ranef} and \code{x$lmer}/\code{x$glmer} in
#'       \code{lmerb()}/\code{glmerb()}).}
#'     \item{\code{fit_ref}}{The calibration reference for Block~2
#'       (\code{fixef}/\eqn{\tau^2_k}) and the per-group Block~1 inputs
#'       (\code{sd_tau}, BLUP coefficients): identical to \code{mer_fit} when
#'       \code{dispformula = ~1}; otherwise an equivalent \code{glmmTMB}
#'       fit with the same \code{dispformula} (Gaussian models only). See
#'       \code{calibration_source}.}
#'     \item{\code{dispersion_fit}}{\code{NULL} when \code{dispformula = ~1};
#'       otherwise the same \code{glmmTMB} object as \code{fit_ref}.}
#'     \item{\code{sigma2_group}}{\code{NULL} when \code{dispformula = ~1};
#'       otherwise a named length-\eqn{J} vector of per-group observation-level
#'       dispersion read from \code{fit_ref}'s dispersion linear predictor
#'       (\code{glmmTMB::predict(fit_ref, type = "disp")}, aggregated by
#'       group). Diagnostic only -- not the value fed to the sampler; compare
#'       against \code{ing_prior_measurement_group}'s \code{sigma2_hat}.}
#'     \item{\code{calibration_source}}{\code{"lme4"} or \code{"glmmTMB"};
#'       which reference fit produced \code{fit_ref} (and therefore
#'       \code{prior_list}, \code{ing_prior}, \code{sd_tau}, and
#'       \code{ing_prior_measurement_group}).}
#'     \item{\code{dispersion_ranef}}{Scalar \eqn{\sigma^2} for Gaussian models
#'       only; \code{NULL} otherwise. Always derived from \code{mer_fit}
#'       (pooled), independent of \code{dispformula}.}
#'     \item{\code{Sigma_ranef}}{Diagonal RE covariance matrix (Block~1).}
#'     \item{\code{prior_list}}{Named Block~2 prior list per RE coefficient.}
#'     \item{\code{ing_prior}}{Named per-component list of the prospective
#'       \code{dIndependent_Normal_Gamma} calibration: Gamma precision-prior
#'       \code{shape} \eqn{= (n_0 + 1 + p_k)/2} and \code{rate}
#'       \eqn{= \hat\tau^2_k (n_0 + p_k - 1)/2} (the glmbayesCore default
#'       calibration with \eqn{n_0 =} \code{n_prior_dispersion}; since
#'       \code{rate} \eqn{= \hat\tau^2_k (\code{shape} - 1)}, the implied
#'       inverse-Gamma prior on \eqn{\tau^2_k} has mean exactly
#'       \eqn{\hat\tau^2_k}), and the default \eqn{\tau^2_k} truncation
#'       window \code{disp_lower} / \code{disp_upper}: the
#'       \eqn{(1-\mathrm{max\_disp\_perc})} / \eqn{\mathrm{max\_disp\_perc}}
#'       quantiles of the \emph{limiting posterior}
#'       \eqn{\Gamma((J+1)/2,\; \hat\tau^2_k (J-1)/2)} -- the weak-prior
#'       (\eqn{n_0 \to 0}) limit of the Block~2 posterior Gamma for the
#'       precision (glmbayesCore Chapter A12, Theorem 2; inverted to a
#'       \eqn{\tau^2} interval).  This window is identical for all
#'       \eqn{n_0}, covers \eqn{\ge} \eqn{2 \times \mathrm{max\_disp\_perc} - 1}
#'       of the exact posterior for every prior strength, and keeps the
#'       envelope sampler's cost stable as priors weaken; see
#'       \code{inst/ING_TRUNCATION_WINDOW.md}.  Used by
#'       \code{\link[=pfamily_list.lmebayes_prior_setup]{pfamily_list}()} when
#'       \code{ptypes = "dIndependent_Normal_Gamma"}; ignored for
#'       \code{dNormal} priors.}
#'     \item{\code{ing_prior_measurement}}{Gaussian models only: prospective
#'       \code{dGamma()} \code{dispersion_ranef} calibration for Block~1 ING
#'       (observation \eqn{\sigma^2} shared across all group levels):
#'       mean-matched \code{shape} / \code{rate} with
#'       \eqn{n_{\mathrm{prior}} = \mathrm{pwt\_measurement}/(1-\mathrm{pwt\_measurement})\times n},
#'       \eqn{p = p_{\mathrm{re}}}, and \eqn{\hat\sigma^2} =
#'       \code{dispersion_ranef} (same ING algebra as \code{ing_prior} for
#'       \eqn{\tau^2_k}; requires \code{pwt_measurement} \eqn{\le 0.5}), plus
#'       \code{disp_lower} / \code{disp_upper} as the central
#'       \eqn{2 \times \mathrm{max\_disp\_perc} - 1} prior-mass interval from
#'       the same calibrated \code{shape}/\code{rate}.  Pass these fields to
#'       \code{\link{dGamma}()}.}
#'     \item{\code{ing_prior_measurement_group}}{Gaussian models only, and only
#'       when \code{dispformula} requests per-group dispersion: named list
#'       (one entry per group level) of per-group \code{dGamma()} density
#'       calibration (\code{sigma2_hat}, \code{shape_ING}, \code{rate},
#'       \code{rate_gamma}, \code{n_prior}, \code{n_j}, \code{n_combined},
#'       \ldots).  \code{NULL} when \code{dispformula = ~1}.  Used by
#'       \code{\link{dGamma_list}()} for both the \code{dGamma()} density and
#'       its mean-matched truncation bounds.  This is the Chapter A12 Part~VI
#'       "full marginal" calibration (\code{mu_j -> mu_bar}, integrating out
#'       \code{fixef} against Block~2's own \code{dNormal} prior; see
#'       \code{inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md} Part VI), falling
#'       back to the Part~I ("classical") calibration -- with a
#'       \code{\link[base]{warning}} -- if Part VI's calibration fails (e.g.\
#'       a random-effect component's Block~2 hyper-design is not
#'       intercept-only).}
#'     \item{\code{ing_prior_measurement_group_classical}}{Gaussian models
#'       only, and only when \code{dispformula} requests per-group
#'       dispersion: the Part~I calibration (same fields as
#'       \code{ing_prior_measurement_group}), kept unconditionally as the
#'       fallback for -- and for comparison against --
#'       \code{ing_prior_measurement_group} above.  \code{NULL} when
#'       \code{dispformula = ~1}. Not consumed by \code{dGamma_list()}.}
#'     \item{\code{dispformula}}{The \code{dispformula} supplied.}
#'   }
#' @details
#' \strong{Why default calibration depends on classical estimates.}
#' \code{Prior_Setup_lmebayes} scales Block~2 covariances from
#' \code{vcov(fit_ref)} by \eqn{(1-\mathrm{pwt})/\mathrm{pwt}} and plugs in
#' RE variances from the full reference fit, where \code{fit_ref} is the
#' pooled \code{lmer}/\code{glmer} fit when \code{dispformula = ~1}, or an
#' equivalent \code{glmmTMB} fit with the same \code{dispformula} when
#' per-group dispersion is requested (see \code{calibration_source} and
#' \code{mer_fit} above).  By default the global intercept
#' prior mean comes from a random-intercept-only null fit; all other prior
#' means are zero (\code{effects_source = "null_effects"}).  This requires:
#' \enumerate{
#'   \item Converged reference \code{lmer}/\code{glmer} fit from
#'     \code{\link{model_setup}} on the full formula (and a random-intercept-only
#'     null fit when \code{intercept_source = "null_model"}).
#'     Fits with \code{lme4} \code{checkConv} failures (e.g.\ large
#'     \code{max|grad|}) are rejected.
#'   \item Every \code{X_hyper[[k]]} column maps to a \code{fixef(fit_ref)} term.
#'   \item Each RE variance \eqn{\tau^2_k} from the reference fit is strictly positive.
#' }
#' @seealso \code{\link{model_setup}}, \code{\link[glmbayesCore]{Prior_Setup}},
#'   \code{\link[glmbayesCore]{build_mu_all}}
#' @export
Prior_Setup_lmebayes <- function(formula,
                                 data,
                                 family = gaussian(),
                                 pwt    = 0.01,
                                 pwt_dispersion = NULL,
                                 n_prior_dispersion = NULL,
                                 pwt_measurement = NULL,
                                 n_prior_measurement = NULL,
                                 dispformula = ~1,
                                 max_disp_perc = 0.99,
                                 intercept_source = c("null_model", "full_model"),
                                 effects_source   = c("null_effects", "full_model")) {

  intercept_source <- match.arg(intercept_source)
  effects_source   <- match.arg(effects_source)

  if (!is.numeric(max_disp_perc) || length(max_disp_perc) != 1L ||
      is.na(max_disp_perc) || max_disp_perc <= 0.5 || max_disp_perc >= 1) {
    stop(
      "'max_disp_perc' must be a scalar in (0.5, 1).",
      call. = FALSE
    )
  }

  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }
  if (is.numeric(pwt)) {
    if (length(pwt) != 1L || is.na(pwt) || pwt <= 0 || pwt >= 1) {
      stop(
        "'pwt' must be a scalar in (0, 1) or a list with one element per ",
        "random-effect component.",
        call. = FALSE
      )
    }
  } else if (!is.list(pwt)) {
    stop(
      "'pwt' must be a scalar in (0, 1) or a list with one element per ",
      "random-effect component.",
      call. = FALSE
    )
  }
  if (!is.null(pwt_dispersion) && !is.null(n_prior_dispersion)) {
    stop(
      "Supply at most one of 'pwt_dispersion' and 'n_prior_dispersion'.",
      call. = FALSE
    )
  }
  if (!is.null(pwt_measurement) && !is.null(n_prior_measurement)) {
    stop(
      "Supply at most one of 'pwt_measurement' and 'n_prior_measurement'.",
      call. = FALSE
    )
  }
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  if (!inherits(family, "family") || is.null(family$family)) {
    stop("'family' must be a family object.", call. = FALSE)
  }

  is_gaussian <- identical(family$family, "gaussian")
  mer_label   <- if (is_gaussian) "lmer" else "glmer"

  if (is_gaussian) {
    ctrl <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  } else {
    # Explicit glmerControl() can crash some lme4 builds; let glmer defaults apply.
    ctrl <- NULL
  }

  design <- model_setup(
    formula = formula,
    data = data,
    family = family,
    control = ctrl
  )

  dispformula_kind <- .lmebayes_prior_setup_dispformula_kind(
    dispformula, design$group_name
  )
  if (identical(dispformula_kind, "group") && !is_gaussian) {
    stop(
      "'dispformula' must be ~1 for family = ", family$family,
      "() (no observation-level dispersion parameter for per-group ",
      "measurement-dispersion calibration).",
      call. = FALSE
    )
  }
  if (identical(dispformula_kind, "pooled") && is_gaussian &&
      (length(pwt_measurement) > 1L || length(n_prior_measurement) > 1L)) {
    stop(
      "'pwt_measurement'/'n_prior_measurement' has more than one value ",
      "(per-group), but dispformula = ~1 (pooled); use ",
      "dispformula = ~", design$group_name, " to calibrate per-group ",
      "measurement-dispersion priors for dGamma_list().",
      call. = FALSE
    )
  }

  glm_est <- .lmebayes_block_glm_estimable(
    y       = design$y,
    groups  = design$groups,
    Z       = design$Z,
    re_rank = design$re_rank,
    family  = family
  )
  design$re_estimable <- glm_est$re_estimable
  design$re_glm_check <- glm_est$re_glm_check

  ## Full-rank status is a per-group DESIGN CHECK only (reported by print();
  ## groups with rank-deficient Z_j are still fully used below).  The lme4
  ## reference fit is always fit by model_setup() (backward-compat: it backs
  ## x$lmer/x$glmer in lmerb()/glmerb() and the pooled dispersion_ranef
  ## scalar, regardless of dispformula).
  mer_fit <- if (is_gaussian) design$lmer_fit else design$glmer_fit
  if (is.null(mer_fit)) {
    stop(
      "model_setup() did not return a reference ", mer_label, " fit.",
      call. = FALSE
    )
  }

  mer_issues <- .lmebayes_mer_convergence_issues(
    mer_fit, sprintf("%s (full formula)", mer_label)
  )
  if (length(mer_issues) > 0L) {
    stop(
      "Prior_Setup_lmebayes() requires converged ", mer_label,
      " reference fits:\n  - ",
      paste(mer_issues, collapse = "\n  - "),
      "\n\nRevise the model or supply hyperpriors manually without ",
      "Prior_Setup_lmebayes().",
      call. = FALSE
    )
  }

  dispersion_ranef <- if (is_gaussian) design$residual_var else NULL

  ## Calibration reference: when dispformula requests per-group dispersion,
  ## Block~2 (fixef/tau^2_k) and the per-group Block~1 inputs (sd_tau, BLUP
  ## coefficients) come from an equivalent glmmTMB fit instead of the pooled
  ## lme4 fit above, so the heteroscedastic structure that dispformula
  ## requests is reflected in every calibrated quantity, not only in
  ## dGamma_list()'s per-group densities. See
  ## inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md.
  calibration_source <- "lme4"
  sigma2_group_ref   <- NULL
  if (is_gaussian && identical(dispformula_kind, "group")) {
    fit_ref <- .lmebayes_fit_glmmtmb_reference(
      formula     = formula,
      data        = data,
      family      = family,
      dispformula = dispformula,
      REML        = TRUE
    )
    tmb_issues <- .lmebayes_glmmtmb_convergence_issues(
      fit_ref, "glmmTMB (full formula, per-group dispersion)"
    )
    if (length(tmb_issues) > 0L) {
      stop(
        "Prior_Setup_lmebayes() requires a converged glmmTMB reference fit ",
        "for dispformula = ~", design$group_name, ":\n  - ",
        paste(tmb_issues, collapse = "\n  - "),
        "\n\nRevise the model or supply hyperpriors manually without ",
        "Prior_Setup_lmebayes().",
        call. = FALSE
      )
    }
    vc_ref   <- extract_glmmtmb_variance_components(
      fit_ref, design$re_coef_names, design$group_name
    )
    tau2_vec <- vc_ref$vcov_re
    sigma2_group_ref <- .lmebayes_glmmtmb_group_sigma2(
      fit_ref, design$group_name, levels(design$groups)
    )
    calibration_source <- "glmmTMB"
  } else {
    fit_ref  <- mer_fit
    tau2_vec <- design$vcov_re
  }
  ref_label <- if (identical(calibration_source, "glmmTMB")) "glmmTMB" else mer_label

  p_re        <- length(design$re_coef_names)
  Sigma_ranef <- diag(unname(tau2_vec), nrow = p_re, ncol = p_re)
  dimnames(Sigma_ranef) <- list(design$re_coef_names, design$re_coef_names)

  fe   <- .lmebayes_reference_fixef(fit_ref)
  V_fe <- .lmebayes_reference_vcov(fit_ref)

  fe_name_for <- function(k, col) {
    if (k == "(Intercept)") {
      if (col %in% names(fe)) col else NA_character_
    } else if (col == "(Intercept)") {
      if (k %in% names(fe)) k else NA_character_
    } else {
      cand <- c(paste0(col, ":", k), paste0(k, ":", col))
      hit  <- cand[cand %in% names(fe)]
      if (length(hit)) hit[1L] else NA_character_
    }
  }

  re_names  <- design$re_coef_names
  tau_tol   <- sqrt(.Machine$double.eps)
  re_issues <- character(0)

  for (k in re_names) {
    X_k    <- design$X_hyper[[k]]
    cols_k <- colnames(X_k)
    fe_nms <- vapply(cols_k, fe_name_for, character(1L), k = k)
    miss_idx <- is.na(fe_nms) | !fe_nms %in% names(fe)

    if (any(miss_idx)) {
      if (k != "(Intercept)" &&
          length(cols_k) == 1L &&
          identical(cols_k, "(Intercept)")) {
        re_issues <- c(
          re_issues,
          sprintf(
            paste0(
              "%s: random slope has no fixed main effect in ", ref_label,
              " (add '%s' to the fixed part of the formula)"
            ),
            k, k
          )
        )
      } else {
        expected_fe <- vapply(seq_along(cols_k), function(i) {
          col <- cols_k[i]
          if (k == "(Intercept)") {
            col
          } else if (col == "(Intercept)") {
            k
          } else {
            paste0(col, ":", k)
          }
        }, character(1L))
        re_issues <- c(
          re_issues,
          sprintf(
            "%s: no %s fixed effect for %s",
            k, ref_label,
            paste(expected_fe[miss_idx], collapse = ", ")
          )
        )
      }
    }

    tau2_k <- unname(tau2_vec[[k]])
    if (is.na(tau2_k) || tau2_k <= tau_tol) {
      re_issues <- c(
        re_issues,
        sprintf(
          paste0(
            "%s: random-effect variance is zero or on the boundary ",
            "(singular fit); group-level variation is not identified"
          ),
          k
        )
      )
    }
  }

  if (length(re_issues) > 0L) {
    stop(
      "Prior_Setup_lmebayes() cannot calibrate default hyperpriors:\n  - ",
      paste(re_issues, collapse = "\n  - "),
      "\n\nRevise the formula (e.g. add a fixed main effect for each random ",
      "slope and avoid RE terms with zero estimated variance), or supply ",
      "hyperpriors manually without Prior_Setup_lmebayes().",
      call. = FALSE
    )
  }

  pwt_list <- .lmebayes_resolve_pwt(pwt, design)

  J_groups <- nlevels(design$groups)
  disp     <- .lmebayes_resolve_disp_prior(
    pwt_dispersion     = pwt_dispersion,
    n_prior_dispersion = n_prior_dispersion,
    J                  = J_groups,
    re_names           = re_names,
    pwt_list           = pwt_list
  )

  # ---- null intercept model for global intercept prior mean ----------------
  fe_null <- fe
  if (intercept_source == "null_model") {
    resp_nm  <- all.vars(formula)[1L]
    grp_nm   <- design$group_name
    null_formula <- stats::as.formula(
      paste(resp_nm, "~ 1 + (1 |", grp_nm, ")", sep = "")
    )
    if (identical(calibration_source, "glmmTMB")) {
      ## Match the heteroscedastic reference so the global intercept prior
      ## mean is calibrated consistently with fit_ref above.
      null_fit <- .lmebayes_fit_glmmtmb_reference(
        formula     = null_formula,
        data        = data,
        family      = family,
        dispformula = dispformula,
        REML        = TRUE
      )
      null_issues <- .lmebayes_glmmtmb_convergence_issues(
        null_fit, "glmmTMB (null intercept model, per-group dispersion)"
      )
    } else {
      null_fit <- if (is_gaussian) {
        lme4::lmer(null_formula, data = data, REML = TRUE, control = ctrl)
      } else if (is.null(ctrl)) {
        lme4::glmer(null_formula, family = family, data = data)
      } else {
        lme4::glmer(null_formula, family = family, data = data, control = ctrl)
      }
      null_issues <- .lmebayes_mer_convergence_issues(
        null_fit, sprintf("%s (null intercept model)", mer_label)
      )
    }
    if (length(null_issues) > 0L) {
      stop(
        "Prior_Setup_lmebayes() requires a converged ", ref_label,
        " random-intercept-only null fit for intercept_source = \"null_model\":\n  - ",
        paste(null_issues, collapse = "\n  - "),
        "\n\nUse intercept_source = \"full_model\" or revise the model.",
        call. = FALSE
      )
    }
    fe_null <- .lmebayes_reference_fixef(null_fit)
  }

  prior_list <- stats::setNames(
    lapply(re_names, function(k) {
      X_k    <- design$X_hyper[[k]]
      cols_k <- colnames(X_k)
      p_k    <- length(cols_k)
      tau2_k <- tau2_vec[[k]]

      fe_nms <- vapply(cols_k, fe_name_for, character(1L), k = k)
      fe_idx <- fe_nms

      mu_fixef <- vapply(seq_len(p_k), function(i) {
        col <- cols_k[i]
        if (identical(k, "(Intercept)") && identical(col, "(Intercept)")) {
          if (intercept_source == "null_model") {
            if (!("(Intercept)" %in% names(fe_null))) {
              stop(
                "Null intercept model did not return an (Intercept) fixed effect.",
                call. = FALSE
              )
            }
            unname(fe_null["(Intercept)"])
          } else {
            unname(fe[fe_nms[i]])
          }
        } else if (effects_source == "null_effects") {
          0
        } else {
          unname(fe[fe_nms[i]])
        }
      }, numeric(1L))
      names(mu_fixef) <- cols_k

      ## Elementwise scaling sqrt(s_i) * sqrt(s_j) with s_i = (1-w_i)/w_i;
      ## reduces to the scalar (1-pwt)/pwt factor when all weights are equal.
      sc_k <- sqrt((1 - pwt_list[[k]]) / pwt_list[[k]])
      Sigma_fixef <- V_fe[fe_idx, fe_idx, drop = FALSE] * outer(sc_k, sc_k)
      dimnames(Sigma_fixef) <- list(cols_k, cols_k)

      list(
        mu_fixef         = mu_fixef,
        Sigma_fixef      = Sigma_fixef,
        dispersion_fixef = tau2_k
      )
    }),
    re_names
  )

  pwt_out <- if (is.numeric(pwt)) pwt else pwt_list

  ## Prospective dIndependent_Normal_Gamma calibration per component (used
  ## only when pfamily_list(ptypes = "dIndependent_Normal_Gamma") is chosen):
  ## Gamma precision prior shape/rate from the glmbayesCore default
  ## calibration (compute_gaussian_prior() with k = 1):
  ##   shape_ING = (n0 + 1 + p_k)/2,  b_0 = tau2_k * (n0 + p_k - 1)/2.
  ## Since b_0 = tau2_k * (shape_ING - 1), the implied inverse-Gamma prior on
  ## tau^2_k has mean exactly tau2_k for every n0 and p_k.
  ##
  ## The tau^2 truncation window (disp_lower / disp_upper) uses the
  ## *limiting posterior* of glmbayesCore Chapter A12, Theorem 2 -- the
  ## weak-prior (n0 -> 0) limit of the Block 2 posterior Gamma:
  ##   a_inf = (J + 1)/2,  b_inf = tau2_k * (J - 1)/2
  ## (so b_inf/(a_inf - 1) = tau2_k: mean-matched, like the prior).  The
  ## window uses max_disp_perc (default 0.99) for the central
  ## (2*max_disp_perc - 1) mass interval.  Quantiles of the *prior* would
  ## stretch without bound as n0 -> 0 (posterior coverage -> 100%, envelope
  ## acceptance -> 0); the limiting-posterior window instead has coverage >=
  ## (2*max_disp_perc - 1) of the exact posterior for every n0 (the finite-n0
  ## posterior is strictly more concentrated than the limit), is identical for
  ## all n0, and keeps the envelope sampler's candidates-per-draw roughly
  ## constant as priors weaken.  See inst/ING_TRUNCATION_WINDOW.md.
  ## Stored here so print() can display the window and pfamily_list()
  ## consumes one source of truth.
  ing_prior <- stats::setNames(
    lapply(re_names, function(k) {
      n0_k    <- unname(disp$n_prior_dispersion[[k]])
      p_k     <- length(prior_list[[k]]$mu_fixef)
      tau2_k  <- unname(prior_list[[k]]$dispersion_fixef)
      shape_k <- (n0_k + 1) / 2 + p_k / 2
      rate_k  <- tau2_k * (n0_k + p_k - 1) / 2
      win_k <- .lmebayes_ing_limiting_posterior_window(tau2_k, J_groups,
                                                       max_disp_perc)
      list(
        shape         = shape_k,
        rate          = rate_k,
        disp_lower    = win_k$disp_lower,
        disp_upper    = win_k$disp_upper,
        max_disp_perc = max_disp_perc
      )
    }),
    re_names
  )

  ing_prior_measurement <- NULL
  meas <- NULL
  if (is_gaussian) {
    n_obs <- length(design$y)
    pwt_meas_vector <- !is.null(pwt_measurement) && length(pwt_measurement) > 1L
    if (pwt_meas_vector && !is.numeric(pwt_measurement)) {
      stop(
        "'pwt_measurement' vector must be numeric.",
        call. = FALSE
      )
    }
    meas <- .lmebayes_resolve_measurement_disp_prior(
      pwt_measurement     = if (pwt_meas_vector) NULL else pwt_measurement,
      n_prior_measurement = n_prior_measurement,
      n_obs               = n_obs
    )
    ing_prior_measurement <- .lmebayes_calibrate_ing_prior_measurement(
      design           = design,
      dispersion_ranef = dispersion_ranef,
      n_prior          = meas$n_prior_measurement,
      max_disp_perc    = max_disp_perc
    )
  }

  group_levels <- levels(design$groups)
  block_formula <- .lmebayes_block_formula_from_re(formula, re_names)
  sd_tau_out <- if (is_gaussian) {
    stats::setNames(sqrt(unname(tau2_vec)), re_names)
  } else {
    NULL
  }

  meas_group <- if (is_gaussian) {
    n_j <- as.integer(table(design$groups))
    names(n_j) <- group_levels
    .lmebayes_resolve_measurement_disp_prior_group(
      pwt_measurement     = pwt_measurement,
      n_prior_measurement = NULL,
      n_j                 = n_j,
      group_levels        = group_levels
    )
  } else {
    NULL
  }

  ## Part I ("classical"): per-group mu_j / Sigma_j from the group's own
  ## within-group null-model regression (A12 3.3.4 rate = S_marg). Always
  ## computed -- both as the fallback below, and as the reference column in
  ## the comparison prints. See inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md Part I.
  ing_prior_measurement_group_classical <- if (is_gaussian && identical(dispformula_kind, "group")) {
    .lmebayes_calibrate_ing_prior_measurement_group(
      design           = design,
      data             = data,
      block_formula    = block_formula,
      sd_tau           = sd_tau_out,
      pwt_group        = meas_group$pwt_measurement,
      n_prior_group    = meas_group$n_prior_measurement,
      group_levels     = group_levels,
      family           = family,
      intercept_source = intercept_source,
      effects_source   = effects_source
    )
  } else {
    NULL
  }

  ## Part VI ("full marginal"): as Part I, but also integrates out fixef
  ## against Block~2's own fixef dNormal prior (mu_j -> mu_bar, Sigma_j ->
  ## Sigma_j + Omega; both bar_mu/Omega read off prior_list, population-level
  ## and shared across groups). This is now the production calibration fed to
  ## dGamma_list() (density *and* truncation bounds, since both derive from
  ## sigma2_hat/rate on this same list) -- see
  ## inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md Part VI. Falls back to Part I
  ## ("classical") with a warning if it fails (e.g. a RE component's Block~2
  ## hyper-design isn't intercept-only, which Part VI doesn't yet support).
  ing_prior_measurement_group_full_marginal <- if (
    is_gaussian && identical(dispformula_kind, "group") &&
      !is.null(ing_prior_measurement_group_classical)
  ) {
    tryCatch(
      .lmebayes_calibrate_ing_prior_measurement_group_full_marginal(
        design           = design,
        data             = data,
        block_formula    = block_formula,
        sd_tau           = sd_tau_out,
        pwt_group        = meas_group$pwt_measurement,
        n_prior_group    = meas_group$n_prior_measurement,
        group_levels     = group_levels,
        prior_list       = prior_list,
        family           = family,
        intercept_source = intercept_source,
        effects_source   = effects_source
      ),
      error = function(e) {
        warning(
          "Part VI full-marginal per-group dispersion calibration failed; ",
          "dGamma_list() will fall back to the Part I (classical) ",
          "calibration for this fit: ", conditionMessage(e),
          call. = FALSE
        )
        NULL
      }
    )
  } else {
    NULL
  }

  ## SWAP POINT: production ing_prior_measurement_group (consumed by
  ## dGamma_list() for both the dGamma() density and its mean-matched
  ## truncation bounds). To revert to the Part I ("classical") behavior,
  ## change the right-hand side to
  ## `ing_prior_measurement_group_classical`.
  ing_prior_measurement_group <- if (!is.null(ing_prior_measurement_group_full_marginal)) {
    ing_prior_measurement_group_full_marginal
  } else {
    ing_prior_measurement_group_classical
  }

  ## Dev diagnostic: print A12 3.3.4 rate vs A12 3.3.5 rate_gamma (both from
  ## the Part I calibration; unaffected by the Part VI swap above) and, when
  ## available, the Part I vs Part VI comparison actually driving the swap.
  if (is_gaussian && identical(dispformula_kind, "group") &&
      !is.null(ing_prior_measurement_group_classical)) {
    .lmebayes_print_ing_prior_measurement_group_compare(
      existing = ing_prior_measurement_group_classical
    )

    if (!is.null(ing_prior_measurement_group_full_marginal)) {
      .lmebayes_print_ing_prior_measurement_group_full_marginal_compare(
        existing      = ing_prior_measurement_group_classical,
        full_marginal = ing_prior_measurement_group_full_marginal
      )
    }
  }

  pwt_measurement_out <- if (is_gaussian) {
    if (!is.null(pwt_measurement) && length(pwt_measurement) > 1L) {
      w <- meas_group$pwt_measurement
      attr(w, "source") <- meas_group$source
      w
    } else {
      w <- meas$pwt_measurement
      attr(w, "source") <- meas$source
      w
    }
  } else {
    NULL
  }
  n_prior_measurement_out <- if (is_gaussian) {
    np <- meas$n_prior_measurement
    attr(np, "source") <- meas$source
    np
  } else {
    NULL
  }

  structure(
    list(
      formula            = formula,
      family             = family,
      pwt                = pwt_out,
      pwt_dispersion     = disp$pwt_dispersion,
      n_prior_dispersion = disp$n_prior_dispersion,
      pwt_measurement    = pwt_measurement_out,
      n_prior_measurement = n_prior_measurement_out,
      pwt_measurement_group = if (is_gaussian) meas_group$pwt_measurement else NULL,
      n_prior_measurement_group = if (is_gaussian) meas_group$n_prior_measurement else NULL,
      max_disp_perc      = max_disp_perc,
      dispformula        = dispformula,
      intercept_source   = intercept_source,
      effects_source     = effects_source,
      data               = data,
      block_formula      = block_formula,
      sd_tau             = sd_tau_out,
      design             = design,
      mer_fit            = mer_fit,
      fit_ref            = fit_ref,
      dispersion_fit     = if (identical(calibration_source, "glmmTMB")) fit_ref else NULL,
      sigma2_group       = sigma2_group_ref,
      calibration_source = calibration_source,
      dispersion_ranef   = dispersion_ranef,
      Sigma_ranef        = Sigma_ranef,
      prior_list            = prior_list,
      ing_prior             = ing_prior,
      ing_prior_measurement = ing_prior_measurement,
      ing_prior_measurement_group = ing_prior_measurement_group,
      ing_prior_measurement_group_classical = ing_prior_measurement_group_classical
    ),
    class = "lmebayes_prior_setup"
  )
}
## Classify 'dispformula' as "pooled" (~1) or "group" (~<group_name>, matching
## the random-effects grouping factor exactly).  Gates the (fragile, per-group
## within-group regression) .lmebayes_calibrate_ing_prior_measurement_group()
## call: it now runs only for "group", instead of unconditionally for every
## Gaussian model.  Mirrors lmebayes:::.lmebayes_validate_dispformula(), which
## performs the analogous check against the resolved dispersion_ranef mode at
## lmerb()/glmerb() time; the two dispformula arguments are independent and
## must be kept consistent by the caller.
#' @keywords internal
#' @noRd
.lmebayes_prior_setup_dispformula_kind <- function(dispformula, group_name) {
  if (!inherits(dispformula, "formula") || length(dispformula) != 2L) {
    stop(
      "'dispformula' must be a one-sided formula, either ~1 (pooled) or ~",
      group_name, " (per-group).",
      call. = FALSE
    )
  }
  vars <- all.vars(dispformula)
  if (length(vars) == 0L) {
    return("pooled")
  }
  if (length(vars) == 1L && identical(vars, group_name)) {
    return("group")
  }
  stop(
    "'dispformula' must be ~1 (pooled) or ~", group_name,
    " (per-group, matching the random-effects grouping factor); got ",
    deparse(dispformula), ".",
    call. = FALSE
  )
}

## Resolve 'pwt' (scalar or list) into a canonical named list with one named
## numeric vector per random-effect component, ordered like X_hyper[[k]].
#' @keywords internal
#' @noRd
.lmebayes_resolve_pwt <- function(pwt, design) {

  re_names <- design$re_coef_names
  p_re     <- length(re_names)

  check_range <- function(v, what) {
    if (!is.numeric(v) || anyNA(v) || any(v <= 0) || any(v >= 1)) {
      stop(sprintf("%s must be numeric with all values in (0, 1).", what),
           call. = FALSE)
    }
  }

  if (is.numeric(pwt)) {
    check_range(pwt, "'pwt'")
    return(stats::setNames(
      lapply(re_names, function(k) {
        cols_k <- colnames(design$X_hyper[[k]])
        stats::setNames(rep(as.numeric(pwt), length(cols_k)), cols_k)
      }),
      re_names
    ))
  }

  if (length(pwt) != p_re) {
    stop(sprintf(
      paste0("'pwt' list has length %d but there are %d random-effect ",
             "components (%s)."),
      length(pwt), p_re, paste(re_names, collapse = ", ")
    ), call. = FALSE)
  }

  nms <- names(pwt)
  if (!is.null(nms) && any(nzchar(nms))) {
    if (!setequal(nms, re_names)) {
      stop(
        "Names of 'pwt' must match the random-effect coefficient names: ",
        paste(re_names, collapse = ", "),
        call. = FALSE
      )
    }
    pwt <- pwt[re_names]
  } else {
    names(pwt) <- re_names
  }

  out <- stats::setNames(vector("list", p_re), re_names)
  for (k in re_names) {
    cols_k <- colnames(design$X_hyper[[k]])
    p_k    <- length(cols_k)
    v      <- pwt[[k]]
    what   <- sprintf("'pwt[[\"%s\"]]'", k)
    check_range(v, what)

    if (length(v) == 1L) {
      v <- rep(as.numeric(v), p_k)
    } else if (length(v) == p_k) {
      vn <- names(v)
      if (!is.null(vn) && any(nzchar(vn))) {
        if (!setequal(vn, cols_k)) {
          stop(sprintf(
            "Names of %s must match the Block 2 predictors: %s.",
            what, paste(cols_k, collapse = ", ")
          ), call. = FALSE)
        }
        v <- v[cols_k]
      }
      v <- as.numeric(v)
    } else {
      stop(sprintf(
        "%s must have length 1 or %d (one value per Block 2 predictor).",
        what, p_k
      ), call. = FALSE)
    }
    out[[k]] <- stats::setNames(v, cols_k)
  }
  out
}

## Resolve the Block 2 dispersion-prior weights into mutually consistent
## per-component vectors: n_k = J * w_k / (1 - w_k)  <=>  w_k = n_k / (n_k + J).
#' @keywords internal
#' @noRd
.lmebayes_resolve_disp_prior <- function(pwt_dispersion,
                                         n_prior_dispersion,
                                         J,
                                         re_names,
                                         pwt_list) {

  p_re <- length(re_names)

  expand <- function(x, what) {
    if (is.list(x)) {
      ok <- vapply(
        x, function(e) is.numeric(e) && length(e) == 1L && !is.na(e),
        logical(1L)
      )
      if (!all(ok)) {
        stop(sprintf(
          "'%s' list elements must each be a single numeric value.", what
        ), call. = FALSE)
      }
      nms <- names(x)
      x <- vapply(x, as.numeric, numeric(1L))
      names(x) <- nms
    }
    if (!is.numeric(x) || anyNA(x)) {
      stop(sprintf("'%s' must be numeric without missing values.", what),
           call. = FALSE)
    }
    if (length(x) == 1L) {
      x <- rep(unname(x), p_re)
    } else if (length(x) == p_re) {
      nms <- names(x)
      if (!is.null(nms) && any(nzchar(nms))) {
        if (!setequal(nms, re_names)) {
          stop(sprintf(
            "Names of '%s' must match the random-effect coefficient names: %s.",
            what, paste(re_names, collapse = ", ")
          ), call. = FALSE)
        }
        x <- x[re_names]
      }
    } else {
      stop(sprintf(
        paste0("'%s' must have length 1 or %d (one value per random-effect ",
               "component)."),
        what, p_re
      ), call. = FALSE)
    }
    stats::setNames(as.numeric(x), re_names)
  }

  if (!is.null(pwt_dispersion)) {
    w <- expand(pwt_dispersion, "pwt_dispersion")
    if (any(w <= 0) || any(w >= 1)) {
      stop("'pwt_dispersion' values must be in (0, 1).", call. = FALSE)
    }
    n   <- J * w / (1 - w)
    src <- "user-supplied (pwt_dispersion)"
  } else if (!is.null(n_prior_dispersion)) {
    n <- expand(n_prior_dispersion, "n_prior_dispersion")
    if (any(n <= 0) || any(!is.finite(n))) {
      stop("'n_prior_dispersion' values must be positive and finite.",
           call. = FALSE)
    }
    w   <- n / (n + J)
    src <- "user-supplied (n_prior_dispersion)"
  } else {
    ## Default: derive from the coefficient pwt (mean across predictors per
    ## component), keeping the dispersion prior consistent with the
    ## coefficient prior strength.  This was briefly replaced by a fixed
    ## 0.2 when the ING tau^2 truncation window came from *prior*
    ## quantiles (weak priors widened the window and collapsed envelope
    ## acceptance); the window now uses limiting-posterior quantiles
    ## independent of n0 (inst/ING_TRUNCATION_WINDOW.md), so weak
    ## dispersion priors no longer carry a computational penalty.
    w <- vapply(re_names, function(k) mean(pwt_list[[k]]), numeric(1L))
    n   <- J * w / (1 - w)
    src <- "derived from pwt"
  }

  w <- stats::setNames(w, re_names)
  n <- stats::setNames(n, re_names)
  attr(w, "source") <- src
  attr(n, "source") <- src
  list(pwt_dispersion = w, n_prior_dispersion = n, source = src)
}

#' @rdname Prior_Setup_lmebayes
#' @method print lmebayes_prior_setup
#' @param x Object of class \code{"lmebayes_prior_setup"}.
#' @param digits Number of decimal places for numeric output. Default 4.
#' @param ... Ignored.
#' @return \code{x} invisibly.
#' @export
print.lmebayes_prior_setup <- function(x, digits = 4L, ...) {

  re_names <- x$design$re_coef_names
  n_fr     <- sum(x$design$re_rank)
  n_all    <- nlevels(x$design$groups)

  disp_src <- attr(x$pwt_dispersion, "source")

  cat("Call: Prior_Setup_lmebayes()\n\n")
  cat(sprintf("  family           : %s (%s link)\n",
              x$family$family, x$family$link))
  cat(sprintf("  intercept_source : %s\n",
              if (!is.null(x$intercept_source)) x$intercept_source else "full_model"))
  cat(sprintf("  effects_source   : %s\n",
              if (!is.null(x$effects_source)) x$effects_source else "full_model"))
  cat(sprintf("  dispformula      : %s\n",
              if (!is.null(x$dispformula)) deparse(x$dispformula) else "~1"))
  cat(sprintf("  calibration_source: %s  (fixef/tau2_k/sd_tau reference fit)\n",
              if (!is.null(x$calibration_source)) x$calibration_source else "lme4"))
  if (!is.null(x$dispersion_ranef)) {
    cat(sprintf(
      "  dispersion_ranef : %.4f  (sigma2, fixed from all %d %s)\n",
      x$dispersion_ranef, n_all, x$design$group_name
    ))
    ing_m <- x$ing_prior_measurement
    if (!is.null(ing_m)) {
      mdp <- if (!is.null(ing_m$max_disp_perc)) ing_m$max_disp_perc else 0.99
      cat(sprintf(
        paste0(
          "  ING sigma^2 window: [%.4g, %.4g]  ",
          "(%.4g/%.4g quantiles from calibrated shape/rate; ",
          "upper/sigma2 = %.3g)\n"
        ),
        ing_m$disp_lower, ing_m$disp_upper,
        1 - mdp, mdp,
        ing_m$disp_upper / unname(ing_m$sigma2_hat)
      ))
      cat(sprintf(
        paste0(
          "  ING sigma^2 shape, rate : %.4g, %.4g  ",
          "(mean-matched IG; E[sigma^2] = rate/(shape-1) = %.4g)\n"
        ),
        ing_m$shape, ing_m$rate,
        ing_m$rate / (ing_m$shape - 1)
      ))
      cat(sprintf(
        "  ING sigma^2 n_prior   : %.4g  (= n * pwt_measurement / (1 - pwt_measurement); p_re = %d)\n",
        ing_m$n_prior, ing_m$p_re
      ))
      if (!is.null(x$pwt_measurement)) {
        meas_src <- attr(x$pwt_measurement, "source")
        pwt_disp <- if (length(x$pwt_measurement) == 1L) {
          sprintf("%.4g", x$pwt_measurement)
        } else {
          paste(
            sprintf("%s=%.4g", names(x$pwt_measurement), x$pwt_measurement),
            collapse = ", "
          )
        }
        cat(sprintf(
          "  pwt_measurement : %s  [%s]\n",
          pwt_disp,
          if (is.null(meas_src)) "unknown" else meas_src
        ))
      }
    }
    ing_grp <- x$ing_prior_measurement_group
    if (!is.null(ing_grp)) {
      guard_df <- data.frame(
        group     = names(ing_grp),
        n_j       = vapply(ing_grp, `[[`, 0, "n_j"),
        n_prior   = vapply(ing_grp, `[[`, 0, "n_prior"),
        sigma2_hat = vapply(ing_grp, `[[`, 0, "sigma2_hat"),
        pwt_group = vapply(ing_grp, `[[`, 0, "pwt_group"),
        stringsAsFactors = FALSE
      )
      cat("\n--- Per-group Block~1 sigma^2 calibration (dGamma_list) ---\n")
      print(round(guard_df, digits))
    }
  } else {
    cat("  dispersion_ranef : NULL  (no observation-level dispersion)\n")
  }
  cat(sprintf(
    "  Full-rank groups (algebraic Z_j): %d of %d %s  (design check only)\n",
    n_fr, n_all, x$design$group_name
  ))
  if (identical(x$family$family, "binomial") &&
      !is.null(x$design$re_estimable)) {
    n_est <- sum(x$design$re_estimable[x$design$re_rank])
    cat(sprintf(
      paste0(
        "  Full-rank with glm MLE          : %d of %d full-rank ",
        "(%d of %d total %s)\n"
      ),
      n_est, n_fr, n_est, n_all, x$design$group_name
    ))
  }
  cat("\n")

  cat("--- Sigma_ranef (diagonal RE covariance) ---\n")
  print(round(x$Sigma_ranef, digits))

  cat("\n--- prior_list: mu_fixef / Sigma_fixef / dispersion_fixef (Block 2) ---\n")
  for (nm in re_names) {
    pl <- x$prior_list[[nm]]
    cat(sprintf("\n  [%s]\n", nm))
    pwt_k <- if (is.numeric(x$pwt)) x$pwt else x$pwt[[nm]]
    pwt_str <- if (length(unique(pwt_k)) == 1L) {
      sprintf("%.4g", pwt_k[1L])
    } else {
      paste(sprintf("%s=%.4g", names(pwt_k), pwt_k), collapse = ", ")
    }
    cat(sprintf("  pwt             : %s\n", pwt_str))
    if (!is.null(x$pwt_dispersion)) {
      cat(sprintf(
        "  pwt_disp        : %.4g  [%s]\n",
        x$pwt_dispersion[[nm]],
        if (is.null(disp_src)) "unknown" else disp_src
      ))
    }
    if (!is.null(x$n_prior_dispersion)) {
      cat(sprintf(
        "  n_prior_disp    : %.4g  (= J * pwt_disp / (1 - pwt_disp))\n",
        x$n_prior_dispersion[[nm]]
      ))
    }
    cat("  mu_fixef:\n")
    print(round(pl$mu_fixef, digits))
    cat("  Sigma_fixef:\n")
    print(round(pl$Sigma_fixef, digits))
    cat(sprintf(
      "  dispersion_fixef: %.4f  (RE variance tau^2_k; Block 2 scale)\n",
      pl$dispersion_fixef))
    ing_k <- x$ing_prior[[nm]]
    if (!is.null(ing_k)) {
      mdp_k <- if (!is.null(ing_k$max_disp_perc)) ing_k$max_disp_perc else 0.99
      cat(sprintf(
        "  ING tau^2 window: [%.4g, %.4g]  (%.4g/%.4g limiting-posterior quantiles; upper/tau2 = %.3g)\n",
        ing_k$disp_lower, ing_k$disp_upper,
        1 - mdp_k, mdp_k,
        ing_k$disp_upper / unname(pl$dispersion_fixef)
      ))
      cat(sprintf(
        "  ING shape, rate : %.4g, %.4g  (Gamma prior on 1/tau^2_k; used only with ptypes = \"dIndependent_Normal_Gamma\")\n",
        ing_k$shape, ing_k$rate
      ))
    }
  }

  invisible(x)
}

#' Per-group classical glm MLE existence for binomial Block-1 design
#'
#' For each group with algebraically full-rank \code{Z_j}, fits
#' \code{glm(y ~ Z_j - 1, family = binomial)} and marks the group estimable
#' when all coefficients and \code{vcov} entries are finite.  Non-binomial
#' families return \code{re_estimable = re_rank} (no glm check).
#' @noRd
.lmebayes_block_glm_estimable <- function(y, groups, Z, re_rank, family) {
  g_levs <- names(re_rank)
  if (is.null(g_levs)) {
    g_levs <- levels(groups)
    names(re_rank) <- g_levs
  }

  re_estimable <- stats::setNames(rep(FALSE, length(g_levs)), g_levs)

  if (!identical(family$family, "binomial")) {
    re_estimable[re_rank] <- TRUE
    return(list(
      re_estimable = re_estimable,
      re_glm_check = NULL
    ))
  }

  y <- as.numeric(y)
  g_chr <- as.character(groups)
  fr_levs <- g_levs[re_rank]

  rows_out <- vector("list", length(fr_levs))

  for (ii in seq_along(fr_levs)) {
    lev <- fr_levs[ii]
    rows <- which(g_chr == lev)
    y_j  <- y[rows]
    X_j  <- Z[rows, , drop = FALSE]
    n_j  <- length(y_j)
    p_j  <- ncol(X_j)

    estimable <- FALSE
    note      <- character(0)

    if (n_j < 2L) {
      note <- "fewer than 2 observations"
    } else if (length(unique(y_j)) < 2L) {
      note <- "single outcome level"
    } else {
      df_j <- data.frame(y = y_j, X_j, check.names = FALSE)
      fit <- tryCatch(
        suppressWarnings(
          stats::glm(
            y ~ . - 1,
            data    = df_j,
            family  = family,
            control = stats::glm.control(maxit = 50L)
          )
        ),
        error = function(e) e
      )
      if (inherits(fit, "error")) {
        note <- conditionMessage(fit)
      } else {
        cf <- stats::coef(fit)
        if (length(cf) != p_j) {
          note <- sprintf(
            "glm returned %d coefficient(s), expected %d",
            length(cf), p_j
          )
        } else if (!isTRUE(fit$rank == p_j)) {
          note <- sprintf("rank-deficient glm fit (rank %d, expected %d)",
                          fit$rank, p_j)
        } else if (any(is.na(cf))) {
          note <- "NA coefficient(s)"
        } else if (any(!is.finite(cf))) {
          note <- "non-finite coefficient(s) (possible separation)"
        } else {
          V_ok <- tryCatch({
            V <- stats::vcov(fit)
            is.matrix(V) && all(is.finite(V))
          }, error = function(e) FALSE)
          if (!isTRUE(V_ok)) {
            note <- "vcov not finite (possible separation)"
          } else {
            estimable <- TRUE
          }
        }
      }
    }

    re_estimable[[lev]] <- estimable
    rows_out[[ii]] <- data.frame(
      group     = lev,
      n         = n_j,
      p         = p_j,
      re_rank   = TRUE,
      estimable = estimable,
      note      = if (length(note)) paste(note, collapse = "; ") else "",
      stringsAsFactors = FALSE
    )
  }

  re_glm_check <- if (length(rows_out)) {
    do.call(rbind, rows_out)
  } else {
    data.frame(
      group = character(0),
      n = integer(0),
      p = integer(0),
      re_rank = logical(0),
      estimable = logical(0),
      note = character(0),
      stringsAsFactors = FALSE
    )
  }

  list(re_estimable = re_estimable, re_glm_check = re_glm_check)
}
