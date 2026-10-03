#' Prior setup for generalized linear mixed models
#'
#' Calibrates priors for the level-2 fixed effects (\code{fixef}) of a
#' hierarchical mixed model using the reference \code{lmer}/\code{glmer} fits
#' on \strong{all} groups (from \code{\link{model_setup}}).  Per-group design
#' rank (\code{groupef.rank}) and estimability (\code{groupef.estimable} --
#' classical-\code{glm} MLE existence for \code{family =
#' binomial()}/\code{poisson()}/\code{Gamma()}, residual degrees of freedom
#' for \code{family = gaussian()}) are diagnostic checks computed by
#' \code{\link{model_setup}} and do not subset the data or the reference
#' \code{glmer} fit.  Random-effect variances are
#' treated as fixed at their mixed-model estimates.  The returned object
#' organizes calibrated priors as two independent blocks:
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
#'   Block~1 dispersion Gamma calibration is omitted (analogous to
#'   \code{\link[glmbayesCore]{Prior_Setup}} for flat GLMs).
#' @param REML Logical scalar passed to \code{\link{model_setup}} (and thus
#'   to the Gaussian reference \code{lmer} / \code{glmmTMB} fits).  Default
#'   \code{TRUE}.  Ignored for non-Gaussian families (\code{glmer} has no
#'   REML criterion).  Affects the classical estimates used for prior
#'   calibration only.
#' @param control Optional \code{\link[lme4]{lmerControl}} /
#'   \code{\link[lme4]{glmerControl}} object passed to
#'   \code{\link{model_setup}}.  When \code{NULL} (default), Gaussian models
#'   use \code{lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))}
#'   and non-Gaussian models leave \code{glmer}'s own control default.
#' @param start Optional starting values passed to \code{\link{model_setup}}
#'   for the reference \code{lmer}/\code{glmer} fit.  Default \code{NULL}.
#' @param verbose Integer scalar passed to \code{\link{model_setup}} for the
#'   reference \code{lmer}/\code{glmer} optimizer.  Default \code{0L}.
#' @param subset,na.action Optional arguments passed through to
#'   \code{\link{model_setup}} (same meanings as in \code{\link[lme4]{lmer}} /
#'   \code{\link[lme4]{glmer}}).  \code{model_setup()} applies them to both
#'   design extraction (\code{\link[lme4]{lFormula}}) and the reference
#'   \code{lmer}/\code{glmer}/\code{glmmTMB} fits, so sampler matrices
#'   (\code{design$y}/\code{D}/\code{W}) and classical calibration use the
#'   same rows.
#' @param weights,offset Optional arguments passed through to
#'   \code{\link{model_setup}} (same meanings as in \code{lmer}/\code{glmer}).
#'   Stored on \code{design$weights} / \code{design$offset} as length-\code{n}
#'   vectors aligned with \code{design$y} (after \code{subset}/
#'   \code{na.action}), for Phase~1 availability to samplers.
#'
#'   \strong{Not fully active for mixed-model sampling yet.} They affect the
#'   classical reference fit used for prior calibration. Low-level
#'   \code{two_block_rNormal_reg()} already accepts and echoes
#'   \code{prior.weights}/\code{offset2} when passed explicitly; high-level
#'   mixed routes that still hard-code unit weights / zero offset (and ICM /
#'   prior-weight algebra that assume unit weights) do not consume
#'   \code{design$weights}/\code{design$offset} yet. Prefer incorporating an
#'   offset in the formula or waiting until those paths are wired before
#'   relying on these for the Bayesian sampler.
#' @param contrasts Optional contrast list passed to
#'   \code{\link{model_setup}} (as in \code{lmer}).  Applied to both design
#'   extraction and the reference fit so factor coding in \code{W} matches
#'   \code{fixef}/\code{vcov}.
#' @param dispformula One-sided formula selecting the Block~1 measurement-
#'   dispersion structure: \code{~1} (default, pooled) or \code{~<group_name>}
#'   (matching the random-effects grouping factor exactly, requesting
#'   per-group dispersion).  Passed through to \code{\link{model_setup}},
#'   which fits and stores the \code{glmmTMB} reference
#'   (\code{design$glmmTMB_fit}) reused below as \code{fit_ref}.  Gaussian
#'   models only compute the per-group \code{group.ing_prior}
#'   calibration (a within-group regression fit for every group level, used
#'   only by \code{\link{dGamma_list}()}) when \code{dispformula} is the
#'   group formula; \code{dispformula = ~1} skips it entirely, so
#'   \code{group.dispersion.pwt} / \code{group.dispersion.nprior} must then be scalar
#'   (pooled).  Mirrors the \code{dispformula} argument on
#'   \code{lmerb()}/\code{glmerb()} in \pkg{lmebayes}, which gates the
#'   analogous choice of sampler route; the two are independent arguments
#'   that must be kept consistent by the caller.
#' @param pop.pwt Prior weight(s) in \eqn{(0, 1)} for the population-level
#'   (Block~2, \code{fixef}) coefficients.  Either a \strong{scalar} (applied
#'   to every random-effect component and every Block~2 predictor), or a
#'   \strong{list with one element per random-effect component} (named
#'   with the RE coefficient names in any order, or unnamed positional).
#'   Each list element is a scalar (recycled over that component's Block~2
#'   predictors) or a vector of length \eqn{p_k} (optionally named with the
#'   predictor column names of \code{X_hyper[[k]]}, reordered to match).
#'   The prior covariance for each \code{fixef_k} block is scaled relative to
#'   \code{vcov(fit_ref)} following the \code{\link[glmbayesCore]{Prior_Setup}}
#'   convention: \eqn{(1-\mathrm{pop.pwt})/\mathrm{pop.pwt}} for a scalar, and
#'   elementwise
#'   \eqn{\sqrt{(1-\mathrm{pop.pwt}_i)/\mathrm{pop.pwt}_i}\,
#'        \sqrt{(1-\mathrm{pop.pwt}_j)/\mathrm{pop.pwt}_j}} for vectors.  At
#'   most one of \code{pop.pwt} (explicitly supplied), \code{pop.nprior},
#'   and \code{pop.sd} may be supplied.  When \code{pop.pwt} is \code{NULL}
#'   (the default), dimension-adaptive defaults are chosen separately for
#'   each random-effect component (see \code{pop.pwt_default_low}/
#'   \code{pop.pwt_default_high}).
#' @param pop.pwt_default_low Scalar in \eqn{(0, 1)}, default \code{0.01}.
#'   Used when \code{pop.pwt} is \code{NULL} and a random-effect component's
#'   Block~2 predictor count \eqn{p_k = ncol(X\_hyper[[k]])} is strictly less
#'   than \code{14} (mirrors \code{\link[glmbayesCore]{Prior_Setup}}
#'   \code{pwt_default_low}).
#' @param pop.pwt_default_high Scalar in \eqn{(0, 1)}, default \code{0.05}.
#'   Used when \code{pop.pwt} is \code{NULL} and \eqn{p_k \ge 14} for that
#'   component (mirrors \code{\link[glmbayesCore]{Prior_Setup}}
#'   \code{pwt_default_high}).
#' @param pop.nprior Optional \emph{absolute} effective prior sample
#'   size(s), in group units (\eqn{J = }\code{nlevels(group)}), for the
#'   population-level (Block~2, \code{fixef}) coefficients -- an
#'   alternative to \code{pop.pwt}, on the same \eqn{n_i = J\,w_i/(1-w_i)}
#'   duality as \code{pop.dispersion.pwt}/\code{pop.dispersion.nprior}.
#'   Same shape as \code{pop.pwt}: a \strong{scalar}, or a \strong{list with
#'   one element per random-effect component} (named or unnamed
#'   positional), each element \code{NULL} (keep the \code{pop.pwt}-derived
#'   default for that component), a scalar, or a vector of length
#'   \eqn{p_k} (optionally named with the predictor column names of
#'   \code{X_hyper[[k]]}, reordered to match).  Converted internally to an
#'   effective weight \eqn{w_i = \mathrm{pop.nprior}_i /
#'   (\mathrm{pop.nprior}_i + J)}, used in place of \code{pop.pwt} for that
#'   component's \code{Sigma} scaling.  At most one of \code{pop.pwt}
#'   (explicitly supplied), \code{pop.nprior}, and \code{pop.sd} may be
#'   supplied.  The returned \code{pop.nprior} field (see \code{Value}
#'   below) is always populated with the effective prior sample size
#'   actually implied (derived from whichever of \code{pop.pwt}/
#'   \code{pop.sd} was used when \code{pop.nprior} is not supplied), in the
#'   same shape this argument expects -- call once without
#'   \code{pop.nprior} to inspect the defaults, then feed a modified copy
#'   of \code{ps$pop.nprior} back in to override.
#' @param pop.sd Optional prior standard deviation(s), on the coefficient
#'   scale of \code{vcov(fit_ref)}, for the population-level (Block~2,
#'   \code{fixef}) coefficients -- an alternative to \code{pop.pwt}.  Same
#'   shape as \code{pop.pwt}: a \strong{scalar}, or a \strong{list with one
#'   element per random-effect component} (named or unnamed positional),
#'   each element \code{NULL} (keep the \code{pop.pwt}-derived default for
#'   that component -- a \emph{partial} override across components), a
#'   scalar, or a vector of length \eqn{p_k} (optionally named with the
#'   predictor column names of \code{X_hyper[[k]]}, reordered to match).
#'   Converted internally to an effective weight
#'   \eqn{w_i = V_{ii} / (V_{ii} + \mathrm{pop.sd}_i^2)}, where \eqn{V_{ii}}
#'   is the corresponding diagonal entry of \code{vcov(fit_ref)}, and used
#'   in place of \code{pop.pwt} for that component's \code{Sigma}
#'   scaling.  At most one of \code{pop.pwt} (explicitly supplied),
#'   \code{pop.nprior}, and \code{pop.sd} may be supplied.  The returned
#'   \code{pop.sd} field (see \code{Value} below) is always populated with
#'   the effective standard deviations actually used
#'   (\eqn{\sqrt{\mathrm{diag}(\Sigma_{\mathrm{fixef},k})}}, derived from
#'   whichever of \code{pop.pwt}/\code{pop.nprior} was used when
#'   \code{pop.sd} is not supplied), in the same shape this argument
#'   expects -- call once without \code{pop.sd} to inspect the defaults,
#'   then feed a modified copy of \code{ps$pop.sd} back in to override.
#' @param pop.dispersion.pwt Optional \emph{relative} prior weight(s) in
#'   \eqn{(0, 1)} for the Block~2 dispersion (precision) prior, decoupled from
#'   \code{pop.pwt}.  A scalar, or a list / numeric vector with one value per
#'   random-effect component (named or positional).  Converted internally to
#'   an effective prior sample size \eqn{n_k = J\,w_k/(1-w_k)} where \eqn{J}
#'   is the number of groups.  At most one of \code{pop.dispersion.pwt} and
#'   \code{pop.dispersion.nprior} may be supplied; when neither is, the
#'   value is derived from \code{pop.pwt} (the mean across a component's
#'   predictors), keeping the Block~2 \eqn{\tau^2_k} dispersion prior
#'   consistent with the coefficient prior strength.  Weak values carry no
#'   computational penalty for \code{dIndependent_Normal_Gamma} sampling: the
#'   \eqn{\tau^2} truncation window comes from limiting-posterior quantiles
#'   independent of the prior strength (see \code{pop.ing_prior} below).
#' @param pop.dispersion.nprior Optional \emph{absolute} effective prior
#'   sample size(s) (in group units) for the Block~2 \eqn{\tau^2_k} dispersion
#'   prior.  A positive scalar, or a list / numeric vector with one value per
#'   random-effect component (named or positional).  See
#'   \code{pop.dispersion.pwt}.
#' @param pop.max_disp_perc Scalar in \eqn{(0.5, 1)}, default
#'   \code{0.99}, or a named/positional length-\eqn{p_{\mathrm{re}}} vector
#'   (one value per random-effect component).  Same tail-probability role as
#'   \code{group.max_disp_perc}, but for the Block~2 \eqn{\tau^2_k}
#'   truncation window stored in \code{pop.ing_prior}.  Only operationally
#'   consumed when the sampler's dispersion \code{pfamily} is
#'   \code{dIndependent_Normal_Gamma} (ignored for \code{dNormal}).
#' @param pop.intercept_source Character string controlling the prior mean for
#'   the global intercept hyperparameter \code{(Intercept)::(Intercept)} only.
#'   One of \code{"null_model"} (default) or \code{"full_model"}.  When
#'   \code{"null_model"}, the prior mean is taken from a random-intercept-only
#'   reference fit \code{y ~ 1 + (1 | group)} that omits all fixed-effect
#'   predictors (analogous to \code{\link[glmbayesCore]{Prior_Setup}} with
#'   \code{intercept_source = "null_model"}).  When \code{"full_model"}, the
#'   full-model MLE intercept is used.
#' @param pop.effects_source Character string controlling the prior mean for
#'   all other Block~2 hyperparameters (including population-mean slopes
#'   stored as \code{(Intercept)} columns in non-intercept RE components, and
#'   any non-intercept columns in \code{X_hyper}).  One of
#'   \code{"null_effects"} (default) or \code{"full_model"}.  When
#'   \code{"null_effects"}, prior means are set to zero.  When
#'   \code{"full_model"}, full-model MLE values are used.
#' @param pop.mu Optional prior mean override(s) for the population-level
#'   (Block~2, \code{fixef}) coefficients -- a finer-grained alternative to
#'   \code{pop.intercept_source}/\code{pop.effects_source}; prefer those
#'   first, and reach for \code{pop.mu} only when a specific mean needs to
#'   differ from either derived default.  Same shape as \code{pop.pwt}: a
#'   \strong{scalar}, or a \strong{list with one element per random-effect
#'   component} (named or unnamed positional), each element \code{NULL}
#'   (keep the \code{pop.intercept_source}/\code{pop.effects_source}-derived
#'   default mean for that component -- a \emph{partial} override across
#'   components), a scalar, or a vector of length \eqn{p_k} (optionally
#'   named with the predictor column names of \code{X_hyper[[k]]}, reordered
#'   to match).  Replaces \code{mu} wholesale for the components
#'   supplied; components left \code{NULL} keep their usual
#'   \code{pop.intercept_source}/\code{pop.effects_source} default.  The
#'   returned \code{pop.mu} field (see \code{Value} below) is always
#'   populated with the values actually used (derived defaults when
#'   \code{pop.mu} is not supplied), in the same shape this argument
#'   expects -- call once without \code{pop.mu} to inspect the defaults,
#'   then feed a modified copy of \code{ps$pop.mu} back in to override.
#' @param group.dispersion.pwt Optional relative prior weight(s) in \eqn{(0, 1)} for
#'   the Block~1 observation \eqn{\sigma^2} Gamma prior (Gaussian models only),
#'   decoupled from \code{pop.pwt} (Block~2 fixef) and from
#'   \code{pop.dispersion.pwt} (Block~2 \eqn{\tau^2_k}).  Either a
#'   \strong{scalar} converted to \eqn{n_{\mathrm{prior}} = w/(1-w)\times n}
#'   on the total observation count \eqn{n} (pooled \code{group.ing_prior}),
#'   or a \strong{named or positional vector of length \eqn{J}}
#'   (\eqn{J =} \code{nlevels(groups)}) with one weight per group level for
#'   per-group \code{dGamma_list()} calibration
#'   (\eqn{n_{\mathrm{prior},j} = w_j/(1-w_j)\times n_j}).  When
#'   \code{group.dispersion.pwt} is a vector, the pooled \code{group.ing_prior}
#'   continues to use the default \code{w = 0.01} on total \eqn{n} (unless
#'   \code{group.dispersion.nprior} is supplied explicitly).  At most one of
#'   \code{group.dispersion.pwt} and \code{group.dispersion.nprior} may be supplied; when neither
#'   is, scalar \code{group.dispersion.pwt = 0.01} is used for both pooled and
#'   per-group paths.
#' @param group.dispersion.nprior Optional positive scalar: absolute effective
#'   prior sample size for the Block~1 \eqn{\sigma^2} prior (observation units,
#'   not groups).  See \code{group.dispersion.pwt}.
#' @param group.alpha_target Scalar in \eqn{(0, 1)}, default
#'   \code{0.01}, or \code{NULL} to disable.  Gaussian models with per-group
#'   \code{dispformula} only: target percentage of posterior \code{beta_j}
#'   draws allowed to fall outside the per-group Mahalanobis/log-concavity
#'   ellipsoid of \code{inst/BLOCK_GIBBS_ERGODICITY_ING.md} Section 16.6.
#'   This is a different quantity from \code{group.dispersion.pwt}: \code{group.dispersion.pwt}
#'   is the Gamma prior \emph{weight} (an input to the prior), while
#'   \code{group.alpha_target} is a target \emph{violation rate} (an
#'   output being aimed at).  When non-\code{NULL} and \code{dispformula}
#'   requests per-group dispersion, \code{Prior_Setup_GLMM()} searches,
#'   for each group, for the smallest \code{group.dispersion.pwt} driving the
#'   predicted violation rate down to \code{group.alpha_target},
#'   \strong{floored at} the \code{group.dispersion.pwt} resolved from the
#'   \code{group.dispersion.pwt}/\code{group.dispersion.nprior} arguments above (so
#'   calibration only ever sharpens, never loosens, the per-group prior;
#'   supplying \code{group.dispersion.pwt} explicitly still guarantees a minimum).
#'   Active by default as a safeguard against poorly-behaved per-group
#'   priors; set to \code{NULL} to opt out and use \code{group.dispersion.pwt}/
#'   \code{group.dispersion.nprior} directly, unmodified.  Silently ignored
#'   (not an error) when left at its default on a model that does not
#'   support per-group calibration (e.g.\ \code{dispformula = ~1}); it is
#'   only an error to \emph{explicitly} supply a non-\code{NULL} value in
#'   that case.  See \code{group.pwt_calibration} below.
#' @param group.max_disp_perc Scalar in \eqn{(0.5, 1)}, default
#'   \code{0.99}, or a named/positional length-\eqn{J} vector (one value per
#'   group level, Gaussian models with per-group \code{dispformula} only).
#'   Tail probability used to compute the Block~1 \eqn{\sigma^2} truncation
#'   window stored in \code{group.ing_prior} -- the pooled scalar default
#'   when a vector is supplied here (\code{dispformula = ~1}), or the
#'   resolved per-group value (per group, per-group \code{dispformula}).
#'   The window is the central
#'   \eqn{2 \times \mathrm{group.max\_disp\_perc} - 1} mass interval:
#'   \code{disp_lower} = \eqn{(1-\mathrm{group.max\_disp\_perc})}
#'   quantile and \code{disp_upper} = \eqn{\mathrm{group.max\_disp\_perc}}
#'   quantile of the relevant Gamma (precision) distribution, inverted to the
#'   dispersion scale.  Tighter values (e.g.\ \code{0.95}) shrink the
#'   truncation window and typically improve Block~1 acceptance rates at the
#'   cost of slightly less envelope coverage; looser values (e.g.\ \code{0.999})
#'   widen the window.  Passed to \code{\link[glmbayesCore]{dGamma}()} as
#'   \code{max_disp_perc}.  A per-group vector lets outlier groups be given a
#'   tighter (or looser) window than the rest; see \code{dGamma_list()}.
#'
#' @return Object of class \code{"Prior_Setup_GLMM"} with fields:
#'   \describe{
#'     \item{\code{formula}}{Model formula.}
#'     \item{\code{family}}{Family object.}
#'     \item{\code{pop.pwt}}{Prior weight(s) used: the scalar as supplied, the
#'       canonical named list of per-predictor weight vectors when a list was
#'       supplied, or the dimension-adaptive default list (one weight per
#'       Block~2 predictor, per RE component) when \code{pop.pwt} was
#'       \code{NULL}.}
#'     \item{\code{pop.nprior}}{Named list with one named numeric vector per
#'       random-effect component -- the effective prior sample size (group
#'       units) actually implied, \eqn{n_i = J\,w_i/(1-w_i)}, whether
#'       derived from \code{pop.pwt}/\code{pop.sd} (the usual case) or
#'       overridden via \code{pop.nprior}.  Always present, and in the same
#'       shape the \code{pop.nprior} argument expects, so it can be
#'       inspected from a first call and fed back in (after editing) as
#'       \code{pop.nprior} to a second call.}
#'     \item{\code{pop.sd}}{Named list with one named numeric vector per
#'       random-effect component -- the effective coefficient-scale prior
#'       standard deviations actually used
#'       (\eqn{\sqrt{\mathrm{diag}(\Sigma_{\mathrm{fixef},k})}}), whether
#'       derived from \code{pop.pwt}/\code{pop.nprior} (the usual case) or
#'       overridden via \code{pop.sd}.  Always present, and in the same
#'       shape the \code{pop.sd} argument expects, so it can be inspected
#'       from a first call and fed back in (after editing) as \code{pop.sd}
#'       to a second call.}
#'     \item{\code{pop.mu}}{Named list with one named numeric vector per
#'       random-effect component -- the prior means actually used
#'       (\code{mu}), whether derived from
#'       \code{pop.intercept_source}/\code{pop.effects_source} (the usual
#'       case) or overridden via \code{pop.mu}.  Always present, and in the
#'       same shape the \code{pop.mu} argument expects, so it can be
#'       inspected from a first call and fed back in (after editing) as
#'       \code{pop.mu} to a second call.}
#'     \item{\code{pop.dispersion.pwt}}{Named per-component vector of relative
#'       dispersion prior weights (always present; consistent with
#'       \code{pop.dispersion.nprior} via \eqn{w_k = n_k/(n_k + J)}).}
#'     \item{\code{pop.dispersion.nprior}}{Named per-component vector of
#'       effective prior sample sizes for the Block~2 dispersion prior
#'       (always present; used by
#'       \code{\link[=pfamily_list.Prior_Setup_GLMM]{pfamily_list}()} to
#'       calibrate \code{dIndependent_Normal_Gamma} components).}
#'     \item{\code{group.dispersion.pwt}}{Gaussian models only: scalar or length-\eqn{J}
#'       vector of relative prior weights for Block~1 \eqn{\sigma^2} -- a
#'       scalar unless a per-group vector was supplied or
#'       \code{group.alpha_target} calibration was active, in which case it
#'       is the resolved length-\eqn{J} vector (see \code{dGamma_list()}).}
#'     \item{\code{group.dispersion.nprior}}{Gaussian models only: effective prior
#'       sample size for Block~1 \eqn{\sigma^2}, on the same pooled-vs.-vector
#'       basis as \code{group.dispersion.pwt} (a scalar on the observation scale, or a
#'       named length-\eqn{J} vector of per-group \eqn{n_{\mathrm{prior},j}}
#'       reflecting the \code{group.alpha_target}-calibrated \code{group.dispersion.pwt}
#'       when calibration was active; see \code{group.pwt_calibration}).}
#'     \item{\code{group.alpha_target}}{The value used, or \code{NULL}
#'       if calibration was not active.}
#'     \item{\code{group.pwt_calibration}}{\code{NULL} unless
#'       calibration was active; otherwise a per-group diagnostics
#'       \code{data.frame} (\code{group}, \code{n_j}, \code{pwt_floor},
#'       \code{pct_outside_before}, \code{w_star}, \code{clipped_at_ceiling},
#'       \code{w_final}, \code{floor_binds}, \code{pct_outside_after}) from
#'       the search described under \code{group.alpha_target}.}
#'     \item{\code{group.max_disp_perc}}{Gaussian models only: scalar or
#'       length-\eqn{J} vector of the tail probability used for the Block~1
#'       \eqn{\sigma^2} truncation window (\code{group.ing_prior}) -- the
#'       scalar as supplied, or the resolved length-\eqn{J} vector when a
#'       vector was supplied.}
#'     \item{\code{pop.max_disp_perc}}{Named length-\eqn{p_{\mathrm{re}}}
#'       vector of resolved per-component values used for the Block~2
#'       \eqn{\tau^2_k} truncation windows in \code{pop.ing_prior} (always
#'       present).}
#'     \item{\code{pop.intercept_source}}{How the global intercept prior mean
#'       was chosen (\code{"null_model"} or \code{"full_model"}).}
#'     \item{\code{pop.effects_source}}{How other Block~2 prior means were
#'       chosen (\code{"null_effects"} or \code{"full_model"}).}
#'     \item{\code{block_formula}}{Within-group Block~1 formula: response
#'       regressed on the random-coefficient structure only (columns of
#'       \code{design$groupef.names}); level-2 hyper covariates are excluded.
#'       Used by \code{dGamma_list()}.}
#'     \item{\code{group.tau_sd}}{Named vector \code{sqrt(Psi)} from the reference
#'       fit; shared population RE standard deviations for per-group calibration.}
#'     \item{\code{data}}{Data frame passed to \code{Prior_Setup_GLMM()}
#'       (reference for \code{dGamma_list()} diagnostics).}
#'     \item{\code{design}}{Full \code{\link{model_setup}} object (all groups),
#'       including length-\code{n} \code{weights} and \code{offset} vectors
#'       aligned with \code{design$y}.}
#'     \item{\code{mer_fit}}{Reference \code{lmer}/\code{glmer} fit on all
#'       groups (the full-formula fit from \code{\link{model_setup}}), always
#'       present regardless of \code{dispformula} (backs
#'       \code{group.dispersion} and \code{x$lmer}/\code{x$glmer} in
#'       \code{lmerb()}/\code{glmerb()}).}
#'     \item{\code{fit_ref}}{The calibration reference for Block~2
#'       (\code{fixef}/\eqn{\tau^2_k}) and the per-group Block~1 inputs
#'       (\code{group.tau_sd}, BLUP coefficients): identical to \code{mer_fit} when
#'       \code{dispformula = ~1}; otherwise an equivalent \code{glmmTMB}
#'       fit with the same \code{dispformula} (Gaussian models only). See
#'       \code{calibration_source}.}
#'     \item{\code{group.dispersion.fit}}{\code{NULL} when \code{dispformula = ~1};
#'       otherwise the same \code{glmmTMB} object as \code{fit_ref}.}
#'     \item{\code{group.dispersion.ref}}{\code{NULL} when \code{dispformula = ~1};
#'       otherwise a named length-\eqn{J} vector of per-group observation-level
#'       dispersion read from \code{fit_ref}'s dispersion linear predictor
#'       (\code{glmmTMB::predict(fit_ref, type = "disp")}, aggregated by
#'       group). Diagnostic only -- not the value fed to the sampler; compare
#'       against \code{group.ing_prior}'s \code{sigma2_hat}.}
#'     \item{\code{calibration_source}}{\code{"lme4"} or \code{"glmmTMB"};
#'       which reference fit produced \code{fit_ref} (and therefore
#'       \code{pop.prior_list}, \code{pop.ing_prior}, \code{group.tau_sd}, and
#'       \code{group.ing_prior}).}
#'     \item{\code{group.dispersion}}{For Gaussian models: resolved Block~1
#'       observation \eqn{\sigma^2} used to mean-match the Gamma prior.
#'       When \code{dispformula = ~1}, the pooled A12 §3.3.4 aggregate
#'       \eqn{\hat\sigma^2 =
#'       (\sum_j S_{\mathrm{marg},j})/(n - J p_{\mathrm{re}})} (Part VI
#'       \code{Omega_j} included); classical lmer residual in
#'       \code{attr(., "classical")}. When \code{dispformula} is per-group,
#'       a named length-\eqn{J} vector of per-group
#'       \code{group.ing_prior[[j]]$sigma2_hat}. \code{attr(., "source")}
#'       records provenance.  \code{NULL} for non-Gaussian families.
#'       Pass to \code{rlmerb()}/\code{rglmerb()} as
#'       \code{group.dispersion} for fixed-\eqn{\sigma^2} routes.}
#'     \item{\code{group.Sigma}}{Diagonal RE covariance matrix (Block~1),
#'       classical reference-fit \eqn{\Psi} diagonals.}
#'     \item{\code{pop.dispersion.ref}}{Named vector of classical reference-fit
#'       \eqn{\hat\tau^2_k} (VarCorr / \code{Psi}) before the A12 §3.3.4
#'       hyper-regression remargin; compare to
#'       \code{pop.prior_list[[k]]$dispersion}.}
#'     \item{\code{pop.prior_list}}{Named Block~2 prior list per RE coefficient.
#'       For Gaussian models \code{dispersion} is the A12 §3.3.4 marginal
#'       hyper-regression center for \eqn{\tau^2_k} (not the raw VarCorr
#'       plug-in; see \code{pop.dispersion.ref}).}
#'     \item{\code{pop.ing_prior}}{Named per-component list of the prospective
#'       \code{dIndependent_Normal_Gamma} calibration: Gamma precision-prior
#'       \code{shape} \eqn{= (n_0 + 1 + p_k)/2} and \code{rate}
#'       \eqn{= \hat\tau^2_k (n_0 + p_k - 1)/2} (the lmebayesCore default
#'       calibration with \eqn{n_0 =} \code{pop.dispersion.nprior}; since
#'       \code{rate} \eqn{= \hat\tau^2_k (\code{shape} - 1)}, the implied
#'       inverse-Gamma prior on \eqn{\tau^2_k} has mean exactly
#'       \eqn{\hat\tau^2_k}), where \eqn{\hat\tau^2_k} is the A12 §3.3.4
#'       marginal hyper-regression center in
#'       \code{pop.prior_list[[k]]$dispersion}, and the default \eqn{\tau^2_k}
#'       truncation window \code{disp_lower} / \code{disp_upper}: the
#'       \eqn{(1-\mathrm{pop.max\_disp\_perc})} /
#'       \eqn{\mathrm{pop.max\_disp\_perc}}
#'       quantiles of the \emph{limiting posterior}
#'       \eqn{\Gamma((J+1)/2,\; \hat\tau^2_k (J-1)/2)} -- the weak-prior
#'       (\eqn{n_0 \to 0}) limit of the Block~2 posterior Gamma for the
#'       precision (lmebayesCore Chapter A12, Theorem 2; inverted to a
#'       \eqn{\tau^2} interval).  This window is identical for all
#'       \eqn{n_0}, covers \eqn{\ge} \eqn{2 \times \mathrm{pop.max\_disp\_perc} - 1}
#'       of the exact posterior for every prior strength, and keeps the
#'       envelope sampler's cost stable as priors weaken; see
#'       \code{inst/ING_TRUNCATION_WINDOW.md}.  Used by
#'       \code{\link[=pfamily_list.Prior_Setup_GLMM]{pfamily_list}()} when
#'       \code{ptypes = "dIndependent_Normal_Gamma"}; ignored for
#'       \code{dNormal} priors.}
#'     \item{\code{group.ing_prior}}{Gaussian models only; \code{NULL}
#'       for non-Gaussian models. Shape depends on \code{dispformula}:
#'       \describe{
#'         \item{\code{dispformula = ~1} (pooled)}{A single list: prospective
#'           \code{dGamma()} \code{group.dispersion} calibration for Block~1
#'           ING (observation \eqn{\sigma^2} shared across all group
#'           levels) -- mean-matched \code{shape} / \code{rate} with
#'           \eqn{n_{\mathrm{prior}} = \mathrm{group.dispersion.pwt}/(1-\mathrm{group.dispersion.pwt})\times n},
#'           \eqn{p = p_{\mathrm{re}}}, and \eqn{\hat\sigma^2} =
#'           \code{group.dispersion} (pooled A12 §3.3.4 aggregate of
#'           per-group \eqn{S_{\mathrm{marg},j}} by default; same ING algebra
#'           as \code{pop.ing_prior} for \eqn{\tau^2_k}; requires
#'           \code{group.dispersion.pwt} \eqn{\le 0.5}),
#'           plus \code{disp_lower} / \code{disp_upper} as the central
#'           \eqn{2 \times \mathrm{group.max\_disp\_perc} - 1} prior-mass
#'           interval from the same calibrated \code{shape}/\code{rate}
#'           (always the pooled/scalar \code{group.max_disp_perc} value,
#'           even when a per-group vector was supplied). Pass these fields
#'           to \code{\link[glmbayesCore]{dGamma}()}.}
#'         \item{\code{dispformula} requesting per-group dispersion}{A named
#'           list (one entry per group level) of per-group \code{dGamma()}
#'           density calibration (\code{sigma2_hat}, \code{shape_ING},
#'           \code{rate}, \code{rate_gamma}, \code{n_prior}, \code{n_j},
#'           \code{n_combined}, \code{disp_lower}, \code{disp_upper},
#'           \code{max_disp_perc}, \code{omega_j}, \ldots).
#'           \code{shape_ING}/\code{rate} already fold in the Part VI
#'           model-derived \code{Omega_j} (fixed-effect/\code{gamma}
#'           uncertainty about \code{b_j}'s prior mean, propagated through
#'           group \code{j}'s hyper-design row and
#'           \code{pop.prior_list[[k]]$Sigma}), so this marginal
#'           integrates out both random effects (\code{b_j}) and fixed
#'           effects (\code{gamma}); see
#'           \code{inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md} Part VI.
#'           \code{disp_lower}/\code{disp_upper} are \eqn{(1-w_j)}/\eqn{w_j}
#'           quantiles of \code{Gamma(shape_ING + n_j/2, rate_post)}, where
#'           \code{rate_post} is \code{sigma2_hat} mean-matched at that same
#'           inflated shape (\eqn{n_j/2} added only for the window, so it
#'           reflects the posterior spread the sampler's own envelope
#'           machinery actually draws \eqn{\sigma^2_j} from -- see
#'           \code{EnvelopeDispersionBuild.cpp}'s \code{shape2 = Shape + n_w/2}
#'           fallback -- while \code{shape_ING}/\code{rate} themselves, the
#'           values actually fed to the sampler as the Gamma prior, are
#'           unchanged); \eqn{w_j} is this group's own resolved
#'           \code{group.max_disp_perc} value (stored per-group as
#'           \code{max_disp_perc}). Used directly by
#'           \code{\link{dGamma_list}()}.}
#'       }
#'       Which shape is returned is fully determined by \code{dispformula};
#'       see \code{\link{dGamma_list}()} for how consumers tell the two
#'       apart.}
#'     \item{\code{dispformula}}{The \code{dispformula} supplied.}
#'   }
#' @details
#' \strong{Naming convention.} Top-level fields use \code{pop.*} for Block~2
#' (population / \eqn{\gamma}, \eqn{\tau^2_k}) and \code{group.*} for Block~1
#' (within-group / observation \eqn{\sigma^2}, RE covariance \eqn{\Psi}).
#' Nested entries in \code{pop.prior_list[[k]]} use glmbayes-style
#' \code{mu}, \code{Sigma}, and \code{dispersion} (where \code{dispersion}
#' is the plug-in RE variance \eqn{\tau^2_k}).
#'
#' \strong{Why default calibration depends on classical estimates.}
#' \code{Prior_Setup_GLMM} scales Block~2 covariances from
#' \code{vcov(fit_ref)} by \eqn{(1-\mathrm{pwt})/\mathrm{pwt}} and plugs in
#' RE variances from the full reference fit, where \code{fit_ref} is the
#' pooled \code{lmer}/\code{glmer} fit when \code{dispformula = ~1}, or an
#' equivalent \code{glmmTMB} fit with the same \code{dispformula} when
#' per-group dispersion is requested (see \code{calibration_source} and
#' \code{mer_fit} above).  By default the global intercept
#' prior mean comes from a random-intercept-only null fit; all other prior
#' means are zero (\code{pop.effects_source = "null_effects"}).  This requires:
#' \enumerate{
#'   \item Converged reference \code{lmer}/\code{glmer} fit from
#'     \code{\link{model_setup}} on the full formula (and a random-intercept-only
#'     null fit when \code{pop.intercept_source = "null_model"}).
#'     Fits with \code{lme4} \code{checkConv} failures (e.g.\ large
#'     \code{max|grad|}) are rejected.
#'   \item Every \code{X_hyper[[k]]} column maps to a \code{fixef(fit_ref)} term.
#'   \item Each RE variance \eqn{\tau^2_k} from the reference fit is strictly positive.
#' }
#' @seealso \code{\link{model_setup}}, \code{\link{Prior_Setup}},
#'   \code{build_mu_all()},
#'   \code{\link{pfamily_list}}, \code{\link{dGamma_list}}
#' @example inst/examples/Ex_Prior_Setup_GLMM.R
#' @export
Prior_Setup_GLMM <- function(formula,
                                 data,
                                 family = gaussian(),
                                 REML = TRUE,
                                 control = NULL,
                                 start = NULL,
                                 verbose = 0L,
                                 subset,
                                 weights,
                                 na.action,
                                 offset,
                                 contrasts = NULL,
                                 dispformula = ~1,
                                 pop.pwt    = NULL,
                                 pop.pwt_default_low  = 0.01,
                                 pop.pwt_default_high = 0.05,
                                 pop.nprior = NULL,
                                 pop.sd     = NULL,
                                 pop.dispersion.pwt = NULL,
                                 pop.dispersion.nprior = NULL,
                                 pop.max_disp_perc = 0.99,
                                 pop.intercept_source = c("null_model", "full_model"),
                                 pop.effects_source   = c("null_effects", "full_model"),
                                 pop.mu = NULL,
                                 group.dispersion.pwt = NULL,
                                 group.dispersion.nprior = NULL,
                                 group.alpha_target = 0.01,
                                 group.max_disp_perc = 0.99) {

  intercept_source <- match.arg(pop.intercept_source)
  effects_source   <- match.arg(pop.effects_source)

  ## Internal aliases: the calibration logic below is unchanged and keeps
  ## using these short local names throughout; only the public argument
  ## names (pop.*/group.*) and the returned field names follow the
  ## pop./group. convention (Block 2/population vs. Block 1/per-group).
  pwt                       <- pop.pwt
  pwt_dispersion            <- pop.dispersion.pwt
  n_prior_dispersion        <- pop.dispersion.nprior
  pwt_measurement           <- group.dispersion.pwt
  n_prior_measurement       <- group.dispersion.nprior
  alpha_target_measurement  <- group.alpha_target
  max_disp_perc_measurement <- group.max_disp_perc
  max_disp_perc_dispersion  <- pop.max_disp_perc

  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula.", call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.", call. = FALSE)
  }
  check_pwt_scalar <- function(x, what) {
    if (!is.numeric(x) || length(x) != 1L || is.na(x) || x <= 0 || x >= 1) {
      stop(sprintf("'%s' must be a scalar in (0, 1).", what), call. = FALSE)
    }
  }
  check_pwt_scalar(pop.pwt_default_low,  "pop.pwt_default_low")
  check_pwt_scalar(pop.pwt_default_high, "pop.pwt_default_high")
  if (!is.null(pwt)) {
    if (is.numeric(pwt)) {
      check_pwt_scalar(pwt, "pop.pwt")
    } else if (!is.list(pwt)) {
      stop(
        "'pop.pwt' must be a scalar in (0, 1), NULL, or a list with one element per ",
        "random-effect component.",
        call. = FALSE
      )
    }
  }
  pwt_alt_supplied <- c(!is.null(pop.pwt), !is.null(pop.sd), !is.null(pop.nprior))
  if (sum(pwt_alt_supplied) > 1L) {
    stop(
      "Supply at most one of 'pop.pwt', 'pop.sd', and 'pop.nprior'.",
      call. = FALSE
    )
  }
  if (!is.null(pwt_dispersion) && !is.null(n_prior_dispersion)) {
    stop(
      "Supply at most one of 'pop.dispersion.pwt' and 'pop.dispersion.nprior'.",
      call. = FALSE
    )
  }
  if (!is.null(pwt_measurement) && !is.null(n_prior_measurement)) {
    stop(
      "Supply at most one of 'group.dispersion.pwt' and 'group.dispersion.nprior'.",
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

  if (is.null(control)) {
    if (is_gaussian) {
      control <- lme4::lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
      )
    }
    ## Non-Gaussian: leave NULL so glmer uses its own defaults (explicit
    ## glmerControl() can crash some lme4 builds).
  }
  ctrl <- control
  ## Evaluate subset against data (as model.frame / lmer do) before
  ## do.call(model_setup, ...): forcing the formal or re-evaluating a
  ## language object via do.call would look up column names in the wrong
  ## environment. weights/offset stay force-evaluated (concrete vectors).
  cl <- match.call()

  ms_args <- list(
    formula     = formula,
    data        = data,
    family      = family,
    REML        = REML,
    control     = control,
    verbose     = verbose,
    contrasts   = contrasts,
    dispformula = dispformula
  )
  if (!is.null(start)) {
    ms_args$start <- start
  }
  if (!is.null(cl$subset)) {
    ms_args$subset <- eval(cl$subset, envir = data, enclos = parent.frame())
  }
  if (!missing(weights)) {
    ms_args$weights <- weights
  }
  if (!missing(na.action)) {
    ms_args$na.action <- na.action
  }
  if (!missing(offset)) {
    ms_args$offset <- offset
  }
  design <- do.call(model_setup, ms_args)

  ## model_setup() already validated 'dispformula' (including the
  ## per-group-requires-gaussian check) and, when it requests per-group
  ## dispersion for a Gaussian model, fit and stored the glmmTMB reference
  ## as design$glmmTMB_fit -- reused as fit_ref below instead of fitting a
  ## second time here.
  dispformula_kind <- .lmebayes_dispformula_kind(
    dispformula, design$group_name
  )
  if (identical(dispformula_kind, "pooled") && is_gaussian &&
      (length(pwt_measurement) > 1L || length(n_prior_measurement) > 1L ||
       length(max_disp_perc_measurement) > 1L)) {
    stop(
      "'group.dispersion.pwt'/'group.dispersion.nprior'/'group.max_disp_perc' ",
      "has more than one value (per-group), but dispformula = ~1 (pooled); ",
      "use dispformula = ~", design$group_name, " to calibrate per-group ",
      "measurement-dispersion priors for dGamma_list().",
      call. = FALSE
    )
  }

  ## design$groupef.estimable / design$groupef.glm_check are already
  ## populated by model_setup() (classical-glm MLE existence check for
  ## binomial/poisson/Gamma, residual-df-for-dispersion check for gaussian;
  ## other families mirror design$groupef.rank) -- no need to recompute
  ## them here.

  ## Full-rank status is a per-group DESIGN CHECK only (reported by print();
  ## groups with rank-deficient Z_j are still fully used below).  The lme4
  ## reference fit is always fit by model_setup() (backward-compat: it backs
  ## x$lmer/x$glmer in lmerb()/glmerb() and the pooled group_dispersion
  ## scalar, regardless of dispformula).
  mer_fit <- if (is_gaussian) design$lmer else design$glmer
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
      "Prior_Setup_GLMM() requires converged ", mer_label,
      " reference fits:\n  - ",
      paste(mer_issues, collapse = "\n  - "),
      "\n\nRevise the model or supply hyperpriors manually without ",
      "Prior_Setup_GLMM().",
      call. = FALSE
    )
  }

  group_levels <- levels(design$group)

  ## Placeholder until marginal centers overwrite below (pooled S_marg or
  ## per-group sigma2_hat). Classical lmer residual kept as diagnostic.
  group_dispersion <- if (is_gaussian) {
    val <- as.numeric(design$dispersion)
    attr(val, "source") <- "lmer reference fit"
    val
  } else {
    NULL
  }

  ## Calibration reference: when dispformula requests per-group dispersion,
  ## Block~2 (fixef/tau^2_k) and the per-group Block~1 inputs (sd_tau, BLUP
  ## coefficients) come from an equivalent glmmTMB fit instead of the pooled
  ## lme4 fit above, so the heteroscedastic structure that dispformula
  ## requests is reflected in every calibrated quantity, not only in
  ## dGamma_list()'s per-group densities. See
  ## inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md.
  calibration_source <- "lme4"
  group_dispersion_ref   <- NULL
  if (is_gaussian && identical(dispformula_kind, "group")) {
    fit_ref <- design$glmmTMB_fit
    if (is.null(fit_ref)) {
      stop(
        "model_setup() did not return a glmmTMB reference fit ",
        "(design$glmmTMB_fit) for dispformula = ~", design$group_name, ".",
        call. = FALSE
      )
    }
    tmb_issues <- .lmebayes_glmmtmb_convergence_issues(
      fit_ref, "glmmTMB (full formula, per-group dispersion)"
    )
    if (length(tmb_issues) > 0L) {
      stop(
        "Prior_Setup_GLMM() requires a converged glmmTMB reference fit ",
        "for dispformula = ~", design$group_name, ":\n  - ",
        paste(tmb_issues, collapse = "\n  - "),
        "\n\nRevise the model or supply hyperpriors manually without ",
        "Prior_Setup_GLMM().",
        call. = FALSE
      )
    }
    vc_ref   <- extract_glmmtmb_variance_components(
      fit_ref, design$groupef.names, design$group_name
    )
    tau2_vec <- vc_ref$Psi
    group_dispersion_ref <- .lmebayes_glmmtmb_group_sigma2(
      fit_ref, design$group_name, levels(design$group)
    )
    calibration_source <- "glmmTMB"
  } else {
    fit_ref  <- mer_fit
    tau2_vec <- design$Psi
  }
  ref_label <- if (identical(calibration_source, "glmmTMB")) "glmmTMB" else mer_label

  p_re        <- length(design$groupef.names)
  group.Sigma <- diag(unname(tau2_vec), nrow = p_re, ncol = p_re)
  dimnames(group.Sigma) <- list(design$groupef.names, design$groupef.names)

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

  re_names  <- design$groupef.names
  tau_tol   <- sqrt(.Machine$double.eps)
  re_issues <- character(0)

  for (k in re_names) {
    X_k    <- design$W[[k]]
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
      "Prior_Setup_GLMM() cannot calibrate default hyperpriors:\n  - ",
      paste(re_issues, collapse = "\n  - "),
      "\n\nRevise the formula (e.g. add a fixed main effect for each random ",
      "slope and avoid RE terms with zero estimated variance), or supply ",
      "hyperpriors manually without Prior_Setup_GLMM().",
      call. = FALSE
    )
  }

  if (is.null(pwt)) {
    pwt <- .lmebayes_default_pwt_list(
      design,
      pwt_default_low  = pop.pwt_default_low,
      pwt_default_high = pop.pwt_default_high
    )
    if (is.null(pop.sd) && is.null(pop.nprior)) {
      for (k in re_names) {
        p_k <- length(colnames(design$W[[k]]))
        w_k <- unname(pwt[[k]][1L])
        label <- if (p_k >= 14L) "high-d" else "low-d"
        message(
          "Using default pop.pwt = ", w_k, " for ", k,
          " (", label, " default, p_k = ", p_k, ")."
        )
      }
    }
  }

  pwt_list <- .lmebayes_resolve_pwt(pwt, design)
  mu_override_list     <- .lmebayes_resolve_mu_override(pop.mu, design)
  nprior_override_list <- .lmebayes_resolve_nprior_override(pop.nprior, design)
  sd_override_list     <- .lmebayes_resolve_sd_override(pop.sd, design)

  J_groups <- nlevels(design$group)
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
    ## Same REML / control / optional mer args as the full-model reference
    ## so intercept calibration sees the same row set and fitting options.
    null_opt <- list()
    if (!is.null(start)) {
      null_opt$start <- start
    }
    if (!is.null(cl$subset)) {
      ## Same concrete index as the full-model reference / model_setup.
      null_opt$subset <- ms_args$subset
    }
    if (!missing(weights)) {
      null_opt$weights <- weights
    }
    if (!missing(na.action)) {
      null_opt$na.action <- na.action
    }
    if (!missing(offset)) {
      null_opt$offset <- offset
    }
    if (!is.null(contrasts)) {
      null_opt$contrasts <- contrasts
    }

    if (identical(calibration_source, "glmmTMB")) {
      ## Match the heteroscedastic reference so the global intercept prior
      ## mean is calibrated consistently with fit_ref above.
      null_fit <- do.call(
        .lmebayes_fit_glmmtmb_reference,
        c(
          list(
            formula     = null_formula,
            data        = data,
            family      = family,
            dispformula = dispformula,
            REML        = REML
          ),
          null_opt
        )
      )
      null_issues <- .lmebayes_glmmtmb_convergence_issues(
        null_fit, "glmmTMB (null intercept model, per-group dispersion)"
      )
    } else {
      null_args <- c(
        list(formula = null_formula, data = data),
        null_opt
      )
      if (is_gaussian) {
        null_args$REML <- REML
        if (!is.null(ctrl)) {
          null_args$control <- ctrl
        }
        null_fit <- do.call(lme4::lmer, null_args)
      } else {
        null_args$family <- family
        if (!is.null(ctrl)) {
          null_args$control <- ctrl
        }
        null_fit <- do.call(lme4::glmer, null_args)
      }
      null_issues <- .lmebayes_mer_convergence_issues(
        null_fit, sprintf("%s (null intercept model)", mer_label)
      )
    }
    if (length(null_issues) > 0L) {
      stop(
        "Prior_Setup_GLMM() requires a converged ", ref_label,
        " random-intercept-only null fit for intercept_source = \"null_model\":\n  - ",
        paste(null_issues, collapse = "\n  - "),
        "\n\nUse intercept_source = \"full_model\" or revise the model.",
        call. = FALSE
      )
    }
    fe_null <- .lmebayes_reference_fixef(null_fit)
  }

  prior_raw <- stats::setNames(
    lapply(re_names, function(k) {
      X_k    <- design$W[[k]]
      cols_k <- colnames(X_k)
      p_k    <- length(cols_k)
      tau2_k <- tau2_vec[[k]]

      fe_nms <- vapply(cols_k, fe_name_for, character(1L), k = k)
      fe_idx <- fe_nms

      mu <- vapply(seq_len(p_k), function(i) {
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
      names(mu) <- cols_k

      ## 'pop.mu' override: replaces the intercept_source/effects_source-
      ## derived default wholesale for this component only; components left
      ## NULL in 'pop.mu' keep the default computed above.
      if (!is.null(mu_override_list[[k]])) {
        mu <- mu_override_list[[k]]
      }

      ## 'pop.nprior'/'pop.sd' overrides: convert to an effective weight and
      ## substitute it for pwt_list[[k]] in the scaling below (at most one of
      ## 'pop.pwt' (explicit), 'pop.nprior', and 'pop.sd' may be supplied,
      ## enforced above); a component left NULL in the supplied alternative
      ## keeps its 'pop.pwt'-derived weight.
      pwt_k <- pwt_list[[k]]
      if (!is.null(nprior_override_list[[k]])) {
        pwt_k <- nprior_override_list[[k]] / (nprior_override_list[[k]] + J_groups)
      } else if (!is.null(sd_override_list[[k]])) {
        V_diag_k <- diag(V_fe[fe_idx, fe_idx, drop = FALSE])
        pwt_k    <- V_diag_k / (V_diag_k + sd_override_list[[k]]^2)
      }
      ## 'V_fe[fe_idx, fe_idx]' is indexed by *fixef* names (fe_idx, e.g.
      ## "Days" for the population-mean-slope column of a non-intercept RE
      ## component), not by the hyper-design column names (cols_k, e.g.
      ## "(Intercept)"); R's elementwise arithmetic above takes names
      ## positionally from its operands, so pwt_k can inherit the wrong
      ## (fixef) names.  Re-align explicitly to cols_k -- the convention
      ## every other per-component vector here (mu, Sigma
      ## dimnames, pwt_used) uses.
      names(pwt_k) <- cols_k

      ## Elementwise scaling sqrt(s_i) * sqrt(s_j) with s_i = (1-w_i)/w_i;
      ## reduces to the scalar (1-pwt)/pwt factor when all weights are equal.
      sc_k <- sqrt((1 - pwt_k) / pwt_k)
      Sigma <- V_fe[fe_idx, fe_idx, drop = FALSE] * outer(sc_k, sc_k)
      ## V_fe (vcov(fit_ref)) is not always exactly symmetric to floating-point
      ## precision (sandwich-type covariance estimators, glmmTMB's optimizer,
      ## etc.); force exact symmetry so downstream dNormal()/
      ## dIndependent_Normal_Gamma() (which validate Sigma with isSymmetric())
      ## do not spuriously reject a submatrix whose asymmetry is pure
      ## floating-point noise.
      Sigma <- (Sigma + t(Sigma)) / 2
      dimnames(Sigma) <- list(cols_k, cols_k)

      list(
        prior = list(
          mu         = mu,
          Sigma      = Sigma,
          dispersion = tau2_k
        ),
        pwt_used = pwt_k
      )
    }),
    re_names
  )
  prior_list    <- lapply(prior_raw, `[[`, "prior")
  pwt_used_list <- lapply(prior_raw, `[[`, "pwt_used")

  ## 'pop.mu'/'pop.nprior'/'pop.sd' outputs: the values actually used for
  ## every component, in the same list-of-named-vectors shape the
  ## corresponding argument expects -- always populated (derived defaults
  ## when the argument was not supplied), so a caller can inspect
  ## ps$pop.mu / ps$pop.nprior / ps$pop.sd from a first call and feed a
  ## modified copy back in as an override.
  pop_mu_out <- stats::setNames(
    lapply(re_names, function(k) prior_list[[k]]$mu),
    re_names
  )

  ## 'pop.nprior' output: effective prior sample size (group units) implied
  ## by the weight actually used for each component/column,
  ## n_i = J * w_i / (1 - w_i) -- the same duality
  ## .lmebayes_resolve_disp_prior() uses for pop.dispersion.pwt/nprior.
  pop_nprior_out <- stats::setNames(
    lapply(re_names, function(k) {
      w <- pwt_used_list[[k]]
      stats::setNames(J_groups * w / (1 - w), names(w))
    }),
    re_names
  )

  pop_sd_out <- stats::setNames(
    lapply(re_names, function(k) sqrt(diag(prior_list[[k]]$Sigma))),
    re_names
  )

  pwt_out <- if (is.numeric(pwt)) pwt else pwt_list

  block_formula <- .lmebayes_block_formula_from_re(formula, re_names)
  ## Classical VarCorr / Psi plug-ins (diagnostic); ING centers below use
  ## A12 §3.3.4 marginal RSS instead.
  pop_dispersion_ref <- stats::setNames(as.numeric(tau2_vec), re_names)
  group_tau_sd_out <- if (is_gaussian) {
    stats::setNames(sqrt(unname(tau2_vec)), re_names)
  } else {
    NULL
  }

  ## Block~2 tau^2_k center: compute_gaussian_prior on the hyper-regression
  ## b·k ~ N(W_k gamma_k, tau^2_k I_J) with prior gamma_k ~ N(mu_k, Sigma_k).
  ## Replaces the classical Psi_kk plug-in in prior_list$dispersion /
  ## pop.ing_prior while keeping pop.dispersion.ref = Psi.
  if (is_gaussian) {
    for (k in re_names) {
      Y_k <- .lmebayes_reference_b_component(
        fit_ref, design$group_name, k, group_levels
      )
      X_k <- as.matrix(design$W[[k]])
      if (!is.null(rownames(X_k)) && setequal(rownames(X_k), group_levels)) {
        X_k <- X_k[group_levels, , drop = FALSE]
      }
      n0_k <- unname(disp$n_prior_dispersion[[k]])
      cal_k <- .lmebayes_compute_ing_prior_cal_tau2_hyper(
        Y       = Y_k,
        X       = X_k,
        mu      = prior_list[[k]]$mu,
        Sigma   = prior_list[[k]]$Sigma,
        n_prior = n0_k,
        k_arg   = 1
      )
      prior_list[[k]]$dispersion <- as.numeric(cal_k$dispersion)
    }
  }

  ## Prospective dIndependent_Normal_Gamma calibration per component (used
  ## only when pfamily_list(ptypes = "dIndependent_Normal_Gamma") is chosen):
  ## Gamma precision prior shape/rate from the lmebayesCore default
  ## calibration (compute_gaussian_prior() with k = 1):
  ##   shape_ING = (n0 + 1 + p_k)/2,  b_0 = tau2_k * (n0 + p_k - 1)/2.
  ## Since b_0 = tau2_k * (shape_ING - 1), the implied inverse-Gamma prior on
  ## tau^2_k has mean exactly tau2_k for every n0 and p_k.
  ##
  ## Here tau2_k is the A12 §3.3.4 marginal hyper-regression center stored in
  ## prior_list[[k]]$dispersion (classical Psi is in pop.dispersion.ref).
  ##
  ## The tau^2 truncation window (disp_lower / disp_upper) uses the
  ## *limiting posterior* of lmebayesCore Chapter A12, Theorem 2 -- the
  ## weak-prior (n0 -> 0) limit of the Block 2 posterior Gamma:
  ##   a_inf = (J + 1)/2,  b_inf = tau2_k * (J - 1)/2
  ## (so b_inf/(a_inf - 1) = tau2_k: mean-matched, like the prior).  The
  ## window uses max_disp_perc_dispersion (default 0.99, scalar or one value
  ## per random-effect component) for the central
  ## (2*max_disp_perc_dispersion - 1) mass interval.  Quantiles of the
  ## *prior* would stretch without bound as n0 -> 0 (posterior coverage ->
  ## 100%, envelope acceptance -> 0); the limiting-posterior window instead
  ## has coverage >= (2*max_disp_perc_dispersion - 1) of the exact posterior
  ## for every n0 (the finite-n0 posterior is strictly more concentrated
  ## than the limit), is identical for all n0, and keeps the envelope
  ## sampler's candidates-per-draw roughly constant as priors weaken.  See
  ## inst/ING_TRUNCATION_WINDOW.md.  Stored here so print() can display the
  ## window and pfamily_list() consumes one source of truth.
  mdp_dispersion <- .lmebayes_expand_scalar_or_vector(
    max_disp_perc_dispersion, re_names, "pop.max_disp_perc"
  )
  ing_prior <- stats::setNames(
    lapply(re_names, function(k) {
      n0_k    <- unname(disp$n_prior_dispersion[[k]])
      p_k     <- length(prior_list[[k]]$mu)
      tau2_k  <- unname(prior_list[[k]]$dispersion)
      shape_k <- (n0_k + 1) / 2 + p_k / 2
      rate_k  <- tau2_k * (n0_k + p_k - 1) / 2
      mdp_k   <- unname(mdp_dispersion[[k]])
      win_k <- .lmebayes_ing_limiting_posterior_window(tau2_k, J_groups,
                                                       mdp_k)
      list(
        shape         = shape_k,
        rate          = rate_k,
        disp_lower    = win_k$disp_lower,
        disp_upper    = win_k$disp_upper,
        max_disp_perc = mdp_k
      )
    }),
    re_names
  )

  ing_prior_measurement <- NULL
  meas <- NULL
  mdp_meas_vector       <- length(max_disp_perc_measurement) > 1L
  mdp_measurement_pooled <- NULL
  if (is_gaussian) {
    n_obs <- length(design$y)
    pwt_meas_vector <- !is.null(pwt_measurement) && length(pwt_measurement) > 1L
    if (pwt_meas_vector && !is.numeric(pwt_measurement)) {
      stop(
        "'group.dispersion.pwt' vector must be numeric.",
        call. = FALSE
      )
    }
    meas <- .lmebayes_resolve_measurement_disp_prior(
      pwt_measurement     = if (pwt_meas_vector) NULL else pwt_measurement,
      n_prior_measurement = n_prior_measurement,
      n_obs               = n_obs
    )

    ## Mirrors pwt_measurement's own pooled-vs-vector split: a per-group
    ## max_disp_perc_measurement vector leaves the *pooled* ing_prior_measurement
    ## window at the package default (0.99) -- the vector only ever tailors
    ## the per-group ing_prior_measurement_group calibration below.
    if (mdp_meas_vector) {
      mdp_measurement_pooled <- 0.99
    } else {
      mdp_measurement_pooled <- as.numeric(max_disp_perc_measurement)
      if (length(mdp_measurement_pooled) != 1L || is.na(mdp_measurement_pooled) ||
          mdp_measurement_pooled <= 0.5 || mdp_measurement_pooled >= 1) {
        stop(
          "'group.max_disp_perc' must be a scalar in (0.5, 1), or a ",
          "length-J vector (one value per group level) for per-group ",
          "calibration.",
          call. = FALSE
        )
      }
    }

    ## Only calibrate the pooled Block~1 Gamma prior when dispformula
    ## actually requests the pooled path. When dispformula requests
    ## per-group dispersion, ing_prior_measurement_group (below) is the one
    ## actually fed to the sampler; a second, unused pooled calibration only
    ## invited the "which one is real" confusion print() used to have to
    ## paper over. meas/mdp_measurement_pooled above are still computed
    ## unconditionally -- they also back the pwt_measurement/
    ## max_disp_perc_measurement fallback values used further down whenever
    ## no vector/calibration override applies, regardless of dispformula.
    if (identical(dispformula_kind, "pooled")) {
      ## Aggregate per-group S_marg,j (A12 §3.3.4 + Part VI).
      classical_sigma2 <- as.numeric(group_dispersion)
      smarg_pool <- .lmebayes_pooled_measurement_smarg_from_groups(
        design           = design,
        data             = data,
        block_formula    = block_formula,
        sd_tau           = group_tau_sd_out,
        prior_list       = prior_list,
        group_levels     = group_levels,
        family           = family,
        intercept_source = intercept_source,
        effects_source   = effects_source
      )
      group_dispersion <- smarg_pool$sigma2_pool
      attr(group_dispersion, "source") <- "pooled S_marg (A12 3.3.4)"
      attr(group_dispersion, "classical") <- classical_sigma2
      ing_prior_measurement <- .lmebayes_calibrate_ing_prior_measurement(
        design           = design,
        group.dispersion = group_dispersion,
        n_prior          = meas$n_prior_measurement,
        max_disp_perc    = mdp_measurement_pooled
      )
    }
  }

  n_j_group <- if (is_gaussian) {
    nj <- as.integer(table(design$group))
    names(nj) <- group_levels
    nj
  } else {
    NULL
  }

  meas_group <- if (is_gaussian) {
    .lmebayes_resolve_measurement_disp_prior_group(
      pwt_measurement     = pwt_measurement,
      n_prior_measurement = NULL,
      n_j                 = n_j_group,
      group_levels        = group_levels
    )
  } else {
    NULL
  }

  mdp_measurement_group <- if (is_gaussian) {
    .lmebayes_expand_scalar_or_vector(
      max_disp_perc_measurement, group_levels, "group.max_disp_perc"
    )
  } else {
    NULL
  }

  ing_prior_measurement_group <- if (is_gaussian && identical(dispformula_kind, "group")) {
    .lmebayes_calibrate_ing_prior_measurement_group(
      design           = design,
      data             = data,
      block_formula    = block_formula,
      sd_tau           = group_tau_sd_out,
      pwt_group        = meas_group$pwt_measurement,
      n_prior_group    = meas_group$n_prior_measurement,
      group_levels     = group_levels,
      prior_list       = prior_list,
      max_disp_perc_group = mdp_measurement_group,
      family           = family,
      intercept_source = intercept_source,
      effects_source   = effects_source
    )
  } else {
    NULL
  }

  ## alpha_target_measurement: search for the smallest per-group
  ## pwt_measurement driving the predicted ellipsoid violation rate down to
  ## alpha_target_measurement (Section 16.6 exact/truncated criterion),
  ## floored at the pwt_measurement resolved above (pass-1/seed). Active by
  ## default; a non-NULL default must stay a true no-op for models that
  ## don't support per-group calibration (e.g. dispformula = ~1) -- only an
  ## *explicit* request on such a model is an error.
  supports_group_calibration <- is_gaussian && identical(dispformula_kind, "group")

  if (!missing(group.alpha_target) && !is.null(alpha_target_measurement) &&
      !supports_group_calibration) {
    stop(
      "'group.alpha_target' requires a Gaussian model with dispformula ",
      "requesting per-group dispersion (dispformula = ~<group_name>).",
      call. = FALSE
    )
  }

  pwt_measurement_calibration <- NULL
  if (!is.null(alpha_target_measurement) && supports_group_calibration) {
    if (!is.numeric(alpha_target_measurement) ||
        length(alpha_target_measurement) != 1L ||
        is.na(alpha_target_measurement) ||
        alpha_target_measurement <= 0 || alpha_target_measurement >= 1) {
      stop("'group.alpha_target' must be a scalar in (0, 1).", call. = FALSE)
    }
    calib <- .lmebayes_calibrate_pwt_measurement_group(
      fit_ref                     = fit_ref,
      design                      = design,
      group_levels                = group_levels,
      groupef.names               = re_names,
      group.Sigma                 = group.Sigma,
      ing_prior_measurement_group = ing_prior_measurement_group,
      alpha_target                = alpha_target_measurement,
      floor_vec                   = meas_group$pwt_measurement
    )
    meas_group <- .lmebayes_resolve_measurement_disp_prior_group(
      pwt_measurement     = calib$pwt_measurement,
      n_prior_measurement = NULL,
      n_j                 = n_j_group,
      group_levels        = group_levels
    )
    meas_group$source <- sprintf(
      "group.alpha_target-calibrated (target = %.4g, floored at group.dispersion.pwt)",
      alpha_target_measurement
    )
    ing_prior_measurement_group <- .lmebayes_calibrate_ing_prior_measurement_group(
      design           = design,
      data             = data,
      block_formula    = block_formula,
      sd_tau           = group_tau_sd_out,
      pwt_group        = meas_group$pwt_measurement,
      n_prior_group    = meas_group$n_prior_measurement,
      group_levels     = group_levels,
      prior_list       = prior_list,
      max_disp_perc_group = mdp_measurement_group,
      family           = family,
      intercept_source = intercept_source,
      effects_source   = effects_source
    )
    pwt_measurement_calibration <- calib$table
  }

  ## group.dispersion.pwt / group.dispersion.nprior share one "which representation is live"
  ## switch: per-group vector/values once a group.dispersion.pwt vector was supplied or
  ## alpha_target calibration ran, pooled (n_obs-scale) scalar otherwise --
  ## so the two fields always describe the same pooled-or-group state
  ## together (dispformula alone is not the switch: e.g. dispformula =
  ## ~group with a scalar group.dispersion.pwt and calibration disabled still reports
  ## the pooled scalar here).
  meas_group_active <- is_gaussian &&
    ((!is.null(pwt_measurement) && length(pwt_measurement) > 1L) ||
       !is.null(pwt_measurement_calibration))
  pwt_measurement_out <- if (is_gaussian) {
    if (meas_group_active) {
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
    if (meas_group_active) {
      np <- meas_group$n_prior_measurement
      attr(np, "source") <- meas_group$source
      np
    } else {
      np <- meas$n_prior_measurement
      attr(np, "source") <- meas$source
      np
    }
  } else {
    NULL
  }
  max_disp_perc_measurement_out <- if (is_gaussian) {
    if (mdp_meas_vector) mdp_measurement_group else mdp_measurement_pooled
  } else {
    NULL
  }

  if (is_gaussian && identical(dispformula_kind, "group") &&
      !is.null(ing_prior_measurement_group)) {
    group_dispersion <- vapply(
      group_levels,
      function(lev) {
        as.numeric(ing_prior_measurement_group[[lev]]$sigma2_hat)
      },
      numeric(1L)
    )
    names(group_dispersion) <- group_levels
    attr(group_dispersion, "source") <- "per-group S_marg (A12 3.3.4)"
  }

  structure(
    list(
      formula                = formula,
      family                 = family,
      pop.pwt                  = pwt_out,
      pop.nprior               = pop_nprior_out,
      pop.sd                   = pop_sd_out,
      pop.mu                   = pop_mu_out,
      group.dispersion         = group_dispersion,
      pop.dispersion.pwt       = disp$pwt_dispersion,
      pop.dispersion.nprior    = disp$n_prior_dispersion,
      group.dispersion.pwt     = pwt_measurement_out,
      group.dispersion.nprior  = n_prior_measurement_out,
      group.alpha_target       = if (supports_group_calibration) alpha_target_measurement else NULL,
      group.pwt_calibration    = pwt_measurement_calibration,
      group.max_disp_perc      = max_disp_perc_measurement_out,
      pop.max_disp_perc        = mdp_dispersion,
      dispformula              = dispformula,
      pop.intercept_source     = intercept_source,
      pop.effects_source       = effects_source,
      data                     = data,
      block_formula            = block_formula,
      group.tau_sd             = group_tau_sd_out,
      design                   = design,
      mer_fit                  = mer_fit,
      fit_ref                  = fit_ref,
      group.dispersion.fit     = if (identical(calibration_source, "glmmTMB")) fit_ref else NULL,
      group.dispersion.ref     = group_dispersion_ref,
      pop.dispersion.ref       = pop_dispersion_ref,
      calibration_source       = calibration_source,
      group.Sigma              = group.Sigma,
      pop.prior_list       = prior_list,
      pop.ing_prior        = ing_prior,
      group.ing_prior    = if (!is.null(ing_prior_measurement)) {
        ing_prior_measurement
      } else {
        ing_prior_measurement_group
      }
    ),
    class = "Prior_Setup_GLMM"
  )
}

## Structural check distinguishing the two shapes 'group.ing_prior' can take:
## pooled (dispformula = ~1) is a single flat list of scalars (shape, rate,
## disp_lower, ...); per-group (dispformula = ~<group_name>) is a named list
## with one such list per group level, so its *first* element is itself a
## list. NULL/empty input is treated as "not grouped" (pooled shape, or no
## calibration at all).
#' @keywords internal
#' @noRd
.lmebayes_ing_prior_is_grouped <- function(x) {
  !is.null(x) && length(x) > 0L && is.list(x[[1L]])
}

## Dimension-adaptive default Block~2 prior weights: one population model
## per random-effect component k, with p_k = ncol(W[[k]]).  Mirrors
## glmbayesCore::Prior_Setup()'s nvar < 14 vs >= 14 split, but applied
## separately to each RE component rather than to the full model.
#' @keywords internal
#' @noRd
.lmebayes_default_pwt_list <- function(design,
                                       pwt_default_low,
                                       pwt_default_high,
                                       threshold = 14L) {
  re_names <- design$groupef.names
  stats::setNames(
    lapply(re_names, function(k) {
      cols_k <- colnames(design$W[[k]])
      p_k    <- length(cols_k)
      w <- if (p_k >= threshold) pwt_default_high else pwt_default_low
      stats::setNames(rep(as.numeric(w), p_k), cols_k)
    }),
    re_names
  )
}

## Resolve 'pwt' (scalar or list) into a canonical named list with one named
## numeric vector per random-effect component, ordered like X_hyper[[k]].
#' @keywords internal
#' @noRd
.lmebayes_resolve_pwt <- function(pwt, design) {

  re_names <- design$groupef.names
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
        cols_k <- colnames(design$W[[k]])
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
    cols_k <- colnames(design$W[[k]])
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

## Resolve an optional 'pop.mu'/'pop.sd' override (NULL, scalar, or a list
## with one element per random-effect component, entries of which may be
## NULL for a partial override) into a canonical named list with one
## element per random-effect component: either NULL (no override; caller
## keeps its usual default for that component) or a full named numeric
## vector aligned to colnames(design$W[[k]]).  Shares .lmebayes_resolve_pwt()'s
## list-shape validation, but (a) allows per-component NULL entries and
## (b) has no (0, 1) range constraint -- 'check_positive' instead requires
## strictly positive values for 'pop.sd' while leaving 'pop.mu' unconstrained.
#' @keywords internal
#' @noRd
.lmebayes_resolve_override_list <- function(x, design, what, check_positive) {

  re_names <- design$groupef.names
  p_re     <- length(re_names)

  if (is.null(x)) {
    return(stats::setNames(vector("list", p_re), re_names))
  }

  check_values <- function(v, lbl) {
    if (!is.numeric(v) || anyNA(v)) {
      stop(sprintf("%s must be numeric without missing values.", lbl),
           call. = FALSE)
    }
    if (check_positive && any(v <= 0)) {
      stop(sprintf("%s must be strictly positive.", lbl), call. = FALSE)
    }
  }

  if (is.numeric(x)) {
    if (length(x) != 1L) {
      stop(sprintf(
        paste0("'%s' must be a scalar, NULL, or a list with one element ",
               "per random-effect component."),
        what
      ), call. = FALSE)
    }
    check_values(x, sprintf("'%s'", what))
    return(stats::setNames(
      lapply(re_names, function(k) {
        cols_k <- colnames(design$W[[k]])
        stats::setNames(rep(as.numeric(x), length(cols_k)), cols_k)
      }),
      re_names
    ))
  }

  if (!is.list(x)) {
    stop(sprintf(
      paste0("'%s' must be a scalar, NULL, or a list with one element per ",
             "random-effect component."),
      what
    ), call. = FALSE)
  }

  if (length(x) != p_re) {
    stop(sprintf(
      paste0("'%s' list has length %d but there are %d random-effect ",
             "components (%s)."),
      what, length(x), p_re, paste(re_names, collapse = ", ")
    ), call. = FALSE)
  }

  nms <- names(x)
  if (!is.null(nms) && any(nzchar(nms))) {
    if (!setequal(nms, re_names)) {
      stop(sprintf(
        "Names of '%s' must match the random-effect coefficient names: %s.",
        what, paste(re_names, collapse = ", ")
      ), call. = FALSE)
    }
    x <- x[re_names]
  } else {
    names(x) <- re_names
  }

  out <- stats::setNames(vector("list", p_re), re_names)
  for (k in re_names) {
    v <- x[[k]]
    if (is.null(v)) next  # partial override: keep the usual default for 'k'

    cols_k <- colnames(design$W[[k]])
    p_k    <- length(cols_k)
    lbl    <- sprintf("'%s[[\"%s\"]]'", what, k)
    check_values(v, lbl)

    if (length(v) == 1L) {
      v <- rep(as.numeric(v), p_k)
    } else if (length(v) == p_k) {
      vn <- names(v)
      if (!is.null(vn) && any(nzchar(vn))) {
        if (!setequal(vn, cols_k)) {
          stop(sprintf(
            "Names of %s must match the Block 2 predictors: %s.",
            lbl, paste(cols_k, collapse = ", ")
          ), call. = FALSE)
        }
        v <- v[cols_k]
      }
      v <- as.numeric(v)
    } else {
      stop(sprintf(
        "%s must have length 1 or %d (one value per Block 2 predictor).",
        lbl, p_k
      ), call. = FALSE)
    }
    out[[k]] <- stats::setNames(v, cols_k)
  }
  out
}

## Prior-mean override: unconstrained reals.
#' @keywords internal
#' @noRd
.lmebayes_resolve_mu_override <- function(mu, design) {
  .lmebayes_resolve_override_list(mu, design, "pop.mu", check_positive = FALSE)
}

## Prior-sd override: strictly positive (coefficient-scale standard
## deviations), converted to an effective weight by the caller.
#' @keywords internal
#' @noRd
.lmebayes_resolve_sd_override <- function(sd, design) {
  .lmebayes_resolve_override_list(sd, design, "pop.sd", check_positive = TRUE)
}

## Prior-nprior override: strictly positive effective prior sample size(s)
## in group units (n_k = J * w_k / (1 - w_k)), converted to an effective
## weight w_i = n_i / (n_i + J) by the caller (mirrors
## .lmebayes_resolve_disp_prior()'s pwt/nprior duality, but per-column and
## with NULL-partial-override support).
#' @keywords internal
#' @noRd
.lmebayes_resolve_nprior_override <- function(nprior, design) {
  .lmebayes_resolve_override_list(nprior, design, "pop.nprior", check_positive = TRUE)
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

#' @rdname Prior_Setup_GLMM
#' @method print Prior_Setup_GLMM
#' @param x Object of class \code{"Prior_Setup_GLMM"}.
#' @param digits Number of decimal places for numeric output. Default 4.
#' @param ... Ignored.
#' @return \code{x} invisibly.
#' @export
print.Prior_Setup_GLMM <- function(x, digits = 4L, ...) {

  re_names <- x$design$groupef.names
  n_fr     <- sum(x$design$groupef.rank)
  n_all    <- nlevels(x$design$group)

  disp_src <- attr(x$pop.dispersion.pwt, "source")

  cat("Call: Prior_Setup_GLMM()\n")

  ## ---- Setup ---------------------------------------------------------------
  cat("\n--- Setup ---\n")
  cat("  Model family, prior-mean sources, and reference fit used for calibration.\n\n")
  cat(sprintf("  family                : %s (%s link)\n",
              x$family$family, x$family$link))
  cat(sprintf("  pop.intercept_source  : %s\n",
              if (!is.null(x$pop.intercept_source)) x$pop.intercept_source else "full_model"))
  cat(sprintf("  pop.effects_source    : %s\n",
              if (!is.null(x$pop.effects_source)) x$pop.effects_source else "full_model"))
  cat(sprintf("  dispformula           : %s\n",
              if (!is.null(x$dispformula)) deparse(x$dispformula) else "~1"))
  cat(sprintf("  calibration_source    : %s  (fixef / tau^2_k / sd_tau)\n",
              if (!is.null(x$calibration_source)) x$calibration_source else "lme4"))

  ## ---- Group dispersion ----------------------------------------------------
  cat("\n--- Group dispersion (sigma^2) ---\n")
  if (!is.null(x$group.ing_prior) &&
      .lmebayes_ing_prior_is_grouped(x$group.ing_prior)) {
    cat("  Per-group Gamma prior on sigma^2 is stored in group.ing_prior.\n")
    cat("  Inspect with print(dGamma_list(<Prior_Setup_GLMM object>)).\n")
  } else if (!is.null(x$group.ing_prior)) {
    cat("  Pooled Gamma prior on sigma^2 is stored in group.ing_prior.\n")
    cat("  (dGamma_list() requires a per-group dispformula.)\n")
    if (!is.null(x$group.dispersion)) {
      cat(sprintf(
        "  group.dispersion       : %.4g  [%s]\n",
        as.numeric(x$group.dispersion),
        if (is.null(attr(x$group.dispersion, "source"))) {
          "unknown"
        } else {
          attr(x$group.dispersion, "source")
        }
      ))
      classical <- attr(x$group.dispersion, "classical")
      if (!is.null(classical)) {
        cat(sprintf(
          "  classical lmer sigma^2 : %.4g  (diagnostic)\n",
          as.numeric(classical)
        ))
      }
    }
  } else {
    cat("  No group-dispersion Gamma prior (non-Gaussian or not calibrated).\n")
  }

  ## ---- Design check --------------------------------------------------------
  cat("\n--- Design check (from model_setup) ---\n")
  cat("  Identifiability flags copied from the design object (not re-checked here).\n\n")
  cat(sprintf(
    "  Full-rank groups (algebraic D_j): %d of %d %s\n",
    n_fr, n_all, x$design$group_name
  ))
  if (!is.null(x$design$groupef.glm_check) &&
      !is.null(x$design$groupef.estimable)) {
    n_est <- sum(x$design$groupef.estimable[x$design$groupef.rank])
    cat(sprintf(
      paste0(
        "  Full-rank & estimable           : %d of %d full-rank ",
        "(%d of %d total %s)\n"
      ),
      n_est, n_fr, n_est, n_all, x$design$group_name
    ))
  }

  ## ---- Group-effect covariance ---------------------------------------------
  cat("\n--- Group-effect covariance (group.Sigma) ---\n")
  cat("  Diagonal Psi with entries tau^2_k for each group-effect coefficient.\n\n")
  print(round(x$group.Sigma, digits))

  ## ---- Population priors ---------------------------------------------------
  cat("\n--- Population priors (pop.prior_list) ---\n")
  cat("  One prior list per group-effect coefficient (input to pfamily_list()).\n")
  for (nm in re_names) {
    pl <- x$pop.prior_list[[nm]]
    cat(sprintf("\n  [%s]\n", nm))
    pwt_k <- if (is.numeric(x$pop.pwt)) x$pop.pwt else x$pop.pwt[[nm]]
    pwt_str <- if (length(unique(pwt_k)) == 1L) {
      sprintf("%.4g", pwt_k[1L])
    } else {
      paste(sprintf("%s=%.4g", names(pwt_k), pwt_k), collapse = ", ")
    }
    cat(sprintf("  pop.pwt               : %s\n", pwt_str))
    if (!is.null(x$pop.dispersion.pwt)) {
      cat(sprintf(
        "  pop.dispersion.pwt     : %.4g  [%s]\n",
        x$pop.dispersion.pwt[[nm]],
        if (is.null(disp_src)) "unknown" else disp_src
      ))
    }
    if (!is.null(x$pop.dispersion.nprior)) {
      cat(sprintf(
        "  pop.dispersion.nprior  : %.4g  (= J * pwt_disp / (1 - pwt_disp))\n",
        x$pop.dispersion.nprior[[nm]]
      ))
    }
    cat("\n  mu:\n")
    print(round(pl$mu, digits))
    cat("\n  Sigma:\n")
    print(round(pl$Sigma, digits))
    cat("\n")
    cat(sprintf(
      "  dispersion             : %.4f  (tau^2_k; A12 3.3.4 marginal)\n",
      pl$dispersion))
    if (!is.null(x$pop.dispersion.ref) && nm %in% names(x$pop.dispersion.ref)) {
      cat(sprintf(
        "  pop.dispersion.ref     : %.4f  (classical VarCorr / Psi)\n",
        unname(x$pop.dispersion.ref[[nm]])
      ))
    }
    ing_k <- x$pop.ing_prior[[nm]]
    if (!is.null(ing_k)) {
      mdp_k <- if (!is.null(ing_k$max_disp_perc)) ing_k$max_disp_perc else 0.99
      cat(sprintf(
        paste0(
          "  ING tau^2 window       : [%.4g, %.4g]  ",
          "(%.4g/%.4g limiting-posterior quantiles; upper/tau2 = %.3g)\n"
        ),
        ing_k$disp_lower, ing_k$disp_upper,
        1 - mdp_k, mdp_k,
        ing_k$disp_upper / unname(pl$dispersion)
      ))
      cat(sprintf(
        paste0(
          "  ING shape, rate        : %.4g, %.4g  ",
          "(Gamma on 1/tau^2_k; for ptypes = \"dIndependent_Normal_Gamma\")\n"
        ),
        ing_k$shape, ing_k$rate
      ))
    }
  }

  invisible(x)
}
