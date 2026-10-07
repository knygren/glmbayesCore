#' Simulation Functions for Linear Mixed Models
#'
#' @description
#' Posterior draws for Bayesian linear mixed models, given design matrices
#' \code{y}, \code{D}, \code{W} and population priors \code{pfamily_list}.
#' The model is stated under \dQuote{Model and notation} below.
#'
#' These are the engine-level entry points, normally reached through
#' \code{\link{rlmerb}} (or \code{lmerb()} in \pkg{lmebayes}), which is also
#' where the returned object is documented field by field. Each returned
#' draw is one independent replicate chain's state at its final sweep,
#' except on the fully-known-variance route, where draws are exact.
#'
#' @section Model and notation:
#' For group \eqn{j = 1, \ldots, J} (\eqn{J = }\code{length(group_levels)}),
#' with \eqn{n_j} observations in group \eqn{j} and \eqn{p_{re} = }
#' \code{ncol(D)} group-varying coefficients:
#'
#' \strong{Likelihood (group level):}
#' \deqn{y_j \mid \beta_j, \sigma_j^2 \ \sim\ N\bigl(D_j\beta_j,\ \sigma_j^2 I_{n_j}\bigr)}
#' \strong{Hierarchical prior (population level):}
#' \deqn{\beta_j \mid \gamma, \Psi \ \sim\ N\bigl(\mathcal{W}_j\gamma,\ \Psi\bigr)}
#'
#' The full, non-centered coefficient vector \eqn{\beta_j} -- not a mean-zero
#' deviation -- appears directly in the likelihood. This is a \emph{centered
#' parameterization} in the sense of \insertCite{LindleySmith1972}{glmbayesCore}:
#' \eqn{\beta_j} is the group parameter, \eqn{\gamma} the population
#' (hyper-mean) parameter, and \eqn{\Psi} the population covariance.
#' \eqn{\mathcal{W}_j} is block-diagonal across the \eqn{p_{re}}
#' coefficient dimensions,
#' \deqn{\mathcal{W}_j = \mathrm{blockdiag}\bigl(W_{1j}, \ldots, W_{p_{re}j}\bigr),
#'   \qquad W_{kj} \in \mathbb{R}^{1 \times q_k},}
#' and each coefficient \eqn{k} may have its own number of level-2 predictors
#' \eqn{q_k} (no padding needed, unlike a 3-index array).
#'
#' \strong{Correspondence to these functions' arguments/return values:}
#' \describe{
#'   \item{\eqn{D} (row-stacked \eqn{D_j})}{\code{D}: the \eqn{l_2 \times p_{re}}
#'     level-1 design matrix, one row per observation; \eqn{D_j} is the
#'     submatrix \code{D[group == } \eqn{j}\code{, ]}.}
#'   \item{\eqn{\mathcal{W}_j}}{Assembled on demand from \code{W}: block
#'     \eqn{W_{kj}} is row \eqn{j} of \code{W[[k]]}, so \code{W} \emph{is}
#'     \code{list(}\eqn{W_1, \ldots, W_{p_{re}}}\code{)}, one
#'     \eqn{J \times q_k} matrix per RE column \eqn{k}.}
#'   \item{\eqn{\gamma}}{\code{popef}/\code{popef.mode} in the returned object
#'     (one entry per RE column, named by \code{colnames(W[[k]])}).}
#'   \item{\eqn{\beta_j}}{The group-\eqn{j} row of \code{groupef}/
#'     \code{groupef.mode} -- \strong{not} \code{lme4}'s mean-zero \eqn{b_j}
#'     (see below).}
#'   \item{\eqn{\Psi}}{Diagonal, \eqn{\Psi = \mathrm{diag}(\tau^2_1, \ldots,
#'     \tau^2_{p_{re}})}, one \eqn{\tau^2_k} plug-in per \code{pfamily_list}
#'     component (see \code{pfamily_list} below); \eqn{\Psi^{-1}} is the
#'     internally-derived group prior precision \code{P}.}
#'   \item{group index \eqn{j}}{\code{group} (\code{levels(group)} fixes the
#'     order \eqn{j = 1, \ldots, J}).}
#' }
#'
#' \strong{Correspondence to \code{lme4}/Laird-Ware:} marginalizing (substituting
#' the population equation into the group likelihood) recovers the standard
#' \code{lme4} form
#' \deqn{y_j = D_j(\mathcal{W}_j\gamma + u_j) + \varepsilon_j =
#'   \underbrace{(D_j\mathcal{W}_j)}_{X_j}\gamma + \underbrace{D_j}_{Z_j}u_j + \varepsilon_j,}
#' where \eqn{u_j := \beta_j - \mathcal{W}_j\gamma} is the mean-zero deviation
#' \code{lme4} calls \eqn{b_j}. \eqn{D_j} itself has no single \code{lme4}
#' counterpart: it is the common ancestor of both \eqn{X_j} (fixed-effects
#' design, after right-multiplying by \eqn{\mathcal{W}_j}) and \eqn{Z_j}
#' (random-effects design, used as-is). \eqn{\gamma} corresponds to
#' \code{lme4}'s \eqn{\beta} (fixed effects), and \eqn{\Psi} to the per-group
#' block of \code{lme4}'s \eqn{G} matrix (\eqn{G = I_J \otimes \Psi}); the
#' functions documented here never materialize \eqn{u_j}/\eqn{b_j} -- \eqn{\beta_j}
#' is sampled directly against \eqn{y_j} at every group-level Gibbs update (the
#' per-group/shared-\eqn{\gamma} conjugate updates for \eqn{\beta_j} and
#' \eqn{\gamma} are the \eqn{Z_j}/\eqn{b_j}-lettered formulas in
#' \code{\link{lmebayes_posterior_icm}}, with \eqn{Z_j \equiv D_j} and
#' \eqn{b_j \equiv \beta_j} there). When every \eqn{\mathcal{W}_j} is
#' intercept-only (\eqn{q_k = 1} for all \eqn{k}, i.e. no level-2 predictors),
#' the model reduces to a standard random-intercept/random-slope \code{lmer}
#' model: \eqn{\hat\gamma} should then recover \code{lme4}'s fixed effects and
#' \eqn{\hat\beta_j} should track \code{lme4::coef()} (\code{lme4::ranef()}
#' tracks \eqn{\hat\beta_j - \hat\gamma}), giving a natural automated
#' equivalence check.
#'
#' For the full write-up (including the block-diagonal Gibbs full-conditional
#' derivation) see \code{notation.md} in \code{system.file(package =
#' "lmebayesCore")} (source: \code{inst/notation.md}).
#'
#' @section What the prior choices mean:
#' Callers never name a route. It follows from two modelling decisions, and
#' the decisions --- not the routing --- are what matter statistically.
#'
#' \strong{1. Is the residual variance \eqn{\sigma^2} known or estimated?}
#' Supplying a positive scalar (or one value per group) treats it as
#' \strong{known} and conditions on it. Supplying a
#' \code{\link[glmbayesCore]{dGamma}()} prior --- one pooled, or a
#' \code{\link{dGamma_list}} with one per group --- \strong{estimates} it.
#'
#' \strong{2. Are the between-group variances \eqn{\tau^2_k} known or
#' estimated?} A \code{dNormal} component treats \eqn{\tau^2_k} as
#' \strong{known}; a single \code{dIndependent_Normal_Gamma} component
#' anywhere in \code{pfamily_list} makes that variance a parameter to be
#' \strong{estimated}.
#'
#' The cross of the two selects one of four engines:
#' \tabular{lll}{
#'   \strong{\eqn{\sigma^2}} \tab \strong{\eqn{\tau^2_k}} \tab \strong{Engine} \cr
#'   known \tab all known \tab \code{rLMMNormal_reg_known_vcov} \cr
#'   known \tab any estimated \tab \code{rLMMNormal_reg_estimated_vcov} \cr
#'   estimated \tab all known \tab \code{rLMMindepNormalGamma_reg_known_vcov} \cr
#'   estimated \tab any estimated \tab \code{rLMMindepNormalGamma_reg_estimated_vcov} \cr
#' }
#'
#' The consequence that shows up in the output is exactness. When both
#' variances are known the posterior is exactly multivariate normal, so
#' draws are exact, \code{popef.mode} is the true posterior \emph{mean}, and
#' \code{m_convergence} is \code{1}. As soon as either variance is
#' estimated the posterior is no longer Gaussian: draws come from the final
#' sweep of a calibrated chain, \code{popef.mode} is a conditional mode
#' rather than a mean, and a short pilot stage runs first.
#'
#' \code{sim_method} is therefore meaningful only in the fully-known case,
#' the one route offering a choice between drawing directly from the closed
#' form (\code{"DEFAULT"}) and running the sweep loop on the same target
#' (\code{"TWO_BLOCK_GIBBS"}, which exists mainly to cross-check the two
#' against each other). Elsewhere the argument is ignored.
#'
#' One input you cannot supply: the group prior precision
#' \eqn{\Psi^{-1}}. It is always derived from \code{pfamily_list}, and a
#' caller-supplied \code{P}/\code{Sigma} is rejected, so it cannot silently
#' disagree with the population priors. See
#' \code{vignette("Chapter-B01", package = "lmebayesCore")} for the engine
#' map and \code{vignette("Chapter-B04", package = "lmebayesCore")} for the
#' exact route.
#'
#' @section Dispatchers:
#' \code{\link{rLMMNormal_reg}} and \code{\link{rLMMindepNormalGamma_reg}}
#' validate inputs and delegate to the appropriate route (or, for the legacy
#' outer-loop \code{rLMMindepNormalGamma_reg}, draw \eqn{\sigma^2} then call
#' \code{rLMMNormal_reg}).
#'
#' @param n Number of stored draws. If \code{length(n) > 1}, the length is used.
#' @param y Response vector of length \code{l2} (\code{= nrow(D)}).
#' @param D Level-1 design matrix (\code{l2 x p_re}); this is \eqn{D}
#'   (row-stacked \eqn{D_j}) in \dQuote{Model and notation} below. Must have
#'   unique, non-empty \code{colnames(D)}: these are the random-effect
#'   coefficient names used to key \code{W} and \code{pfamily_list}
#'   (there is no separate \code{groupef.names} argument to override them).
#' @param group Grouping factor of length \code{l2} (must be a \code{factor};
#'   \code{levels(group)} fixes the row order of \code{groupef} draws -- there
#'   is no separate \code{group_levels} argument. To use a level order/superset
#'   not present in the observed data, construct \code{group} as
#'   \code{factor(observed_group, levels = full_superset)} yourself. The name
#'   used for the grouping column in \code{groupef} (\code{group_name})
#'   is resolved from \code{attr(group, "group_name")} if set, otherwise from
#'   \code{group}'s own variable name via \code{substitute()} -- this only
#'   works when \code{group} is passed as a bare variable (e.g.
#'   \code{group = school_id}); otherwise attach the name yourself via
#'   \code{attr(group, "group_name") <- "school_id"}.
#' @param W Named list of group-level design matrices (\code{J x q_k}),
#'   one per column of \code{D}; this is \eqn{\mathcal{W}} in \dQuote{Model
#'   and notation} below: \code{W[[k]]} is \eqn{W_k} and \code{W[[k]][j, ]}
#'   is the block \eqn{W_{kj}} of \eqn{\mathcal{W}_j}.
#' @param pfamily_list Named list of population \code{pfamily} objects. The
#'   group random-effect prior precision (formerly a separate \code{P}
#'   argument) is always derived internally from \code{pfamily_list}: one
#'   \eqn{\tau^2_k} plug-in per component (fixed \code{dispersion} for
#'   \code{dNormal}, prior mean \eqn{rate/(shape - 1)} for
#'   \code{dIndependent_Normal_Gamma}), assembled into a diagonal precision
#'   matrix -- this is \eqn{\Psi^{-1}} in \dQuote{Model and notation} below.
#'   There is no way to supply a precision inconsistent with
#'   \code{pfamily_list}.
#' @param dispprior_list Observation-dispersion prior:
#'   \code{list(dispersion = sigma2)} for fixed \eqn{\sigma^2} routes
#'   (\code{sigma2} a single positive scalar, pooled across groups, or -- for
#'   \code{rLMMNormal_reg}/\code{rLMMNormal_reg_known_vcov}/
#'   \code{rLMMNormal_reg_estimated_vcov} only -- a numeric vector of length
#'   \code{nlevels(group)} giving one fixed, known dispersion per group,
#'   matched to \code{levels(group)} either positionally or by name);
#'   \code{dGamma()} fields for legacy \code{rLMMindepNormalGamma_reg}; or,
#'   for \code{rLMMindepNormalGamma_reg_known_vcov}/
#'   \code{rLMMindepNormalGamma_reg_estimated_vcov}, one of three shapes for
#'   the ING observation-dispersion prior: (1) a single
#'   \code{\link[glmbayesCore]{dGamma}()} pfamily (pooled \eqn{\sigma^2}
#'   shared across groups) or (2) a named list of \code{dGamma()} pfamilies,
#'   one per group level, as returned by \code{\link{dGamma_list}()}
#'   (per-group \eqn{\sigma^2_j}) -- both preferred, since
#'   \code{shape}/\code{rate}/\code{disp_lower}/\code{disp_upper} are read
#'   straight from each pfamily's own \code{prior_list} -- or (3) the legacy
#'   flat list itself (\code{shape}/\code{rate}/\code{disp_upper}/
#'   \code{disp_lower}, or \code{shape_group}/\code{rate_group}/
#'   \code{disp_lower_group}/\code{disp_upper_group}). In every case,
#'   group-level \code{mu} and \code{Sigma} for the random-effect coefficients
#'   are \emph{not} read from \code{dispprior_list}: \code{mu} is always
#'   \eqn{\mathcal{W}_j\gamma} (recomputed every sweep), and \code{Sigma} is
#'   always derived internally as \code{solve(P)} from \code{pfamily_list}.
#' @param offset,weights Observation offset and prior weights (glmbayes-style
#'   formals: \code{offset = NULL}, \code{weights = 1}). Normalized to length
#'   \code{length(y)} and echoed on the return as \code{offset}/
#'   \code{offset2}/\code{prior.weights}. \strong{Not yet used} by the
#'   mixed-model sampling path (ICM / sweeps still assume unit weights and
#'   zero offset).
#' @param icm_tol,icm_maxit ICM convergence controls for the internal population start.
#' @param tv_tol Total-variation tolerance in \code{(0, 1)} for calibration.
#'   Inner Gibbs sweeps per stored draw (\code{m_convergence}) are derived from
#'   Theorem~3 at the ICM population start; pilot chain counts likewise.
#' @param progbar Show a text progress bar during sampling.
#' @param verbose Print convergence calibration / ICM lines.
#' @param gap_tol,mode_gap_max,diag_sweeps,stage_verbose Pilot-stage controls for
#'   \code{rLMMNormal_reg_estimated_vcov} and ING estimated routes (see route docs).
#'
#' @return An object of class \code{c("<route>", "rLMMNormal_reg", "list")}
#'   (or \code{c("rLMMindepNormalGamma_reg", "list")} for the legacy outer-loop
#'   engine), where \code{<route>} is the specific export that ran (e.g.
#'   \code{"rLMMNormal_reg_known_vcov"}). Components use package
#'   \strong{group}/\strong{population} names (see \file{inst/notation.md}),
#'   in glm/glmbayes-style order:
#'   \describe{
#'     \item{\code{groupef}}{Draws of non-centered group coefficients
#'       \eqn{\beta_j}: a data frame with the grouping column plus one column
#'       per \code{colnames(D)}. Each row is one group in one stored draw.}
#'     \item{\code{groupef.mode}}{\eqn{J \times p_{re}} matrix of ICM (or
#'       exact posterior-mean) group coefficients \eqn{\hat\beta_j}; rows
#'       are \code{levels(group)}, columns are \code{colnames(D)}.
#'       \strong{Not} \code{lme4}'s mean-zero \eqn{u_j}.}
#'     \item{\code{groupef.iters}}{Optional group-level envelope iteration counts.}
#'     \item{\code{popef}}{Named list of \code{n x q_k} matrices of population
#'       coefficient draws \eqn{\gamma_k}.}
#'     \item{\code{popef.mode}, \code{popef.init}}{ICM (or exact) population
#'       point estimates and main-stage starts.}
#'     \item{\code{popef.dispersion}, \code{popef.iters}}{Optional population
#'       per-draw diagnostics when produced by the sampler.}
#'     \item{\code{group.dispersion}}{Observation residual variance
#'       \eqn{\sigma^2}/\eqn{\sigma^2_j} when present (fixed scalar/vector or
#'       draws under a Gamma measurement prior). Optional
#'       \code{group.dispersion.mean} / \code{group.dispersion.iters}.
#'       Distinct from \code{popef.dispersion} (\eqn{\tau^2_k}).}
#'     \item{\code{pfamily_list}, \code{dispprior_list}}{Population
#'       (\code{pfamily_list}) and observation-dispersion (\code{dispprior_list})
#'       priors that were used.}
#'     \item{\code{prior.weights}, \code{offset}, \code{offset2}}{Normalized
#'       copies of the \code{weights}/\code{offset} arguments (glmbayes
#'       naming). Not yet consumed by sampling.}
#'     \item{\code{any_non_normal}}{Whether any population component is not
#'       \code{dNormal}.}
#'     \item{\code{family}, \code{design}, \code{n}}{Likelihood family, matrix
#'       inputs (\code{y}, \code{D}, \code{group}, \code{W},
#'       \code{groupef.names}, \code{group_name}, plus echoed
#'       \code{weights}/\code{offset}), and chain count.}
#'     \item{\code{call}}{Matched call.}
#'     \item{\code{m_convergence}}{Inner Gibbs sweeps per stored draw
#'       (always \code{1L} for the exact-iid known-vcov route).}
#'     \item{\code{convergence_info}}{Calibration details, including
#'       \code{draw_engine}, \code{sim_method_used}, \code{icm_info}, and
#'       optional pilot UB fields.}
#'     \item{\code{pilot}}{When a pilot ran: list with \code{n},
#'       \code{m_convergence}, \code{chisq}, and \code{draws}; otherwise
#'       \code{NULL}.}
#'     \item{\code{sweep_history}}{Main-stage sweep history when collected.}
#'   }
#'
#' @references
#' \insertAllCited{}
#' @importFrom Rdpack reprompt
#' @family simfuncs
#' @seealso \code{rGLMM_reg()} in \pkg{lmebayesCore}, \code{\link{rlmerb}}, \code{\link{print_groupef}},
#'   \code{\link[glmbayesCore]{rindepNormalGamma_reg}},
#'   \code{\link{lmebayes_posterior_icm}}
#' @example inst/examples/Ex_rLMM_reg.R
#' @name rLMM_reg
#' @order 1
NULL

#' Shared matrix-level validation for LMM replicate-chain engines
#'
#' \code{groupef.names} and \code{group_levels} are no longer separate
#' arguments: they are always \code{colnames(D)} and \code{levels(group)}
#' respectively. \code{group_name} must already be resolved by the caller
#' (see \code{\link{.lmebayes_resolve_group_name}}); this function only
#' sanity-checks it.
#' @noRd
.rLMM_validate_matrix_inputs <- function(
    n,
    y,
    D,
    W,
    tv_tol,
    group_name,
    group
) {
  if (length(n) > 1L) n <- length(n)
  n <- as.integer(n[1L])
  if (n < 1L) stop("'n' must be at least 1.", call. = FALSE)

  if (!is.numeric(tv_tol) || length(tv_tol) != 1L ||
      !is.finite(tv_tol) || tv_tol <= 0 || tv_tol >= 1) {
    stop("'tv_tol' must be a single value in (0, 1).", call. = FALSE)
  }

  y <- as.vector(y)
  D <- as.matrix(D)
  l2 <- nrow(D)
  if (length(y) != l2) {
    stop("length(y) must equal nrow(D).", call. = FALSE)
  }

  re_names <- colnames(D)
  if (is.null(re_names) || length(re_names) != ncol(D) || anyNA(re_names) ||
      any(!nzchar(re_names)) || anyDuplicated(re_names)) {
    stop(
      "'D' must have unique, non-empty column names (colnames(D)); ",
      "there is no 'groupef.names' argument to override this.",
      call. = FALSE
    )
  }

  if (!is.factor(group)) {
    stop(
      "'group' must be a factor (wrap with factor(group, levels = ...) ",
      "to control level order or supply a fixed superset of levels); ",
      "there is no 'group_levels' argument to override this.",
      call. = FALSE
    )
  }
  group_levels <- levels(group)
  if (length(group_levels) < 1L) {
    stop("'group' must have at least one level.", call. = FALSE)
  }

  if (is.null(group_name) || !nzchar(group_name)) {
    stop(
      "'group_name' could not be derived and was not supplied; ",
      "pass 'group_name' explicitly.",
      call. = FALSE
    )
  }

  if (!is.list(W) || is.data.frame(W)) {
    stop("'W' must be a list of design matrices.", call. = FALSE)
  }
  if (length(W) != length(re_names)) {
    stop("length(W) must equal ncol(D) = ", length(re_names), ".",
         call. = FALSE)
  }
  if (!setequal(names(W), re_names)) {
    stop(
      "names(W) must match colnames(D): ",
      paste(re_names, collapse = ", "), ".", call. = FALSE
    )
  }
  W <- W[re_names]

  list(
    n             = n,
    y             = y,
    D             = D,
    W             = W,
    tv_tol        = tv_tol,
    re_names      = re_names,
    group_levels  = group_levels,
    group_name    = group_name
  )
}

#' Validate the \code{sim_method} argument shared by LMM reg engines
#'
#' \code{"DEFAULT"} lets each route pick its own engine (exact iid sampling
#' for \code{rLMMNormal_reg_known_vcov}; two-block Gibbs everywhere else,
#' since it is the only engine available there). \code{"TWO_BLOCK_GIBBS"}
#' forces the two-block Gibbs engine even on routes where iid sampling would
#' otherwise be available.
#' @noRd
.rLMM_validate_sim_method <- function(sim_method, fn_name = "rLMMNormal_reg") {
  allowed <- c("DEFAULT", "TWO_BLOCK_GIBBS")
  if (!is.character(sim_method) || length(sim_method) != 1L ||
      is.na(sim_method) || !(sim_method %in% allowed)) {
    stop(
      fn_name, "(): 'sim_method' must be one of ",
      paste(paste0('"', allowed, '"'), collapse = ", "), ".",
      call. = FALSE
    )
  }
  sim_method
}

#' Validate a fixed-\eqn{\sigma^2} \code{prior_list} for the LMM Block~1 route
#'
#' Accepts either a single positive scalar (pooled \eqn{\sigma^2}, the
#' original contract) or, when \code{group_levels} is supplied, a numeric
#' vector of length \code{length(group_levels)} giving one fixed, known
#' dispersion per group. An unnamed length-\code{J} vector is accepted
#' positionally (already in \code{group_levels} order); a named vector must
#' match \code{group_levels} exactly (as a set) and is reordered accordingly.
#' @param group_levels Character vector of group levels, or \code{NULL} to
#'   only allow a pooled scalar (the original behavior).
#' @noRd
.rLMM_validate_fixed_dispersion_prior_list <- function(
    prior_list,
    fn_name = "rLMMNormal_reg",
    group_levels = NULL
) {
  if (!is.list(prior_list) || is.null(prior_list$dispersion)) {
    stop(
      fn_name, "(): 'dispprior_list' must contain 'dispersion' (fixed sigma^2).",
      call. = FALSE
    )
  }
  d <- prior_list$dispersion
  d <- .rLMM_validate_fixed_dispersion_vector(
    d, group_levels = group_levels, fn_name = fn_name,
    what = "prior_list$dispersion"
  )
  extra <- setdiff(
    names(prior_list),
    c("dispersion", "ddef")
  )
  if (length(extra)) {
    stop(
      fn_name, "(): 'dispprior_list' must contain fixed dispersion only; ",
      "unexpected fields: ", paste(extra, collapse = ", "), ".",
      call. = FALSE
    )
  }
  d
}

#' Shared positive-scalar-or-per-group-vector validator for fixed dispersion
#' @noRd
.rLMM_validate_fixed_dispersion_vector <- function(
    d,
    group_levels = NULL,
    fn_name = "rLMMNormal_reg",
    what = "dispersion"
) {
  J <- length(group_levels)
  if (!is.numeric(d) || length(d) < 1L || any(!is.finite(d)) || any(d <= 0)) {
    stop(
      fn_name, "(): '", what, "' must be a single positive number",
      if (J > 0L) paste0(" or a length-", J, " vector of positive per-group values"),
      ".",
      call. = FALSE
    )
  }
  if (length(d) == 1L) {
    return(as.numeric(d))
  }
  if (J < 1L || length(d) != J) {
    stop(
      fn_name, "(): '", what, "' must be a single positive number",
      if (J > 0L) paste0(" or a length-", J, " vector (one per group level: ",
                          paste(group_levels, collapse = ", "), ")"),
      ".",
      call. = FALSE
    )
  }
  nms <- names(d)
  if (!is.null(nms) && any(nzchar(nms))) {
    if (!setequal(nms, group_levels)) {
      stop(
        fn_name, "(): names('", what, "') must match the group levels (",
        paste(group_levels, collapse = ", "), ") exactly.",
        call. = FALSE
      )
    }
    d <- d[group_levels]
  }
  stats::setNames(as.numeric(d), group_levels)
}

#' @noRd
.rLMM_validate_dGamma_dispersion_prior_list <- function(
    prior_list,
    fn_name = "rLMMindepNormalGamma_reg"
) {
  if (!is.list(prior_list)) {
    stop(fn_name, "(): 'dispprior_list' must be a list.", call. = FALSE)
  }
  req <- c("shape", "rate", "beta", "Inv_Dispersion")
  miss <- req[!req %in% names(prior_list)]
  if (length(miss)) {
    stop(
      fn_name, "(): 'dispprior_list' must contain ",
      paste(req, collapse = ", "), " (from dGamma()).",
      call. = FALSE
    )
  }
  if (!isTRUE(prior_list$Inv_Dispersion)) {
    stop(
      fn_name, "(): dGamma() observation-dispersion prior requires ",
      "Inv_Dispersion = TRUE.",
      call. = FALSE
    )
  }
  shape <- prior_list$shape
  rate  <- prior_list$rate
  if (!is.numeric(shape) || length(shape) != 1L || !is.finite(shape) ||
      shape <= 0) {
    stop(fn_name, "(): 'dispprior_list$shape' must be a positive scalar.",
         call. = FALSE)
  }
  if (!is.numeric(rate) || length(rate) != 1L || !is.finite(rate) ||
      rate <= 0) {
    stop(fn_name, "(): 'dispprior_list$rate' must be a positive scalar.",
         call. = FALSE)
  }
  beta <- as.matrix(prior_list$beta)
  if (nrow(beta) != 1L || ncol(beta) != 1L) {
    stop(
      fn_name, "(): 'dispprior_list$beta' must be a 1 x 1 matrix for ",
      "observation-level dispersion.",
      call. = FALSE
    )
  }
  if (!is.null(prior_list$disp_lower) && !is.null(prior_list$disp_upper)) {
    if (prior_list$disp_upper <= prior_list$disp_lower) {
      stop(
        fn_name, "(): 'dispprior_list$disp_upper' must exceed 'disp_lower'.",
        call. = FALSE
      )
    }
  }
  prior_list
}

#' Plug-in observation \eqn{\sigma^2} from a dGamma / ING measurement \code{prior_list}
#' @noRd
.rLMM_dispersion_fix_from_prior_list <- function(
    prior_list,
    fn_name = "rLMM_reg"
) {
  if (!is.list(prior_list)) {
    stop(fn_name, "(): 'dispprior_list' must be a list.", call. = FALSE)
  }
  if (!is.null(prior_list$shape_group) || !is.null(prior_list$rate_group)) {
    shape_group <- as.numeric(prior_list$shape_group)
    rate_group  <- as.numeric(prior_list$rate_group)
    if (length(shape_group) < 1L || length(shape_group) != length(rate_group) ||
        any(!is.finite(shape_group)) || any(shape_group <= 0) ||
        any(!is.finite(rate_group)) || any(rate_group <= 0)) {
      stop(
        fn_name, "(): 'dispprior_list$shape_group' and 'dispprior_list$rate_group' ",
        "must be positive, finite, and of the same length.",
        call. = FALSE
      )
    }
    # Pooled reference (mean of per-group prior means); used only as an ICM
    # plug-in starting point, not as the sampler's actual per-group prior.
    return(mean(shape_group / rate_group))
  }
  if (is.null(prior_list$shape) || is.null(prior_list$rate)) {
    stop(
      fn_name, "(): 'dispprior_list' must contain 'shape' and 'rate' ",
      "(plug-in sigma^2 = shape / rate is derived internally).",
      call. = FALSE
    )
  }
  shape <- as.numeric(prior_list$shape[1L])
  rate  <- as.numeric(prior_list$rate[1L])
  if (!is.finite(shape) || shape <= 0 || !is.finite(rate) || rate <= 0) {
    stop(
      fn_name, "(): 'dispprior_list$shape' and 'dispprior_list$rate' must be positive scalars.",
      call. = FALSE
    )
  }
  as.numeric(shape / rate)
}

#' Observation-level linear predictor from group random effects
#' @noRd
.rLMM_observation_mu <- function(D, group, b_mat, group_levels) {
  g_chr <- as.character(group)
  g_idx <- match(g_chr, group_levels)
  if (anyNA(g_idx)) {
    stop("group levels not found in 'group_levels'.", call. = FALSE)
  }
  b_obs <- b_mat[g_idx, , drop = FALSE]
  rowSums(D * b_obs)
}

#' Group-level random-effect matrix from \code{coefficients} output
#' @noRd
.rLMM_b_matrix_from_coefficients <- function(
    coef_df,
    re_names,
    group_levels,
    group_name
) {
  b_mat <- matrix(
    NA_real_,
    nrow = length(group_levels),
    ncol = length(re_names),
    dimnames = list(group_levels, re_names)
  )
  for (lev in group_levels) {
    rows <- coef_df[[group_name]] == lev
    if (!any(rows)) {
      stop("coefficients missing group level: ", lev, call. = FALSE)
    }
    hit <- which(rows)[1L]
    b_mat[lev, ] <- as.numeric(coef_df[hit, re_names, drop = TRUE])
  }
  b_mat
}

#' ICM start for matrix-level LMM engines
#'
#' Iterated conditional modes for Block~2 hyperparameters at fixed Block~1
#' dispersion; used internally by \code{rLMMNormal_reg_*_vcov}.
#' @noRd
.rLMM_icm_at_start <- function(
    y,
    D,
    group,
    W,
    prior_list_block1,
    pfamily_list,
    re_names,
    group_levels,
    group_name,
    icm_tol,
    icm_maxit,
    verbose,
    engine_label
) {
  design_icm <- list(
    y             = y,
    D             = D,
    group         = factor(group, levels = group_levels),
    W             = W,
    groupef.names = re_names,
    group_name    = group_name
  )
  icm <- .two_block_icm_at_start(
    design       = design_icm,
    prior_list   = prior_list_block1,
    pfamily_list = pfamily_list,
    re_names     = re_names,
    family       = gaussian(),
    tol          = icm_tol,
    maxit        = icm_maxit
  )
  if (isTRUE(verbose)) {
    cat(sprintf(
      "  %s: ICM posterior mean (converged: %s, %d iter, delta = %.2e)\n\n",
      engine_label,
      icm$icm$converged,
      icm$icm$iterations,
      icm$icm$delta
    ))
  }
  icm
}

#' Labels and \code{convergence_info$method} for rate calibration plug-ins
#' @noRd
.rLMM_rate_calibration_meta <- function(
    any_non_normal,
    random_measurement = FALSE
) {
  if (isTRUE(random_measurement) && isTRUE(any_non_normal)) {
    list(
      method = "disp_upper_bound+disp_lower_bound",
      label  = paste0(
        "conservative: measurement disp_upper + ",
        "RE disp_lower plug-ins"
      )
    )
  } else if (isTRUE(random_measurement)) {
    list(
      method = "disp_upper_bound",
      label  = "conservative: measurement disp_upper plug-in"
    )
  } else if (isTRUE(any_non_normal)) {
    list(
      method = "disp_lower_bound",
      label  = "conservative: non-dNormal RE dispersion (disp_lower plug-in)"
    )
  } else {
    list(method = "exact", label = "exact")
  }
}

#' Calibrate inner Gibbs sweeps for matrix-level LMM engines
#'
#' Computes the two-block convergence rate at the chain start, then sets
#' \code{m_convergence} to at least \code{l_for_tv(tv_tol) + 1} (Theorem 3).
#' @noRd
.rLMM_calibrate_m_convergence <- function(
    D,
    group,
    W,
    prior_list_block1,
    pfamily_list,
    group_levels,
    tv_tol,
    any_non_normal,
    random_measurement = FALSE,
    engine_label,
    verbose
) {
  rate <- two_block_rate_from_pfamily_list(
    x                 = D,
    group             = group,
    x_hyper           = W,
    prior_list_block1 = prior_list_block1,
    pfamily_list      = pfamily_list,
    family            = gaussian(),
    group_levels      = group_levels
  )
  m_min <- .two_block_cap_inner_sweeps(
    two_block_l_for_tv(rate, tv_tol, method = "theorem3") + 1L
  )
  m_convergence <- m_min
  calib_meta <- .rLMM_rate_calibration_meta(
    any_non_normal     = any_non_normal,
    random_measurement = random_measurement
  )
  if (isTRUE(verbose)) {
    cat(sprintf(
      paste0(
        "--- %s: convergence calibration [%s]: lambda* = %.4f, ",
        "tv_tol = %g => m_min = %d, using m_convergence = %d ---\n\n"
      ),
      engine_label, calib_meta$label, rate$lambda_star, tv_tol, m_min,
      m_convergence
    ))
  }
  list(
    m_convergence     = m_convergence,
    convergence_info  = list(
      method        = calib_meta$method,
      tv_tol        = tv_tol,
      lambda_star   = rate$lambda_star,
      eigenvalues   = rate$eigenvalues,
      m_min         = m_min,
      m_convergence = m_convergence
    )
  )
}
#' Build a pooled ING measurement prior_list from a single dGamma() pfamily
#'
#' \code{mu}/\code{Sigma} are deliberately not part of the returned list --
#' see \code{.rLMM_validate_ing_measurement_prior_list()}'s header comment
#' for why the caller never supplies them.
#' @noRd
.rLMM_ing_measurement_prior_from_dGamma_pooled <- function(pf, fn_name) {
  if (!identical(pf$pfamily, "dGamma")) {
    stop(
      fn_name, "(): a 'dispprior_list' pfamily object must be dGamma() (pooled ",
      "measurement-dispersion prior); got \"", pf$pfamily, "\".",
      call. = FALSE
    )
  }
  pl <- pf$prior_list
  shape <- as.numeric(pl$shape[1L])
  rate  <- as.numeric(pl$rate[1L])
  if (!is.finite(shape) || shape <= 0) {
    stop(
      fn_name, "(): dGamma() 'dispprior_list$shape' must be a positive scalar.",
      call. = FALSE
    )
  }
  if (!is.finite(rate) || rate <= 0) {
    stop(
      fn_name, "(): dGamma() 'dispprior_list$rate' must be a positive scalar.",
      call. = FALSE
    )
  }
  disp_lower <- pl$disp_lower
  disp_upper <- pl$disp_upper
  if (!is.null(disp_lower) && !is.null(disp_upper) && disp_upper <= disp_lower) {
    stop(
      fn_name, "(): dGamma() 'dispprior_list$disp_upper' must exceed 'disp_lower'.",
      call. = FALSE
    )
  }
  if (is.null(disp_upper) || !is.numeric(disp_upper) || length(disp_upper) != 1L ||
      !is.finite(disp_upper) || disp_upper <= 0) {
    stop(
      fn_name, "(): dGamma() 'dispprior_list$disp_upper' is required ",
      "(conservative measurement-dispersion plug-in for lambda* calibration).",
      call. = FALSE
    )
  }
  out <- list(shape = shape, rate = rate, disp_upper = as.numeric(disp_upper))
  if (!is.null(disp_lower)) out$disp_lower <- as.numeric(disp_lower)
  out
}

#' Build a per-group ING measurement prior_list from a named list of
#' dGamma() pfamilies (one per group level, as returned by dGamma_list())
#' @noRd
.rLMM_ing_measurement_prior_from_dGamma_group <- function(
    prior_list,
    fn_name,
    group_levels
) {
  if (is.null(group_levels)) {
    stop(
      fn_name, "(): 'group_levels' is required to validate a per-group ",
      "measurement dispersion dispprior_list.",
      call. = FALSE
    )
  }
  nms <- names(prior_list)
  if (is.null(nms) || any(!nzchar(nms)) || !setequal(nms, group_levels)) {
    stop(
      fn_name, "(): names(dispprior_list) must match every group level (",
      paste(group_levels, collapse = ", "), ") exactly when 'dispprior_list' is ",
      "a named list of dGamma() pfamilies.",
      call. = FALSE
    )
  }
  shape_group      <- stats::setNames(numeric(length(group_levels)), group_levels)
  rate_group       <- stats::setNames(numeric(length(group_levels)), group_levels)
  disp_lower_group <- stats::setNames(numeric(length(group_levels)), group_levels)
  disp_upper_group <- stats::setNames(numeric(length(group_levels)), group_levels)
  for (lev in group_levels) {
    pf <- prior_list[[lev]]
    if (!identical(pf$pfamily, "dGamma")) {
      stop(
        fn_name, "(): prior_list[[\"", lev, "\"]] pfamily must be dGamma(); ",
        "got \"", pf$pfamily, "\".",
        call. = FALSE
      )
    }
    pl    <- pf$prior_list
    shape <- as.numeric(pl$shape[1L])
    rate  <- as.numeric(pl$rate[1L])
    lo    <- pl$disp_lower
    hi    <- pl$disp_upper
    if (!is.finite(shape) || shape <= 0 || !is.finite(rate) || rate <= 0) {
      stop(
        fn_name, "(): prior_list[[\"", lev, "\"]] must have positive, ",
        "finite 'shape' and 'rate'.",
        call. = FALSE
      )
    }
    if (is.null(lo) || is.null(hi) || !is.numeric(lo) || !is.numeric(hi) ||
        length(lo) != 1L || length(hi) != 1L ||
        !is.finite(lo) || !is.finite(hi) || lo <= 0 || hi <= lo) {
      stop(
        fn_name, "(): prior_list[[\"", lev, "\"]] must have a finite, ",
        "positive 'disp_lower' and 'disp_upper' with disp_upper > disp_lower.",
        call. = FALSE
      )
    }
    shape_group[[lev]]      <- shape
    rate_group[[lev]]       <- rate
    disp_lower_group[[lev]] <- as.numeric(lo)
    disp_upper_group[[lev]] <- as.numeric(hi)
  }
  list(
    shape_group      = shape_group,
    rate_group       = rate_group,
    disp_lower_group = disp_lower_group,
    disp_upper_group = disp_upper_group
  )
}

#' Validate a shared ING measurement prior for per-group Block~1 updates
#'
#' Accepts, in order of preference: (1) a single \code{dGamma()} pfamily
#' (pooled measurement dispersion), (2) a named list of \code{dGamma()}
#' pfamilies, one per \code{group_levels} entry, exactly as returned by
#' \code{\link{dGamma_list}()} (per-group measurement dispersion), or (3) the
#' legacy flat shape -- pooled \code{shape}/\code{rate}/\code{disp_lower}/
#' \code{disp_upper} (scalar), or per-group \code{shape_group}/
#' \code{rate_group}/\code{disp_lower_group}/\code{disp_upper_group} named
#' vectors (one entry per \code{group_levels}) -- kept only for backward
#' compatibility with \code{rlmerb()}'s own glue
#' (\code{.lmebayes_ing_measurement_prior_list()}/
#' \code{.lmebayes_ing_measurement_prior_list_group()}). Pooled vs. per-group
#' is mutually exclusive and determines the dispersion mode used by the
#' Block~1 sweep (see \code{.two_block_block1_ing_one_chain()}).
#'
#' \code{mu}/\code{Sigma} are never read here (and never required), even
#' from the legacy flat shape if present: \code{mu} is always \eqn{W_j\gamma}
#' (recomputed every sweep from \code{pfamily_list}'s own Block~2 hyper-prior,
#' never a fixed input), and \code{Sigma} -- the shared random-effects
#' covariance every group's \eqn{\beta_j} prior uses -- is filled in by the
#' caller (\code{rLMMindepNormalGamma_reg_known_vcov()}/
#' \code{_estimated_vcov()}) as \code{solve(P)} immediately after this
#' function returns, once \code{P} has been derived from \code{pfamily_list}
#' (see those functions' bodies); this guarantees \code{Sigma} can never
#' drift out of sync with \code{pfamily_list}'s own implied precision, which
#' requiring it here (and trusting whatever the caller supplied) could not.
#' @noRd
.rLMM_validate_ing_measurement_prior_list <- function(
    prior_list,
    p_re,
    fn_name = "rLMMindepNormalGamma_reg",
    group_levels = NULL
) {
  if (!is.list(prior_list)) {
    stop(fn_name, "(): 'dispprior_list' must be a list.", call. = FALSE)
  }

  if (inherits(prior_list, "pfamily")) {
    return(.rLMM_ing_measurement_prior_from_dGamma_pooled(prior_list, fn_name))
  }
  if (length(prior_list) > 0L &&
      all(vapply(prior_list, function(x) inherits(x, "pfamily"), logical(1L)))) {
    return(.rLMM_ing_measurement_prior_from_dGamma_group(
      prior_list, fn_name = fn_name, group_levels = group_levels
    ))
  }

  is_group <- !is.null(prior_list$shape_group) || !is.null(prior_list$rate_group)
  if (is_group) {
    if (is.null(group_levels)) {
      stop(
        fn_name, "(): 'group_levels' is required to validate a per-group ",
        "measurement dispersion dispprior_list.",
        call. = FALSE
      )
    }
    req_grp <- c("shape_group", "rate_group", "disp_lower_group", "disp_upper_group")
    miss <- req_grp[!req_grp %in% names(prior_list)]
    if (length(miss)) {
      stop(
        fn_name, "(): per-group 'dispprior_list' must contain ",
        paste(req_grp, collapse = ", "), ".",
        call. = FALSE
      )
    }
    for (nm in req_grp) {
      v <- prior_list[[nm]]
      if (!is.numeric(v) || is.null(names(v)) ||
          !setequal(names(v), group_levels)) {
        stop(
          fn_name, "(): 'dispprior_list$", nm, "' must be a numeric vector ",
          "named by every group level (", paste(group_levels, collapse = ", "),
          ").",
          call. = FALSE
        )
      }
      if (any(!is.finite(v)) || any(v <= 0)) {
        stop(
          fn_name, "(): 'dispprior_list$", nm, "' must be positive and finite ",
          "for every group.",
          call. = FALSE
        )
      }
    }
    if (any(prior_list$disp_upper_group[group_levels] <=
            prior_list$disp_lower_group[group_levels])) {
      stop(
        fn_name, "(): 'dispprior_list$disp_upper_group' must exceed ",
        "'disp_lower_group' for every group.",
        call. = FALSE
      )
    }
    return(prior_list)
  }

  req <- c("shape", "rate")
  miss <- req[!req %in% names(prior_list)]
  if (length(miss)) {
    stop(
      fn_name, "(): 'dispprior_list' must contain ",
      paste(req, collapse = ", "),
      " (from dIndependent_Normal_Gamma()).",
      call. = FALSE
    )
  }
  shape <- prior_list$shape
  rate  <- prior_list$rate
  if (!is.numeric(shape) || length(shape) != 1L || !is.finite(shape) ||
      shape <= 0) {
    stop(fn_name, "(): 'dispprior_list$shape' must be a positive scalar.",
         call. = FALSE)
  }
  if (!is.numeric(rate) || length(rate) != 1L || !is.finite(rate) ||
      rate <= 0) {
    stop(fn_name, "(): 'dispprior_list$rate' must be a positive scalar.",
         call. = FALSE)
  }
  if (!is.null(prior_list$disp_lower) && !is.null(prior_list$disp_upper)) {
    if (prior_list$disp_upper <= prior_list$disp_lower) {
      stop(
        fn_name, "(): 'dispprior_list$disp_upper' must exceed 'disp_lower'.",
        call. = FALSE
      )
    }
  }
  if (is.null(prior_list$disp_upper) || !is.numeric(prior_list$disp_upper) ||
      length(prior_list$disp_upper) != 1L || !is.finite(prior_list$disp_upper) ||
      prior_list$disp_upper <= 0) {
    stop(
      fn_name, "(): 'dispprior_list$disp_upper' is required (conservative ",
      "measurement-dispersion plug-in for lambda* calibration).",
      call. = FALSE
    )
  }
  prior_list
}

#' Block~1 Gaussian prior list (\code{P}, \code{dispersion})
#' @noRd
.rLMM_block1_prior_gaussian <- function(P, dispersion) {
  list(P = P, dispersion = dispersion, ddef = FALSE)
}

#' Conservative measurement \eqn{\sigma^2} plug-in for rate calibration
#' @noRd
.rLMM_measurement_disp_upper_for_rate <- function(ing_prior_list, fn_name) {
  d <- ing_prior_list$disp_upper
  if (is.null(d) || !is.numeric(d) || length(d) != 1L || !is.finite(d) ||
      d <= 0) {
    stop(
      fn_name, "(): 'dispprior_list$disp_upper' is required for lambda* ",
      "calibration when measurement dispersion is random.",
      call. = FALSE
    )
  }
  as.numeric(d)
}

#' Conservative measurement-dispersion inputs for lambda* calibration
#'
#' Scalar \code{shape}/\code{rate} priors plug in a single \code{disp_upper};
#' per-group priors instead build a per-observation precision \code{weights}
#' vector (\code{1 / disp_upper_group[group]}), so that groups with a wider
#' truncation window do not drag down the conservative rate bound for groups
#' with a tighter one. \code{dispersion_scalar} is a pooled fallback, used
#' only to satisfy generic \code{prior_list_block1} validation and as the
#' fixed \code{dispersion} plug-in when \code{weights} is not used downstream.
#' @noRd
.rLMM_measurement_rate_inputs <- function(
    ing_prior_list,
    group,
    group_levels,
    fn_name
) {
  if (!is.null(ing_prior_list$disp_upper_group)) {
    du <- ing_prior_list$disp_upper_group
    miss <- setdiff(group_levels, names(du))
    if (length(miss)) {
      stop(
        fn_name, "(): 'dispprior_list$disp_upper_group' is missing group ",
        "level(s): ", paste(miss, collapse = ", "), ".",
        call. = FALSE
      )
    }
    w <- 1 / as.numeric(du[as.character(group)])
    if (any(!is.finite(w))) {
      stop(
        fn_name, "(): non-finite per-observation weights derived from ",
        "'dispprior_list$disp_upper_group'.",
        call. = FALSE
      )
    }
    return(list(
      weights           = w,
      dispersion_scalar = mean(as.numeric(du))
    ))
  }
  list(
    weights           = NULL,
    dispersion_scalar = .rLMM_measurement_disp_upper_for_rate(ing_prior_list, fn_name)
  )
}

#' Truncated Gamma draw for ING measurement dispersion (prior-only path)
#' @noRd
.rLMM_ing_sample_sigma2 <- function(pl_j) {
  shape <- as.numeric(pl_j$shape[1L])
  rate  <- as.numeric(pl_j$rate[1L])
  lo    <- pl_j$disp_lower
  hi    <- pl_j$disp_upper
  if (is.null(lo)) {
    lo <- as.numeric(stats::qgamma(0.01, shape = shape, rate = rate))
  }
  if (is.null(hi)) {
    hi <- as.numeric(stats::qgamma(0.99, shape = shape, rate = rate))
  }
  F_lo <- stats::pgamma(lo, shape = shape, rate = rate)
  F_hi <- stats::pgamma(hi, shape = shape, rate = rate)
  stats::qgamma(
    F_lo + stats::runif(1L) * (F_hi - F_lo),
    shape = shape,
    rate  = rate
  )
}

#' Prior-only ING draw for one group when \code{Z_j} is rank-deficient
#' @noRd
.rLMM_ing_prior_draw_one_group <- function(pl_j, re_names) {
  sigma2 <- .rLMM_ing_sample_sigma2(pl_j)
  mu     <- as.numeric(pl_j$mu[, 1L])
  names(mu) <- re_names
  Sigma  <- as.matrix(pl_j$Sigma)
  L <- tryCatch(
    chol(Sigma),
    error = function(e) chol(Sigma + 1e-8 * diag(nrow(Sigma)))
  )
  b <- mu + sqrt(sigma2) * as.numeric(crossprod(L, stats::rnorm(length(mu))))
  names(b) <- re_names
  list(coefficients = b, dispersion = sigma2, iters = 1L)
}

#' @noRd
.rLMM_ing_prior_list_subset <- function(pl_j, keep, re_names) {
  keep <- as.integer(keep)
  re_keep <- re_names[keep]
  list(
    mu            = pl_j$mu[re_keep, , drop = FALSE],
    Sigma         = pl_j$Sigma[re_keep, re_keep, drop = FALSE],
    shape         = pl_j$shape,
    rate          = pl_j$rate,
    max_disp_perc = pl_j$max_disp_perc,
    disp_lower    = pl_j$disp_lower,
    disp_upper    = pl_j$disp_upper
  )
}

#' One school/group Block~1 ING draw (wrapper around \code{rindepNormalGamma_reg})
#'
#' Called once per factor level inside \code{.two_block_block1_ing_one_chain}.
#' Returns a **fixed contract** for the batch driver:
#' \code{list(coefficients = named b_j, dispersion = sigma2_j, iters = ...)}.
#'
#' Why not call \code{rindepNormalGamma_reg()} directly in the school loop?
#' \itemize{
#'   \item \code{rindepNormalGamma_reg} expects a full-rank design and returns an
#'     \code{rglmb}-style object (\code{coefficients} matrix, \code{Prior}, etc.).
#'     The sweep driver only needs one named \code{b_j} row plus scalar
#'     \code{sigma2_j} and envelope iteration count.
#'   \item Per-school \code{Z_j} can be rank-deficient (too few pupils, collinear
#'     RE columns). The ING envelope then fails or is undefined on the full
#'     \code{p_re} columns; this helper routes deficient schools to prior-only or
#'     QR-identifiable subsets before calling the sampler.
#'   \item Optional \code{full_rank} (from \code{design$groupef.rank}) skips repeated
#'     QR when rank was computed once at setup.
#' }
#'
#' Three paths (mutually exclusive after rank is resolved):
#' \describe{
#'   \item[Path A ? standard]{Full column rank and \code{n_j >= p_re}: one
#'     \code{rindepNormalGamma_reg(n = 1)} on all RE columns.}
#'   \item[Path B ? prior only]{\code{rank(Z_j) == 0}: no likelihood information;
#'     draw \code{(b_j, sigma2_j)} from the ING prior only
#'     (\code{.rLMM_ing_prior_draw_one_group}).}
#'   \item[Path C ? partial rank]{Some but not all columns identified: ING on the
#'     QR-identifiable columns, then impute the remaining coordinates from the
#'     prior given the drawn \code{sigma2_j}.}
#' }
#'
#' @noRd
.rLMM_ing_one_group_draw <- function(
    y_j,
    Z_j,
    pl_j,
    re_names,
    family,
    full_rank = NULL
) {
  p_re <- length(re_names)
  n_j  <- nrow(Z_j)

  .lmebayes_check_disp_bounds_or_stop(
    pl_j$disp_lower, pl_j$disp_upper, ".rLMM_ing_one_group_draw"
  )

  # Resolve effective rank of this school's design (cached or QR).
  rk   <- if (isTRUE(full_rank)) {
    p_re
  } else if (isFALSE(full_rank)) {
    0L
  } else {
    as.integer(Matrix::rankMatrix(Z_j, method = "qr")[1L])
  }

  # Path A: data identify all RE columns ? direct ING envelope on full Z_j.
  if (rk >= p_re && n_j >= p_re) {
    ing <- glmbayesCore::rindepNormalGamma_reg(
      n            = 1L,
      y            = y_j,
      x            = Z_j,
      prior_list   = pl_j,
      family       = family,
      progbar      = FALSE,
      verbose      = FALSE,
      use_parallel = FALSE
    )
    return(list(
      coefficients = ing$coefficients[1L, , drop = TRUE],
      dispersion   = as.numeric(ing$dispersion[1L]),
      iters        = as.numeric(ing$iters[1L])
    ))
  }

  # Path B: no identifiable column ? skip likelihood; prior draw only.
  if (rk < 1L) {
    return(.rLMM_ing_prior_draw_one_group(pl_j, re_names))
  }

  # Path C: rank 1 .. p_re-1 ? ING on identifiable columns, impute the rest.
  qr_j  <- qr(Z_j)
  keep  <- qr_j$pivot[seq_len(qr_j$rank)]
  pl_sub <- .rLMM_ing_prior_list_subset(pl_j, keep, re_names)
  .lmebayes_check_disp_bounds_or_stop(
    pl_sub$disp_lower, pl_sub$disp_upper, ".rLMM_ing_one_group_draw (Path C)"
  )
  ing <- glmbayesCore::rindepNormalGamma_reg(
    n            = 1L,
    y            = y_j,
    x            = Z_j[, keep, drop = FALSE],
    prior_list   = pl_sub,
    family       = family,
    progbar      = FALSE,
    verbose      = FALSE,
    use_parallel = FALSE
  )

  mu     <- as.numeric(pl_j$mu[, 1L])
  names(mu) <- re_names
  Sigma  <- as.matrix(pl_j$Sigma)
  sigma2 <- as.numeric(ing$dispersion[1L])
  b_full <- mu
  coef_sub <- ing$coefficients[1L, , drop = TRUE]
  b_full[keep] <- as.numeric(coef_sub)
  non <- setdiff(seq_len(p_re), keep)
  if (length(non)) {
    for (k in non) {
      b_full[k] <- mu[k] + sqrt(sigma2 * Sigma[k, k]) * stats::rnorm(1L)
    }
  }
  names(b_full) <- re_names

  list(
    coefficients = b_full,
    dispersion   = sigma2,
    iters        = as.numeric(ing$iters[1L])
  )
}

#' Build a group-level \code{prior_list} for \code{rindepNormalGamma_reg()}
#' @noRd
.rLMM_ing_prior_list_for_group <- function(
    ing_prior_list,
    mu_j,
    re_names
) {
  .lmebayes_check_disp_bounds_or_stop(
    ing_prior_list$disp_lower, ing_prior_list$disp_upper,
    ".rLMM_ing_prior_list_for_group (entry)"
  )
  mu_j <- as.numeric(mu_j)
  if (is.null(names(mu_j)) || !setequal(names(mu_j), re_names)) {
    names(mu_j) <- re_names
  }
  mu_mat <- matrix(mu_j[re_names], ncol = 1L,
                   dimnames = list(re_names, NULL))
  out <- list(
    mu            = mu_mat,
    Sigma         = ing_prior_list$Sigma,
    shape         = ing_prior_list$shape,
    rate          = ing_prior_list$rate,
    max_disp_perc = if (!is.null(ing_prior_list$max_disp_perc)) {
      ing_prior_list$max_disp_perc
    } else {
      0.99
    }
  )
  if (!is.null(ing_prior_list$disp_lower)) {
    out$disp_lower <- ing_prior_list$disp_lower
  }
  if (!is.null(ing_prior_list$disp_upper)) {
    out$disp_upper <- ing_prior_list$disp_upper
  }
  .lmebayes_check_disp_bounds_or_stop(
    out$disp_lower, out$disp_upper, ".rLMM_ing_prior_list_for_group (exit)"
  )
  out
}

#' Build Block~1 ING \code{prior_list} for one chain (mu_all from prep)
#' @noRd
.two_block_block1_ing_prior_list_one_chain <- function(mu_all, ing_prior_list) {
  prior_list <- list(
    mu            = mu_all,
    Sigma         = ing_prior_list$Sigma,
    shape         = ing_prior_list$shape,
    rate          = ing_prior_list$rate,
    max_disp_perc = if (!is.null(ing_prior_list$max_disp_perc)) {
      ing_prior_list$max_disp_perc
    } else {
      0.99
    }
  )
  if (!is.null(ing_prior_list$disp_lower)) {
    prior_list$disp_lower <- ing_prior_list$disp_lower
  }
  if (!is.null(ing_prior_list$disp_upper)) {
    prior_list$disp_upper <- ing_prior_list$disp_upper
  }
  prior_list
}

#' Block~1 one chain: block ING envelope update via \code{rIndepNormalGammaRegBlock}
#' (Centering ? Build ? DispersionBuild ? Sim in one C++ call).
#' @noRd
.two_block_block1_envelope_draw_one_chain <- function(
    y,
    Z,
    groups,
    prior_list,
    p_re,
    re_names,
    group_levels
) {
  nobs <- length(y)
  # Joint sampler (gs^k product-face table ? only feasible for very few groups;
  # for k groups with gs faces each the table has gs^k entries, exponential in k):
  # .rIndepNormalGammaRegBlock_cpp(
  #   n             = 1L,
  #   y             = y,
  #   x             = Z,
  #   block         = groups,
  #   prior_list    = prior_list,
  #   prior_lists   = NULL,
  #   offset        = rep(0, nobs),
  #   wt            = rep(1, nobs),
  #   p_re          = p_re,
  #   n_rss_iter    = 10L,
  #   Gridtype      = 3L,
  #   n_envopt      = -1L,
  #   RSS_ML        = NA_real_,
  #   use_parallel  = TRUE,
  #   use_opencl    = FALSE,
  #   progbar       = FALSE,
  #   verbose       = FALSE,
  #   group_levels  = group_levels,
  #   re_names      = re_names
  # )[c("b", "dispersion_ranef", "iters_mean")]
  .rIndepNormalGammaRegBlockInd_cpp(
    n             = 1L,
    y             = y,
    x             = Z,
    block         = groups,
    prior_list    = prior_list,
    prior_lists   = NULL,
    offset        = rep(0, nobs),
    wt            = rep(1, nobs),
    p_re          = p_re,
    n_rss_iter    = 10L,
    Gridtype      = 3L,
    n_envopt      = -1L,
    RSS_ML        = NA_real_,
    use_parallel  = TRUE,
    use_opencl    = FALSE,
    progbar       = FALSE,
    verbose       = FALSE,
    group_levels  = group_levels,
    re_names      = re_names
  )[c("b", "dispersion_ranef", "iters_mean")]
}

#' Block~1 one chain: independent per-group ING draws for a list of
#' per-group \code{dGamma()} measurement dispersion priors.
#'
#' Reuses \code{.rLMM_ing_one_group_draw()} (rank paths A/B/C) once per group,
#' each with its own \code{shape_group}/\code{rate_group}/\code{disp_lower_group}/
#' \code{disp_upper_group}. Unlike the shared-\eqn{\sigma^2} joint envelope
#' path, groups are conditionally independent given the current Block~2
#' fixef/tau2 (no joint acceptance across groups).
#' @noRd
.two_block_block1_ing_group_draw_one_chain <- function(
    y,
    Z,
    groups,
    ing_prior_list,
    mu_all,
    re_names,
    group_levels,
    family,
    full_rank = TRUE
) {
  p_re  <- length(re_names)
  g_chr <- as.character(groups)

  b_mat <- matrix(
    NA_real_,
    nrow = length(group_levels),
    ncol = p_re,
    dimnames = list(group_levels, re_names)
  )
  dispersion_ranef <- stats::setNames(numeric(length(group_levels)), group_levels)
  iters_j <- numeric(length(group_levels))

  for (j in seq_along(group_levels)) {
    lev  <- group_levels[j]
    rows <- which(g_chr == lev)

    pl_group_j <- list(
      Sigma         = ing_prior_list$Sigma,
      shape         = ing_prior_list$shape_group[[lev]],
      rate          = ing_prior_list$rate_group[[lev]],
      max_disp_perc = ing_prior_list$max_disp_perc,
      disp_lower    = ing_prior_list$disp_lower_group[[lev]],
      disp_upper    = ing_prior_list$disp_upper_group[[lev]]
    )
    .lmebayes_check_disp_bounds_or_stop(
      pl_group_j$disp_lower, pl_group_j$disp_upper,
      sprintf(".two_block_block1_ing_group_draw_one_chain (group '%s')", lev)
    )
    pl_j <- .rLMM_ing_prior_list_for_group(
      ing_prior_list = pl_group_j,
      mu_j           = mu_all[, lev],
      re_names       = re_names
    )
    draw_j <- .rLMM_ing_one_group_draw(
      y_j       = y[rows],
      Z_j       = Z[rows, , drop = FALSE],
      pl_j      = pl_j,
      re_names  = re_names,
      family    = family,
      full_rank = full_rank
    )

    b_mat[lev, ]            <- draw_j$coefficients[re_names]
    dispersion_ranef[[lev]] <- draw_j$dispersion
    iters_j[j]              <- draw_j$iters
  }

  list(
    b                = b_mat,
    dispersion_ranef = dispersion_ranef,
    iters_mean       = mean(iters_j),
    iters_group      = stats::setNames(iters_j, group_levels)
  )
}

#' Block~1 one chain: block ING envelope update for all groups
#'
#' Dispatches on the shape of \code{ing_prior_list}: a shared (pooled)
#' \code{shape}/\code{rate} routes through the joint C++ envelope
#' (\code{.two_block_block1_envelope_draw_one_chain()}, one shared
#' \eqn{\sigma^2} for all groups); a per-group \code{shape_group}/
#' \code{rate_group} (from a list of \code{dGamma()} pfamilies, one per
#' group) routes through \code{.two_block_block1_ing_group_draw_one_chain()}
#' instead (independent \eqn{\sigma^2_j} per group).
#' @noRd
.two_block_block1_ing_one_chain <- function(
    batch,
    i,
    design,
    block1_prior,
    ing_prior_list,
    family,
    ptypes
) {
  prep <- .two_block_block1_prep_one_chain(
    batch              = batch,
    i                  = i,
    design             = design,
    block1_prior       = block1_prior,
    ptypes             = ptypes,
    use_cpp_mu_all     = FALSE,
    use_cpp_prior_tau2 = FALSE
  )
  mu_all <- prep$mu_all
  p_re <- length(batch$re_names)
  re_names <- batch$re_names
  group_levels <- batch$group_levels

  y <- design$y
  Z <- as.matrix(design$D)

  if (!is.null(ing_prior_list$shape_group)) {
    out <- .two_block_block1_ing_group_draw_one_chain(
      y              = y,
      Z              = Z,
      groups         = design$group,
      ing_prior_list = ing_prior_list,
      mu_all         = mu_all,
      re_names       = re_names,
      group_levels   = group_levels,
      family         = family
    )
  } else {
    prior_list <- .two_block_block1_ing_prior_list_one_chain(
      mu_all, ing_prior_list
    )

    out <- .two_block_block1_envelope_draw_one_chain(
      y            = y,
      Z            = Z,
      groups       = design$group,
      prior_list   = prior_list,
      p_re         = p_re,
      re_names     = re_names,
      group_levels = group_levels
    )
  }

  if (isTRUE(getOption("glmbayesCore.debug_block1_ing_levels", FALSE))) {
    cat(sprintf(
      "  [Block1 ING debug] chain %d/%d: envelope sim dispersion=%s\n",
      i, batch$n, paste(signif(out$dispersion_ranef, 6), collapse = ", ")
    ))
    utils::flush.console()
  }

  out
}

#' Block~1 batch: ING per-group updates for all replicate chains
#' @noRd
.two_block_block1_ing_all_chains <- function(
    n,
    fixef,
    tau2,
    b,
    iters_ranef,
    re_names,
    group_levels,
    design,
    block1_prior,
    ing_prior_list,
    family,
    ptypes,
    iters_ranef_group = NULL,
    progbar = FALSE,
    progbar_prefix = "",
    progbar_finish_newline = TRUE
) {
  show_bar <- isTRUE(progbar) && n > 1L
  batch <- list(
    n            = n,
    fixef        = fixef,
    tau2         = tau2,
    re_names     = re_names,
    group_levels = group_levels,
    b            = .two_block_ensure_batch_b_dimnames(b, group_levels, re_names, n),
    iters_ranef  = iters_ranef + 0
  )
  group_mode <- !is.null(ing_prior_list$shape_group)
  dispersion_ranef <- if (group_mode) {
    matrix(
      NA_real_, nrow = n, ncol = length(group_levels),
      dimnames = list(NULL, group_levels)
    )
  } else {
    numeric(n)
  }
  if (group_mode) {
    if (is.null(iters_ranef_group)) {
      iters_ranef_group <- matrix(
        0,
        nrow = n,
        ncol = length(group_levels),
        dimnames = list(NULL, group_levels)
      )
    } else {
      iters_ranef_group <- iters_ranef_group + 0
    }
  }
  debug_b1 <- isTRUE(getOption("glmbayesCore.debug_block1_ing_levels", FALSE))

  for (i in seq_len(n)) {
    if (show_bar) .two_block_progress_bar(i, n, prefix = progbar_prefix)
    if (debug_b1 && nzchar(progbar_prefix)) {
      cat(progbar_prefix, "chain ", i, "/", n, " Block1 ING\n", sep = "")
      utils::flush.console()
    }
    out <- .two_block_block1_ing_one_chain(
      batch          = batch,
      i              = i,
      design         = design,
      block1_prior   = block1_prior,
      ing_prior_list = ing_prior_list,
      family         = family,
      ptypes         = ptypes
    )
    batch$b <- .two_block_batch_b_assign_slice(
      batch$b, i, out$b, use_cpp_b_slice = FALSE
    )
    batch$iters_ranef <- .two_block_batch_iters_ranef_add(
      batch$iters_ranef, i, out$iters_mean, use_cpp_iters_ranef_add = FALSE
    )
    if (group_mode) {
      dispersion_ranef[i, ] <- out$dispersion_ranef[group_levels]
      iters_ranef_group[i, ] <- iters_ranef_group[i, ] +
        out$iters_group[group_levels]
    } else {
      dispersion_ranef[i] <- out$dispersion_ranef
    }
  }
  if (show_bar) {
    .two_block_progress_bar_finish(newline = progbar_finish_newline)
  }

  list(
    b                = .two_block_ensure_batch_b_dimnames(
      batch$b, group_levels, re_names, n
    ),
    iters_ranef      = batch$iters_ranef,
    iters_ranef_group = if (group_mode) iters_ranef_group else NULL,
    dispersion_ranef = dispersion_ranef
  )
}

#' Pack ING sweep outputs (includes per-chain \code{dispersion_ranef} from envelope sim)
#' @noRd
.rGLMM_sweep_save_ing <- function(
    n,
    fixef,
    tau2,
    b,
    iters,
    iters_ranef,
    dispersion_ranef,
    re_names,
    group_levels,
    design,
    collect_block1 = TRUE,
    iters_ranef_group = NULL
) {
  out <- .rGLMM_sweep_save(
    n              = n,
    fixef          = fixef,
    tau2           = tau2,
    b              = b,
    iters          = iters,
    iters_ranef    = iters_ranef,
    re_names       = re_names,
    group_levels   = group_levels,
    design         = design,
    collect_block1 = collect_block1
  )
  out$dispersion_ranef <- dispersion_ranef
  if (!is.null(iters_ranef_group)) {
    out$iters_sigma2_draws <- iters_ranef_group
  }
  out
}

#' Two-block Gibbs sweep with per-group ING Block~1 measurement dispersion
#' @noRd
.rGLMM_sweep_ing_block1 <- function(
    n_chains,
    start_fixef,
    inner_sweeps,
    design,
    block1_prior,
    ing_prior_list,
    pfamily_list,
    family,
    re_names,
    group_levels,
    collect_block1 = TRUE,
    progbar        = FALSE,
    stage_label    = "",
    diag_sweeps    = FALSE,
    fixef_mode     = NULL,
    b_mode         = NULL,
    b_start        = NULL,
    ptypes         = NULL,
    tau2_start     = NULL,
    use_cpp_block2 = TRUE
) {
  if (is.null(ptypes)) {
    ptypes <- vapply(pfamily_list, function(pf) pf$pfamily, character(1))
    names(ptypes) <- re_names
  }

  if (is.null(tau2_start)) {
    tau2_start <- .two_block_tau2_start_from_pfamily(pfamily_list, re_names)
  } else {
    if (is.null(names(tau2_start)) || !setequal(names(tau2_start), re_names)) {
      stop("'tau2_start' must be a named vector with names(re_names).",
           call. = FALSE)
    }
    tau2_start <- as.numeric(tau2_start[re_names])
    names(tau2_start) <- re_names
  }
  if (is.null(b_start)) {
    if (is.null(b_mode)) {
      stop("'b_start' or 'b_mode' required for batch init.", call. = FALSE)
    }
    b_start <- b_mode
  }

  batch <- .rGLMM_sweep_initialize(
    n_chains     = n_chains,
    start_fixef  = start_fixef,
    b_start      = b_start,
    tau2_start   = tau2_start,
    re_names     = re_names,
    group_levels = group_levels
  )

  progbar_use <- isTRUE(progbar)
  sweep_stats <- vector("list", inner_sweeps)
  sweep_cov   <- vector("list", inner_sweeps)
  sweep_disp  <- vector("list", inner_sweeps)
  dispersion_ranef <- numeric(n_chains)
  group_mode <- !is.null(ing_prior_list$shape_group)
  iters_ranef_group <- if (group_mode) {
    matrix(
      0,
      nrow = n_chains,
      ncol = length(group_levels),
      dimnames = list(NULL, group_levels)
    )
  } else {
    NULL
  }

  for (m in seq_len(inner_sweeps)) {
    prefix_b1 <- if (progbar_use) {
      .two_block_progbar_prefix(stage_label, m, inner_sweeps, "Block1")
    } else {
      ""
    }
    prefix_b2 <- if (progbar_use) {
      .two_block_progbar_prefix(stage_label, m, inner_sweeps, "Block2")
    } else {
      ""
    }

    b1 <- .two_block_block1_ing_all_chains(
      n                      = batch$n,
      fixef                  = batch$fixef,
      tau2                   = batch$tau2,
      b                      = batch$b,
      iters_ranef            = batch$iters_ranef,
      re_names               = re_names,
      group_levels           = group_levels,
      design                 = design,
      block1_prior           = block1_prior,
      ing_prior_list         = ing_prior_list,
      family                 = family,
      ptypes                 = ptypes,
      iters_ranef_group      = iters_ranef_group,
      progbar                = progbar_use,
      progbar_prefix         = prefix_b1,
      progbar_finish_newline = FALSE
    )

    batch$b           <- b1$b
    batch$iters_ranef <- b1$iters_ranef
    dispersion_ranef  <- b1$dispersion_ranef
    if (group_mode) {
      iters_ranef_group <- b1$iters_ranef_group
    }

    b2 <- .two_block_block2_all_chains(
      n                      = batch$n,
      b                      = batch$b,
      fixef                  = batch$fixef,
      tau2                   = batch$tau2,
      iters                  = batch$iters,
      re_names               = re_names,
      group_levels           = group_levels,
      design                 = design,
      pfamily_list           = pfamily_list,
      ptypes                 = ptypes,
      use_cpp_block2         = use_cpp_block2,
      progbar                = progbar_use,
      progbar_prefix         = prefix_b2,
      progbar_finish_newline = (m == inner_sweeps)
    )
    batch$fixef <- b2$fixef
    batch$tau2  <- b2$tau2
    batch$iters <- b2$iters

    sweep_stats[[m]] <- .two_block_snapshot_fixef_stats(
      fixef    = batch$fixef,
      re_names = re_names
    )
    sweep_cov[[m]] <- .two_block_snapshot_fixef_cov(
      fixef    = batch$fixef,
      re_names = re_names
    )
    sweep_disp[[m]] <- .two_block_snapshot_disp_stats(
      tau2     = batch$tau2,
      ptypes   = ptypes,
      re_names = re_names
    )
    if (progbar_use && n_chains <= 1L) {
      prefix_sweep <- if (nzchar(stage_label)) {
        sprintf("[%s] sweep %d/%d: ", stage_label, m, inner_sweeps)
      } else {
        sprintf("sweep %d/%d: ", m, inner_sweeps)
      }
      .two_block_progress_bar(m, inner_sweeps, prefix = prefix_sweep)
      .two_block_progress_bar_finish(newline = (m == inner_sweeps))
    }
  }

  out <- .rGLMM_sweep_save_ing(
    n              = batch$n,
    fixef          = batch$fixef,
    tau2           = batch$tau2,
    b              = batch$b,
    iters          = batch$iters,
    iters_ranef    = batch$iters_ranef,
    iters_ranef_group = iters_ranef_group,
    dispersion_ranef = dispersion_ranef,
    re_names       = re_names,
    group_levels   = group_levels,
    design         = design,
    collect_block1 = collect_block1
  )
  out$sweep_history <- .two_block_build_sweep_history(
    stage_label = stage_label,
    sweep_stats = sweep_stats,
    fixef_mode  = fixef_mode,
    re_names    = re_names,
    sweep_cov   = sweep_cov,
    disp_stats  = sweep_disp
  )
  if (isTRUE(diag_sweeps)) {
    print(out$sweep_history)
  }
  invisible(out)
}

#' Format ING sweep output for staged \code{fixef.*} naming
#' @noRd
.rLMM_format_ing_sweep_out <- function(
    sweep_out,
    n,
    re_names,
    group_levels,
    fixef_mode,
    fixef_init
) {
  staged <- .rLMM_format_v2_out(
    v2_out       = sweep_out,
    n            = n,
    re_names     = re_names,
    group_levels = group_levels,
    fixef_mode   = fixef_mode,
    fixef_init   = fixef_init
  )
  staged[["group.dispersion"]] <- sweep_out$dispersion_ranef
  staged[["group.dispersion.mean"]] <- if (is.matrix(sweep_out$dispersion_ranef)) {
    colMeans(sweep_out$dispersion_ranef)
  } else {
    mean(sweep_out$dispersion_ranef)
  }
  if (!is.null(sweep_out$iters_sigma2_draws)) {
    staged[["group.dispersion.iters"]] <- sweep_out$iters_sigma2_draws
  }
  staged$sweep_history         <- sweep_out$sweep_history
  staged
}

#' ING measurement LMM: pilot then main via ING sweep-outer engine
#' @noRd
.rLMMIngNormal_reg_run_with_pilot <- function(
    inp,
    group,
    P,
    ing_prior_list,
    pfamily_list,
    pf_summary,
    icm_tol,
    icm_maxit,
    progbar,
    verbose,
    stage_verbose = FALSE,
    gap_tol       = 0.0196,
    mode_gap_max  = 1.0,
    diag_sweeps   = FALSE,
    any_non_normal = TRUE,
    random_measurement = TRUE,
    engine_label  = "rLMMindepNormalGamma_reg_estimated_vcov",
    result_class  = "rLMMindepNormalGamma_reg_estimated_vcov",
    cl
) {
  re_names         <- inp$re_names
  group_levels     <- inp$group_levels
  group_name       <- inp$group_name
  tv_tol           <- inp$tv_tol
  n                <- inp$n
  n_pilot_arg      <- NULL
  m_convergence_pilot <- NULL
  rate_calibration <- NULL
  collect_block1   <- TRUE
  family           <- gaussian()
  is_gaussian      <- TRUE
  ptypes           <- pf_summary$ptypes

  dispersion_fix <- .rLMM_dispersion_fix_from_prior_list(
    ing_prior_list, fn_name = engine_label
  )

  rate_inputs <- .rLMM_measurement_rate_inputs(
    ing_prior_list, group, group_levels, engine_label
  )
  prior_list_block1_icm <- .rLMM_block1_prior_gaussian(P, dispersion_fix)
  prior_list_block1_rate <- .rLMM_block1_prior_gaussian(P, rate_inputs$dispersion_scalar)
  .two_block_validate_block1_prior(prior_list_block1_icm, family = family)
  .two_block_validate_block1_prior(prior_list_block1_rate, family = family)

  calib_meta <- .rLMM_rate_calibration_meta(
    any_non_normal     = any_non_normal,
    random_measurement = random_measurement
  )

  gap_tol <- .two_block_validate_gap_tol(gap_tol)

  will_pilot <- .two_block_pilot_will_run(
    is_gaussian,
    n_pilot_arg,
    gap_tol,
    tv_tol,
    any_non_normal     = any_non_normal,
    random_measurement = random_measurement
  )
  run_pilot <- will_pilot
  run_ub    <- will_pilot && !is.null(tv_tol)

  if (run_pilot && is.null(m_convergence_pilot)) {
    m_convergence_pilot <- if (!is.null(tv_tol)) {
      NULL
    } else {
      10L
    }
  }

  icm <- .rLMM_icm_at_start(
    y                     = inp$y,
    D                     = inp$D,
    group                 = group,
    W                     = inp$W,
    prior_list_block1     = prior_list_block1_icm,
    pfamily_list          = pfamily_list,
    re_names              = re_names,
    group_levels          = group_levels,
    group_name            = group_name,
    icm_tol               = icm_tol,
    icm_maxit             = icm_maxit,
    verbose               = verbose,
    engine_label          = engine_label
  )
  fixef_mode <- icm$start
  ranef_mode <- icm$b_start
  icm_info   <- icm$icm
  if (isTRUE(verbose)) {
    cat(sprintf(
      "  %s: Block 2 start at lmer tau^2 plug-in (converged: %s, %d iter, delta = %.2e)\n\n",
      engine_label,
      icm_info$converged,
      icm_info$iterations,
      icm_info$delta
    ))
  }
  b_start <- ranef_mode

  design <- list(
    y             = inp$y,
    D             = inp$D,
    group         = factor(group, levels = group_levels),
    W             = inp$W,
    groupef.names = re_names,
    group_name    = group_name,
    groupef.rank  = .lmebayes_groupef_rank_from_Z(
      inp$D, group, group_levels = group_levels
    )
  )

  if (is.null(b_start)) {
    b_start <- matrix(
      0,
      nrow = length(group_levels),
      ncol = length(re_names),
      dimnames = list(group_levels, re_names)
    )
  }

  fixef_mode_ref <- fixef_mode
  b_mode_ref     <- b_start
  progbar_use    <- isTRUE(progbar) || isTRUE(verbose) || isTRUE(stage_verbose)

  rate <- two_block_rate_from_pfamily_list(
    x                 = inp$D,
    group             = group,
    x_hyper           = inp$W,
    prior_list_block1 = prior_list_block1_rate,
    pfamily_list      = pfamily_list,
    weights           = rate_inputs$weights,
    family            = family,
    group_levels      = group_levels
  )

  m_min <- NULL
  if (!is.null(tv_tol)) {
    m_min <- .two_block_cap_inner_sweeps(
    two_block_l_for_tv(rate, tv_tol, method = "theorem3") + 1L
  )
  }

  p_dim            <- sum(vapply(fixef_mode, length, integer(1L)))
  D_max            <- if (!is.null(mode_gap_max)) sqrt(p_dim) * mode_gap_max else 0
  m_pilot_from_gap <- NULL

  if (run_pilot && is.null(m_convergence_pilot) && !is.null(tv_tol)) {
    erf1_inv_tv <- stats::qnorm((tv_tol + 1) / 2) / sqrt(2)
    c_tol       <- erf1_inv_tv * 2 * sqrt(2)
    m_pilot_from_gap <- .two_block_m_pilot_from_gap(rate, D_max, c_tol, m_min)
    m_convergence_pilot <- m_pilot_from_gap
  }

  pilot_plan <- .two_block_resolve_pilot_plan(
    is_gaussian         = is_gaussian,
    n                   = n,
    n_pilot_arg         = n_pilot_arg,
    gap_tol             = gap_tol,
    tv_tol              = tv_tol,
    m_convergence_user  = NULL,
    m_convergence_pilot = m_convergence_pilot,
    rate                = rate,
    p_dim               = p_dim,
    m_min               = m_min,
    any_non_normal      = any_non_normal,
    random_measurement  = random_measurement
  )
  n_pilot        <- pilot_plan$n_pilot
  m_convergence  <- pilot_plan$m_convergence
  pilot_cost_opt <- pilot_plan$pilot_cost_opt
  run_pilot      <- n_pilot > 0L
  run_ub         <- run_pilot && !is.null(tv_tol)

  if (is.null(m_min) && !run_pilot) {
    m_convergence <- 10L
  }

  if (is.null(rate_calibration) && !is.null(tv_tol)) {
    rate_calibration <- list(
      lambda_star = rate$lambda_star,
      eigenvalues = rate$eigenvalues,
      m_min       = m_min
    )
  }

  calib_label <- calib_meta$label

  if (isTRUE(verbose) && !is.null(tv_tol)) {
    cat(sprintf(
      paste0(
        "--- %s: convergence calibration [%s]:\n",
        "    lambda* = %.4f, tv_tol = %g => m_min = %d (mode start), ",
        "main m_convergence = %d ---\n\n"
      ),
      engine_label, calib_label, rate$lambda_star, tv_tol, m_min, m_convergence
    ))
    if (run_pilot && !is.null(mode_gap_max) && !is.null(m_pilot_from_gap)) {
      cat(sprintf(
        paste0(
          "--- %s: pilot sweep calibration [mode_gap_max = %g SD/dim, p = %d, ",
          "D_max = %.4f]:\n    m_min = %d, lambda* = %.4f => ",
          "m_convergence_pilot = %d ---\n\n"
        ),
        engine_label, mode_gap_max, p_dim, D_max, m_min,
        rate$lambda_star, m_convergence_pilot
      ))
    }
    if (run_pilot) {
      .two_block_print_pilot_plan(
        pilot_plan          = pilot_plan,
        n                   = n,
        m_convergence_pilot = m_convergence_pilot,
        rate                = rate,
        tv_tol              = tv_tol,
        p                   = p_dim,
        verbose             = verbose
      )
    }
  }

  convergence_info <- list(
    method              = calib_meta$method,
    tv_tol              = tv_tol,
    gap_tol             = gap_tol,
    n_pilot             = n_pilot,
    n_pilot_source      = pilot_plan$n_pilot_source,
    n_pilot_gap_tol     = pilot_plan$n_pilot_gap_tol,
    lambda_star         = rate$lambda_star,
    eigenvalues         = rate$eigenvalues,
    m_min               = m_min,
    m_certificate       = pilot_plan$m_certificate,
    m_convergence       = m_convergence,
    m_convergence_pilot = if (run_pilot) m_convergence_pilot else NULL,
    mode_gap_max        = if (run_pilot) mode_gap_max else NULL,
    m_pilot_from_gap    = if (run_pilot) m_pilot_from_gap else NULL,
    pilot_cost_opt      = pilot_cost_opt,
    draw_engine         = "rGLMM_sweep_ing_block1_ind"
  )

  m_convergence_used <- m_convergence
  fixef_init         <- fixef_mode
  pilot_res          <- NULL
  pilot_chisq        <- NULL
  pilot_ub           <- NULL
  tau2_start_main    <- if (!is.null(icm) && !is.null(icm$tau2_start)) {
    icm$tau2_start
  } else {
    .two_block_tau2_start_from_pfamily(pfamily_list, re_names)
  }

  sweep_common <- list(
    design         = design,
    block1_prior   = prior_list_block1_icm,
    ing_prior_list = ing_prior_list,
    pfamily_list   = pfamily_list,
    family         = family,
    re_names       = re_names,
    group_levels   = group_levels,
    collect_block1 = collect_block1,
    progbar        = progbar_use,
    fixef_mode     = fixef_mode_ref,
    b_mode         = b_mode_ref,
    b_start        = b_mode_ref,
    ptypes         = ptypes,
    use_cpp_block2 = TRUE
  )

  if (run_pilot) {
    if (isTRUE(verbose)) {
      cat(sprintf(
        "--- %s [sweep-outer]: pilot stage (%d chains; m_convergence_pilot = %d) ---\n\n",
        engine_label, n_pilot, m_convergence_pilot
      ))
    }

    pilot_raw <- do.call(
      .rGLMM_sweep_ing_block1,
      c(
        list(
          n_chains     = n_pilot,
          start_fixef  = fixef_mode,
          inner_sweeps = m_convergence_pilot,
          stage_label  = "pilot",
          diag_sweeps  = isTRUE(diag_sweeps),
          tau2_start   = tau2_start_main
        ),
        sweep_common
      )
    )

    pilot_chisq <- .two_block_pilot_chisq_test(
      fixef_draws = pilot_raw$fixef_draws,
      re_names    = re_names,
      fixef_mode  = fixef_mode,
      n_pilot     = n_pilot
    )

    if (isTRUE(stage_verbose) || isTRUE(verbose)) {
      cat(sprintf(
        "--- %s: pilot vs mode chi-squared test: p = %.4g (df = %d, n_pilot = %d) ---\n\n",
        engine_label,
        pilot_chisq$p_value, pilot_chisq$df, pilot_chisq$n_pilot
      ))
    }

    fixef_init <- .two_block_fixef_colmeans(
      pilot_raw$fixef_draws, re_names, fixef_mode
    )

    tau2_start_main <- .two_block_tau2_start_from_dispersion_draws(
      pilot_raw$dispersion_fixef_draws, re_names
    )

    if (run_ub) {
      pilot_ub <- .two_block_pilot_ub_from_coefficients(
        pilot_coefficients = pilot_raw$coefficients,
        n_pilot            = n_pilot,
        re_names           = re_names,
        group_levels       = group_levels,
        group_name         = group_name,
        x                  = inp$D,
        group              = group,
        x_hyper            = inp$W,
        prior_list         = prior_list_block1_rate,
        pfamily_list       = pfamily_list,
        family             = family,
        tv_tol             = tv_tol,
        dispersion         = rate_inputs$dispersion_scalar,
        weights            = rate_inputs$weights
      )
      if (pilot_ub$m_min_upper > m_convergence_used) {
        m_convergence_used <- pilot_ub$m_min_upper
      }
      convergence_info$lambda_star_upper <- pilot_ub$rate_upper$lambda_star
      convergence_info$eigenvalues_upper <- pilot_ub$max_eigenvalues
      convergence_info$m_min_upper       <- pilot_ub$m_min_upper
      convergence_info$i_max_rate        <- pilot_ub$i_max_rate
      convergence_info$lambda_star_vec   <- pilot_ub$lambda_star_vec
      convergence_info$m_convergence     <- m_convergence_used
    }

    ## Omega-marginalized (Section 16) safeguard: known_vcov engine only
    ## (any_non_normal = FALSE => Lambda fixed/known, the assumption the
    ## Section 16 derivation relies on). When the marginalized lambda_star
    ## is computable on every pilot draw (no Lambda + H_j PD failures) and
    ## stays below the safeguard cutoff, it is COMBINED with the plain
    ## disp_upper plug-in envelope (.two_block_combine_rate_envelopes():
    ## rank-matched pointwise max of the two eigenvalue spectra, then a
    ## single TV certificate from the combined spectrum) rather than
    ## replacing it outright -- neither safeguard's own blind spot (a
    ## fixed, possibly-loose disp_upper plug-in on one side; a finite
    ## pilot sample that may miss a worse unobserved state on the other)
    ## can then make m_convergence_used smaller than either alone would
    ## have required. Otherwise the plain-rate m_convergence_used computed
    ## above is kept and a warning is issued.
    marginal_ub <- NULL
    if (run_ub && !any_non_normal) {
      ## The main stage starts from the pilot MEAN (fixef_init), not the
      ## exact mode, so the override must certify the SAME pilot-start D0
      ## the plain-rate m_convergence above used (two_block_
      ## m_convergence_for_pilot_start()'s D0 = qchisq(pilot_start_tol, p) /
      ## n_pilot) -- D0 = 0 would only certify a mode start and silently
      ## understate the sweeps needed.
      D0_main <- two_block_d0_pilot_start(
        n_pilot = n_pilot, p = p_dim, pilot_start_tol = 0.95
      )
      marginal_ub <- .two_block_pilot_marginal_ub_from_coefficients(
        pilot_coefficients = pilot_raw$coefficients,
        n_pilot            = n_pilot,
        y                  = inp$y,
        D                  = inp$D,
        group              = group,
        x_hyper            = inp$W,
        prior_list_block1  = prior_list_block1_rate,
        pfamily_list       = pfamily_list,
        tv_tol             = tv_tol,
        group_name         = group_name,
        group_levels       = group_levels,
        re_names           = re_names,
        ing_prior_list     = ing_prior_list,
        D0                 = D0_main
      )
      convergence_info$marginal_D0             <- D0_main
      convergence_info$lambda_star_marginal    <- marginal_ub$lambda_star_marginal
      convergence_info$m_min_marginal          <- marginal_ub$m_min_marginal
      convergence_info$marginal_rate_valid     <- marginal_ub$valid
      convergence_info$marginal_pd_fail_groups <- marginal_ub$failing_groups
      convergence_info$marginal_n_skipped      <- marginal_ub$n_skipped
      convergence_info$marginal_cutoff         <- marginal_ub$cutoff
      convergence_info$marginal_fallback_message <- NULL

      if (isTRUE(marginal_ub$valid)) {
        combined_ub <- .two_block_combine_rate_envelopes(
          rate_a = pilot_ub$rate_upper,
          rate_b = marginal_ub$rate_marginal,
          tv_tol = tv_tol,
          D0     = D0_main
        )
        if (combined_ub$m_min_combined > m_convergence_used) {
          m_convergence_used <- combined_ub$m_min_combined
        }
        convergence_info$m_convergence        <- m_convergence_used
        convergence_info$lambda_star_combined <- combined_ub$rate_combined$lambda_star
        convergence_info$eigenvalues_combined <- combined_ub$rate_combined$eigenvalues
        convergence_info$m_min_combined       <- combined_ub$m_min_combined
        if (isTRUE(verbose) || isTRUE(stage_verbose)) {
          cat(sprintf(
            paste0(
              "--- %s: Omega-marginalized safeguard [%d pilot draws, all groups' Lambda + H_j PD]:\n",
              "    lambda_star_marginal = %.4f (< %.2f); combined (rank-matched max) with the\n",
              "    disp_upper plug-in envelope => lambda_star_combined = %.4f => main\n",
              "    m_convergence = %d ---\n\n"
            ),
            engine_label, n_pilot, marginal_ub$lambda_star_marginal,
            marginal_ub$cutoff, combined_ub$rate_combined$lambda_star, m_convergence_used
          ))
        }
      } else {
        reason <- if (length(marginal_ub$failing_groups)) {
          sprintf(
            "outlier groups detected (%s; %d / %d pilot draws skipped)",
            paste(marginal_ub$failing_groups, collapse = ", "),
            marginal_ub$n_skipped, n_pilot
          )
        } else {
          sprintf(
            "lambda_star_marginal = %.4f is at/above the safeguard cutoff (%.2f)",
            marginal_ub$lambda_star_marginal, marginal_ub$cutoff
          )
        }
        fallback_msg <- sprintf(
          paste0(
            "%s - Omega-marginalized convergence bound may not be valid; ",
            "falling back to lambda_star = %.4f (m_convergence = %d)."
          ),
          reason, convergence_info$lambda_star_upper, m_convergence_used
        )
        convergence_info$marginal_fallback_message <- fallback_msg
        warning(engine_label, "(): ", fallback_msg, call. = FALSE)
        if (isTRUE(verbose) || isTRUE(stage_verbose)) {
          cat(sprintf("--- %s: %s ---\n\n", engine_label, fallback_msg))
        }
      }
    }

    if (isTRUE(stage_verbose) && run_ub) {
      .two_block_print_pilot_stage_diagnostics(
        n_pilot            = n_pilot,
        n_main             = n,
        pilot_ub           = pilot_ub,
        rate_calibration   = rate_calibration,
        m_convergence_used = m_convergence_used
      )
    } else if (isTRUE(verbose)) {
      cat(sprintf(
        "--- %s [sweep-outer]: pilot complete; main stage (%d chains; m_convergence = %d) ---\n\n",
        engine_label, n, m_convergence_used
      ))
    }

    pilot_res <- .rLMM_format_ing_sweep_out(
      sweep_out    = pilot_raw,
      n            = n_pilot,
      re_names     = re_names,
      group_levels = group_levels,
      fixef_mode   = fixef_mode,
      fixef_init   = fixef_mode
    )
  } else if (isTRUE(verbose)) {
    cat(sprintf(
      "--- %s [sweep-outer]: main stage (%d chains; m_convergence = %d) ---\n\n",
      engine_label, n, m_convergence_used
    ))
  }

  main_raw <- do.call(
    .rGLMM_sweep_ing_block1,
    c(
      list(
        n_chains     = n,
        start_fixef  = fixef_init,
        inner_sweeps = m_convergence_used,
        stage_label  = "main",
        diag_sweeps  = isTRUE(diag_sweeps),
        tau2_start   = tau2_start_main
      ),
      sweep_common
    )
  )

  draw_engine_args <- c(
    list(
      n_chains     = n,
      start_fixef  = fixef_init,
      inner_sweeps = m_convergence_used,
      stage_label  = "main",
      diag_sweeps  = isTRUE(diag_sweeps),
      tau2_start   = tau2_start_main
    ),
    sweep_common
  )

  main_res <- .rLMM_format_ing_sweep_out(
    sweep_out    = main_raw,
    n            = n,
    re_names     = re_names,
    group_levels = group_levels,
    fixef_mode   = fixef_mode,
    fixef_init   = fixef_init
  )

  .lmebayes_assemble_reg_result(
    staged              = main_res,
    call                = cl,
    m_convergence       = m_convergence_used,
    convergence_info    = convergence_info,
    pfamily_list        = pfamily_list,
    dispprior_list      = ing_prior_list,
    family              = family,
    groupef.mode        = ranef_mode,
    any_non_normal      = any_non_normal,
    design              = design,
    result_class        = result_class,
    parent_class        = "rLMMindepNormalGamma_reg",
    draw_engine         = "rGLMM_sweep_ing_block1_ind",
    sim_method_used     = "TWO_BLOCK_GIBBS",
    icm_info            = icm_info,
    pilot_draws         = if (run_pilot) pilot_res else NULL,
    n_pilot             = if (run_pilot) n_pilot else NULL,
    m_convergence_pilot = if (run_pilot) m_convergence_pilot else NULL,
    pilot_chisq         = if (run_pilot) pilot_chisq else NULL,
    pilot_ub            = if (run_ub) pilot_ub else NULL,
    tv_tol              = if (run_ub) tv_tol else NULL,
    offset              = inp$offset,
    weights             = if (is.null(inp$weights)) 1 else inp$weights
  )
}

#' Shared sampling pipeline for \code{rLMMNormal_reg_*_vcov} routes
#' @noRd
.rLMMNormal_reg_run <- function(
    inp,
    group,
    P,
    dispersion,
    pfamily_list,
    pf_summary,
    icm_tol,
    icm_maxit,
    progbar,
    verbose,
    engine_label,
    any_non_normal,
    draw_engine,
    result_class,
    cl,
    fixef_start = NULL
) {
  prior_list_block1 <- list(
    P          = P,
    dispersion = dispersion,
    ddef       = FALSE
  )
  .two_block_validate_block1_prior(prior_list_block1, family = gaussian())

  icm_info   <- NULL
  ranef_mode <- NULL

  if (is.null(fixef_start)) {
    icm <- .rLMM_icm_at_start(
      y                 = inp$y,
      D                 = inp$D,
      group             = group,
      W                 = inp$W,
      prior_list_block1 = prior_list_block1,
      pfamily_list      = pfamily_list,
      re_names          = inp$re_names,
      group_levels      = inp$group_levels,
      group_name        = inp$group_name,
      icm_tol           = icm_tol,
      icm_maxit         = icm_maxit,
      verbose           = verbose,
      engine_label      = engine_label
    )
    fixef_mode <- icm$start
    ranef_mode <- icm$b_start
    icm_info   <- icm$icm
  } else {
    if (!is.list(fixef_start) || is.null(names(fixef_start))) {
      stop("'fixef_start' must be a named list.", call. = FALSE)
    }
    if (!setequal(names(fixef_start), inp$re_names)) {
      stop("names(fixef_start) must match groupef.names.", call. = FALSE)
    }
    fixef_mode <- fixef_start[inp$re_names]
  }

  calib <- .rLMM_calibrate_m_convergence(
    D                 = inp$D,
    group             = group,
    W                 = inp$W,
    prior_list_block1 = prior_list_block1,
    pfamily_list      = pfamily_list,
    group_levels      = inp$group_levels,
    tv_tol            = inp$tv_tol,
    any_non_normal    = any_non_normal,
    engine_label      = engine_label,
    verbose           = verbose
  )
  m_convergence    <- calib$m_convergence
  convergence_info <- calib$convergence_info
  convergence_info$draw_engine <- draw_engine

  if (is.null(ranef_mode)) {
    ## fixef_start supplied directly (e.g. legacy rLMMindepNormalGamma_reg()
    ## outer loop): no ICM b-mode available. Block~1 draws b from its exact
    ## full conditional given fixef/tau2 each sweep, so any placeholder here
    ## only matters for array shape/dimnames, not for the sampled values.
    ranef_mode <- matrix(
      0,
      nrow = length(inp$group_levels),
      ncol = length(inp$re_names),
      dimnames = list(inp$group_levels, inp$re_names)
    )
  }

  design <- list(
    y             = inp$y,
    D             = inp$D,
    group         = factor(group, levels = inp$group_levels),
    W             = inp$W,
    groupef.names = inp$re_names,
    group_name    = inp$group_name
  )

  ## Consistent with the other sweep-outer engines (.rLMMNormal_reg_run_with_pilot(),
  ## .rGLMM_sweep_ing_block1()): 'verbose' alone also enables the progress bar.
  progbar_use <- isTRUE(progbar) || isTRUE(verbose)

  ## Known vcov: tau2 is fixed at the pfamily plug-in values throughout (all
  ## Block~2 components are dNormal(), so .two_block_block2_all_chains() never
  ## resamples tau2).
  tau2_start <- .two_block_tau2_start_from_pfamily(pfamily_list, inp$re_names)

  raw <- rGLMM_sweep(
    n_chains       = inp$n,
    start_fixef    = fixef_mode,
    inner_sweeps   = m_convergence,
    design         = design,
    block1_prior   = prior_list_block1,
    pfamily_list   = pfamily_list,
    family         = gaussian(),
    re_names       = inp$re_names,
    group_levels   = inp$group_levels,
    collect_block1 = TRUE,
    progbar        = progbar_use,
    stage_label    = "main",
    fixef_mode     = fixef_mode,
    b_mode         = ranef_mode,
    b_start        = ranef_mode,
    tau2_start     = tau2_start
  )

  staged <- .rLMM_format_sweep_out(
    v6_out       = raw,
    n            = inp$n,
    re_names     = inp$re_names,
    group_levels = inp$group_levels,
    fixef_mode   = fixef_mode,
    fixef_init   = fixef_mode
  )

  .lmebayes_assemble_reg_result(
    staged           = staged,
    call             = cl,
    m_convergence    = m_convergence,
    convergence_info = convergence_info,
    pfamily_list     = pfamily_list,
    dispprior_list   = prior_list_block1,
    family           = gaussian(),
    groupef.mode     = ranef_mode,
    any_non_normal   = any_non_normal,
    design           = design,
    result_class     = result_class,
    parent_class     = "rLMMNormal_reg",
    draw_engine      = draw_engine,
    sim_method_used  = "TWO_BLOCK_GIBBS",
    icm_info         = icm_info,
    offset              = inp$offset,
    weights             = if (is.null(inp$weights)) 1 else inp$weights
  )
}

#' Shared exact-iid sampling pipeline for \code{rLMMNormal_reg_known_vcov_iid}
#'
#' Counterpart to \code{.rLMMNormal_reg_run()} that draws directly from
#' the exact multivariate-normal posterior via
#' \code{\link{rLMMNormal_joint_iid}} instead of running ICM + Theorem~3
#' calibration + two-block Gibbs sweeps. No pilot stage, calibration, or
#' \code{icm_tol}/\code{icm_maxit} are needed: every draw is already exact,
#' so \code{m_convergence} is fixed at \code{1L}.
#' @noRd
.rLMMNormal_reg_run_iid <- function(
    inp,
    group,
    dispersion,
    pfamily_list,
    pf_summary,
    progbar,
    verbose,
    engine_label,
    result_class,
    cl
) {
  ## rLMMNormal_joint_iid() derives its own Block~1 P from pfamily_list
  ## internally (and rejects a caller-supplied P/Sigma).
  prior_list_block1 <- list(
    dispersion = dispersion,
    ddef       = FALSE
  )

  ## rLMMNormal_joint_iid() has no 'group_name' formal: attach it to
  ## 'group' itself so .lmebayes_resolve_group_name() picks it up there.
  attr(group, "group_name") <- inp$group_name
  out <- rLMMNormal_joint_iid(
    n                 = inp$n,
    y                 = inp$y,
    x                 = inp$D,
    group             = group,
    x_hyper           = inp$W,
    prior_list_block1 = prior_list_block1,
    pfamily_list      = pfamily_list,
    progbar           = progbar,
    verbose           = verbose
  )

  staged <- .rLMM_format_v2_out(
    v2_out       = out,
    n            = inp$n,
    re_names     = inp$re_names,
    group_levels = inp$group_levels,
    fixef_mode   = out$fixef_mean,
    fixef_init   = out$fixef_mean
  )

  .lmebayes_assemble_reg_result(
    staged           = staged,
    call             = cl,
    m_convergence    = 1L,
    convergence_info = list(
      method        = "exact_iid",
      m_convergence = 1L,
      m_min         = 0L,
      lambda_star   = 0,
      eigenvalues   = NULL
    ),
    pfamily_list     = pfamily_list,
    dispprior_list   = prior_list_block1,
    family           = gaussian(),
    groupef.mode     = out$b_mean,
    any_non_normal   = FALSE,
    design           = list(
      y             = inp$y,
      D             = inp$D,
      group         = factor(group, levels = inp$group_levels),
      W             = inp$W,
      groupef.names = inp$re_names,
      group_name    = inp$group_name
    ),
    result_class     = result_class,
    parent_class     = "rLMMNormal_reg",
    draw_engine      = "rLMMNormal_joint_iid",
    sim_method_used  = "DEFAULT",
    icm_info         = list(converged = TRUE, iterations = 1L, delta = 0),
    offset              = inp$offset,
    weights             = if (is.null(inp$weights)) 1 else inp$weights
  )
}

#' @describeIn rLMM_reg Dispatcher for fixed \eqn{\sigma^2}: routes to
#'   \code{\link{rLMMNormal_reg_known_vcov}} or
#'   \code{\link{rLMMNormal_reg_estimated_vcov}} by population \code{pfamily_list}.
#' @param sim_method Simulation method for \code{rLMMNormal_reg_known_vcov}:
#'   \code{"DEFAULT"} (exact iid draws, see
#'   \code{\link{rLMMNormal_reg_known_vcov_iid}}) or
#'   \code{"TWO_BLOCK_GIBBS"} (two-block Gibbs sweeps, see
#'   \code{\link{rLMMNormal_reg_known_vcov_two_bg}}). Accepted but has no
#'   effect on \code{rLMMNormal_reg_estimated_vcov} (only \code{TWO_BLOCK_GIBBS}
#'   is implemented there): kept as a formal purely so callers may pass the
#'   same argument to either route without branching.
#' @export
rLMMNormal_reg <- function(
    n,
    y,
    D,
    group,
    W,
    pfamily_list,
    dispprior_list,
    offset          = NULL,
    weights         = 1,
    icm_tol         = 1e-10,
    icm_maxit       = 200L,
    tv_tol          = 0.01,
    progbar         = TRUE,
    verbose         = FALSE,
    sim_method      = "DEFAULT"
) {
  cl <- match.call()

  group_name <- .lmebayes_resolve_group_name(
    group, substitute(group), fn_name = "rLMMNormal_reg"
  )

  inp <- .rLMM_validate_matrix_inputs(
    n, y, D, W, tv_tol,
    group_name, group
  )
  inp$offset <- offset
  inp$weights <- weights
  dispersion <- .rLMM_validate_fixed_dispersion_prior_list(
    dispprior_list, group_levels = inp$group_levels
  )
  pfamily_list <- .two_block_validate_pfamily_list(
    pfamily_list, inp$re_names, J = length(inp$group_levels)
  )
  pf_summary <- .two_block_summarize_pfamily_list(pfamily_list)

  route_fn <- if (pf_summary$all_dNormal) {
    rLMMNormal_reg_known_vcov
  } else {
    rLMMNormal_reg_estimated_vcov
  }
  mc <- match.call(expand.dots = FALSE)
  mc[[1L]] <- route_fn
  out <- eval(mc, parent.frame())
  out$call <- cl
  out
}

#' @describeIn rLMM_reg Fixed \eqn{\sigma^2}; all population components
#'   \code{dNormal} (known \eqn{\tau^2_k}). Dispatches on \code{sim_method}:
#'   \code{"DEFAULT"} draws directly from the exact multivariate-normal
#'   posterior via \code{\link{rLMMNormal_reg_known_vcov_iid}} (no Gibbs
#'   sweeps, no burn-in, no autocorrelation between draws);
#'   \code{"TWO_BLOCK_GIBBS"} runs
#'   \code{\link{rLMMNormal_reg_known_vcov_two_bg}} instead (two-block Gibbs
#'   with Theorem~3 rate calibration -- the only option available for
#'   \code{rLMMNormal_reg_estimated_vcov} and the GLMM routes).
#' @export
rLMMNormal_reg_known_vcov <- function(
    n,
    y,
    D,
    group,
    W,
    pfamily_list,
    dispprior_list,
    offset          = NULL,
    weights         = 1,
    icm_tol         = 1e-10,
    icm_maxit       = 200L,
    tv_tol          = 0.01,
    progbar         = TRUE,
    verbose         = FALSE,
    sim_method      = "DEFAULT"
) {
  cl <- match.call()
  fn_name <- "rLMMNormal_reg_known_vcov"
  sim_method <- .rLMM_validate_sim_method(sim_method, fn_name = fn_name)

  route_fn <- if (identical(sim_method, "TWO_BLOCK_GIBBS")) {
    rLMMNormal_reg_known_vcov_two_bg
  } else {
    rLMMNormal_reg_known_vcov_iid
  }
  mc <- match.call(expand.dots = FALSE)
  mc[[1L]] <- route_fn
  mc$sim_method <- NULL
  out <- eval(mc, parent.frame())
  out$call <- cl
  out
}

#' @describeIn rLMM_reg Fixed \eqn{\sigma^2}; all population components
#'   \code{dNormal}; exact iid draws from the closed-form multivariate-normal
#'   posterior via \code{\link{rLMMNormal_joint_iid}}
#'   (\code{sim_method = "DEFAULT"} route of
#'   \code{\link{rLMMNormal_reg_known_vcov}}).
#' @export
rLMMNormal_reg_known_vcov_iid <- function(
    n,
    y,
    D,
    group,
    W,
    pfamily_list,
    dispprior_list,
    offset          = NULL,
    weights         = 1,
    icm_tol         = 1e-10,
    icm_maxit       = 200L,
    tv_tol          = 0.01,
    progbar         = TRUE,
    verbose         = FALSE
) {
  cl <- match.call()
  fn_name <- "rLMMNormal_reg_known_vcov_iid"

  group_name <- .lmebayes_resolve_group_name(
    group, substitute(group), fn_name = fn_name
  )

  inp <- .rLMM_validate_matrix_inputs(
    n, y, D, W, tv_tol,
    group_name, group
  )
  inp$offset <- offset
  inp$weights <- weights
  dispersion <- .rLMM_validate_fixed_dispersion_prior_list(
    dispprior_list, fn_name = fn_name, group_levels = inp$group_levels
  )
  pfamily_list <- .two_block_validate_pfamily_list(
    pfamily_list, inp$re_names, J = length(inp$group_levels)
  )
  pf_summary <- .two_block_summarize_pfamily_list(pfamily_list)
  if (!pf_summary$all_dNormal) {
    stop(
      fn_name, "(): all population components must be dNormal(); ",
      "use rLMMNormal_reg_estimated_vcov() or rLMMNormal_reg().",
      call. = FALSE
    )
  }

  .rLMMNormal_reg_run_iid(
    inp              = inp,
    group            = group,
    dispersion       = dispersion,
    pfamily_list     = pfamily_list,
    pf_summary       = pf_summary,
    progbar          = progbar,
    verbose          = verbose,
    engine_label     = fn_name,
    result_class     = "rLMMNormal_reg_known_vcov",
    cl               = cl
  )
}

#' @describeIn rLMM_reg Fixed \eqn{\sigma^2}; all population components
#'   \code{dNormal} (known \eqn{\tau^2_k}). Exact Theorem~3 rate calibration;
#'   two-block Gibbs sweeps (\code{sim_method = "TWO_BLOCK_GIBBS"} route of
#'   \code{\link{rLMMNormal_reg_known_vcov}}).
#' @export
rLMMNormal_reg_known_vcov_two_bg <- function(
    n,
    y,
    D,
    group,
    W,
    pfamily_list,
    dispprior_list,
    offset          = NULL,
    weights         = 1,
    icm_tol         = 1e-10,
    icm_maxit       = 200L,
    tv_tol          = 0.01,
    progbar         = TRUE,
    verbose         = FALSE
) {
  cl <- match.call()
  fn_name <- "rLMMNormal_reg_known_vcov_two_bg"

  group_name <- .lmebayes_resolve_group_name(
    group, substitute(group), fn_name = fn_name
  )

  inp <- .rLMM_validate_matrix_inputs(
    n, y, D, W, tv_tol,
    group_name, group
  )
  inp$offset <- offset
  inp$weights <- weights
  dispersion <- .rLMM_validate_fixed_dispersion_prior_list(
    dispprior_list, fn_name = fn_name, group_levels = inp$group_levels
  )
  pfamily_list <- .two_block_validate_pfamily_list(
    pfamily_list, inp$re_names, J = length(inp$group_levels)
  )
  P <- .rLMM_P_from_pfamily_list(pfamily_list, inp$re_names)
  pf_summary <- .two_block_summarize_pfamily_list(pfamily_list)
  if (!pf_summary$all_dNormal) {
    stop(
      fn_name, "(): all population components must be dNormal(); ",
      "use rLMMNormal_reg_estimated_vcov() or rLMMNormal_reg().",
      call. = FALSE
    )
  }

  .rLMMNormal_reg_run(
    inp              = inp,
    group            = group,
    P                = P,
    dispersion       = dispersion,
    pfamily_list     = pfamily_list,
    pf_summary       = pf_summary,
    icm_tol          = icm_tol,
    icm_maxit        = icm_maxit,
    progbar          = progbar,
    verbose          = verbose,
    engine_label     = fn_name,
    any_non_normal   = FALSE,
    draw_engine      = "two_block_rNormal_reg_known_vcov",
    result_class     = "rLMMNormal_reg_known_vcov",
    cl               = cl
  )
}

#' @describeIn rLMM_reg Fixed \eqn{\sigma^2}; ING population components
#'   (estimated \eqn{\tau^2_k}). Optional pilot stage; conservative
#'   \code{disp_lower} rate bound. \code{sim_method} is accepted but has no
#'   effect: variance components are estimated here, so the joint posterior
#'   is not exactly Gaussian and only two-block Gibbs sampling is
#'   implemented.
#' @export
rLMMNormal_reg_estimated_vcov <- function(
    n,
    y,
    D,
    group,
    W,
    pfamily_list,
    dispprior_list,
    offset          = NULL,
    weights         = 1,
    icm_tol         = 1e-10,
    icm_maxit       = 200L,
    tv_tol          = 0.01,
    progbar         = TRUE,
    verbose         = FALSE,
    gap_tol         = 0.0196,
    mode_gap_max    = 1.0,
    diag_sweeps     = FALSE,
    stage_verbose   = FALSE,
    sim_method      = "DEFAULT"
) {
  cl <- match.call()
  fn_name <- "rLMMNormal_reg_estimated_vcov"
  .rLMM_validate_sim_method(sim_method, fn_name = fn_name)

  group_name <- .lmebayes_resolve_group_name(
    group, substitute(group), fn_name = fn_name
  )

  inp <- .rLMM_validate_matrix_inputs(
    n, y, D, W, tv_tol,
    group_name, group
  )
  inp$offset <- offset
  inp$weights <- weights
  dispersion <- .rLMM_validate_fixed_dispersion_prior_list(
    dispprior_list, fn_name = fn_name, group_levels = inp$group_levels
  )
  pfamily_list <- .two_block_validate_pfamily_list(
    pfamily_list, inp$re_names, J = length(inp$group_levels)
  )
  P <- .rLMM_P_from_pfamily_list(pfamily_list, inp$re_names)
  pf_summary <- .two_block_summarize_pfamily_list(pfamily_list)
  if (pf_summary$all_dNormal) {
    stop(
      fn_name, "(): at least one population component must not be dNormal(); ",
      "use rLMMNormal_reg_known_vcov() or rLMMNormal_reg().",
      call. = FALSE
    )
  }

  .rLMMNormal_reg_run_with_pilot(
    inp            = inp,
    group          = group,
    P              = P,
    dispersion     = dispersion,
    pfamily_list   = pfamily_list,
    pf_summary     = pf_summary,
    icm_tol        = icm_tol,
    icm_maxit      = icm_maxit,
    progbar        = progbar,
    verbose        = verbose,
    stage_verbose  = stage_verbose,
    gap_tol        = gap_tol,
    mode_gap_max   = mode_gap_max,
    diag_sweeps    = diag_sweeps,
    engine_label   = fn_name,
    result_class   = "rLMMNormal_reg_estimated_vcov",
    cl             = cl
  )
}

#' @describeIn rLMM_reg Legacy outer-loop sampler: draws \eqn{\sigma^2} via
#'   \code{\link[glmbayesCore]{rGamma_reg}}, then calls \code{\link{rLMMNormal_reg}} each
#'   replicate. Prefer \code{\link{rLMMindepNormalGamma_reg_known_vcov}} or
#'   \code{\link{rLMMindepNormalGamma_reg_estimated_vcov}} for the ING
#'   observation-dispersion path used by \code{\link{rlmerb}}.
#' @export
rLMMindepNormalGamma_reg <- function(
    n,
    y,
    D,
    group,
    W,
    pfamily_list,
    dispprior_list,
    offset          = NULL,
    weights         = 1,
    icm_tol         = 1e-10,
    icm_maxit       = 200L,
    tv_tol          = 0.01,
    progbar         = TRUE,
    verbose         = FALSE
) {
  cl <- match.call()

  group_name <- .lmebayes_resolve_group_name(
    group, substitute(group), fn_name = "rLMMindepNormalGamma_reg"
  )

  inp <- .rLMM_validate_matrix_inputs(
    n, y, D, W, tv_tol,
    group_name, group
  )
  inp$offset <- offset
  inp$weights <- weights
  dispprior_list <- .rLMM_validate_dGamma_dispersion_prior_list(dispprior_list)
  dispersion_fix <- .rLMM_dispersion_fix_from_prior_list(
    dispprior_list, fn_name = "rLMMindepNormalGamma_reg"
  )

  pfamily_list <- .two_block_validate_pfamily_list(
    pfamily_list, inp$re_names, J = length(inp$group_levels)
  )
  P <- .rLMM_P_from_pfamily_list(pfamily_list, inp$re_names)
  pf_summary <- .two_block_summarize_pfamily_list(pfamily_list)

  prior_list_block1_cal <- list(
    P          = P,
    dispersion = dispersion_fix,
    ddef       = FALSE
  )
  .two_block_validate_block1_prior(prior_list_block1_cal, family = gaussian())

  icm <- .rLMM_icm_at_start(
    y                     = inp$y,
    D                     = inp$D,
    group                 = group,
    W                     = inp$W,
    prior_list_block1     = prior_list_block1_cal,
    pfamily_list          = pfamily_list,
    re_names              = inp$re_names,
    group_levels          = inp$group_levels,
    group_name            = inp$group_name,
    icm_tol               = icm_tol,
    icm_maxit             = icm_maxit,
    verbose               = verbose,
    engine_label          = "rLMMindepNormalGamma_reg"
  )
  fixef_start <- icm$start
  ranef_mode  <- icm$b_start
  icm_info    <- icm$icm

  calib <- .rLMM_calibrate_m_convergence(
    D                 = inp$D,
    group             = group,
    W                 = inp$W,
    prior_list_block1 = prior_list_block1_cal,
    pfamily_list      = pfamily_list,
    group_levels      = inp$group_levels,
    tv_tol            = inp$tv_tol,
    any_non_normal    = pf_summary$any_non_normal,
    engine_label      = "rLMMindepNormalGamma_reg",
    verbose           = verbose
  )
  m_convergence    <- calib$m_convergence
  convergence_info <- calib$convergence_info
  convergence_info$draw_engine <- "rLMMindepNormalGamma_reg_outer"

  if (is.null(ranef_mode)) {
    ranef_mode <- matrix(
      0,
      nrow = length(inp$group_levels),
      ncol = length(inp$re_names),
      dimnames = list(inp$group_levels, inp$re_names)
    )
  }
  b_mat <- ranef_mode

  n_obs <- length(inp$y)
  x_disp <- matrix(1, nrow = n_obs, ncol = 1L)
  wt     <- rep(1, n_obs)

  fixef_cur <- fixef_start
  re_names  <- inp$re_names
  n         <- inp$n

  fixef_draws <- stats::setNames(
    lapply(re_names, function(k) {
      q_k <- length(fixef_cur[[k]])
      matrix(NA_real_, nrow = n, ncol = q_k,
             dimnames = list(NULL, names(fixef_cur[[k]])))
    }),
    re_names
  )
  dispersion_ranef_draws <- numeric(n)
  coef_parts             <- vector("list", n)
  disp_fixef_draws       <- NULL
  iters_fixef_draws      <- NULL

  if (isTRUE(progbar) && n > 1L) {
    pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  for (i in seq_len(n)) {
    mu_hat <- .rLMM_observation_mu(
      inp$D, group, b_mat, inp$group_levels
    )
    gamma_out <- glmbayesCore::rGamma_reg(
      n           = 1L,
      y           = inp$y,
      x           = x_disp,
      prior_list  = dispprior_list,
      offset      = mu_hat,
      weights     = wt,
      family      = gaussian(),
      progbar     = FALSE,
      verbose     = FALSE
    )
    sigma2_i <- as.numeric(gamma_out$dispersion[1L])
    dispersion_ranef_draws[i] <- sigma2_i

    draw_engine <- if (pf_summary$all_dNormal) {
      "two_block_rNormal_reg_known_vcov"
    } else {
      "two_block_rNormal_reg_estimated_vcov"
    }
    inp_i <- inp
    inp_i$n <- 1L
    lmm_i <- .rLMMNormal_reg_run(
      inp            = inp_i,
      group          = group,
      P              = P,
      dispersion     = sigma2_i,
      pfamily_list   = pfamily_list,
      pf_summary     = pf_summary,
      icm_tol        = icm_tol,
      icm_maxit      = icm_maxit,
      progbar        = FALSE,
      verbose        = FALSE,
      engine_label   = "rLMMindepNormalGamma_reg",
      any_non_normal = pf_summary$any_non_normal,
      draw_engine    = draw_engine,
      result_class   = "rLMMNormal_reg",
      cl             = NULL,
      fixef_start    = fixef_cur
    )

    for (k in re_names) {
      fixef_draws[[k]][i, ] <- lmm_i$popef[[k]][1L, ]
    }
    coef_parts[[i]] <- lmm_i$groupef
    fixef_cur <- lapply(lmm_i$popef, function(m) {
      stats::setNames(m[1L, ], colnames(m))
    })
    b_mat <- .rLMM_b_matrix_from_coefficients(
      lmm_i$groupef,
      re_names,
      inp$group_levels,
      inp$group_name
    )
    if (i == 1L) {
      disp_fixef_draws  <- lmm_i$popef.dispersion
      iters_fixef_draws <- lmm_i$popef.iters
    } else {
      disp_fixef_draws  <- rbind(disp_fixef_draws, lmm_i$popef.dispersion)
      iters_fixef_draws <- rbind(iters_fixef_draws, lmm_i$popef.iters)
    }

    if (isTRUE(progbar) && n > 1L) {
      utils::setTxtProgressBar(pb, i / n)
    }
  }

  coefficients <- do.call(rbind, coef_parts)
  rownames(coefficients) <- NULL

  v2_out <- list(
    fixef_draws            = fixef_draws,
    coefficients           = coefficients,
    dispersion_fixef_draws = disp_fixef_draws,
    iters_fixef_draws      = iters_fixef_draws,
    mu_all_last            = build_mu_all(
      list(
        W             = inp$W,
        groupef.names = re_names,
        group         = factor(group, levels = inp$group_levels)
      ),
      fixef_cur,
      group_levels = inp$group_levels
    )$mu_all
  )

  staged <- .rLMM_format_v2_out(
    v2_out       = v2_out,
    n            = n,
    re_names     = re_names,
    group_levels = inp$group_levels,
    fixef_mode   = fixef_start,
    fixef_init   = fixef_start
  )
  staged[["group.dispersion"]] <- dispersion_ranef_draws
  staged[["group.dispersion.mean"]] <- mean(dispersion_ranef_draws)

  .lmebayes_assemble_reg_result(
    staged           = staged,
    call             = cl,
    m_convergence    = m_convergence,
    convergence_info = convergence_info,
    pfamily_list     = pfamily_list,
    dispprior_list   = dispprior_list,
    family           = gaussian(),
    groupef.mode     = ranef_mode,
    any_non_normal   = pf_summary$any_non_normal,
    design           = list(
      y             = inp$y,
      D             = inp$D,
      group         = factor(group, levels = inp$group_levels),
      W             = inp$W,
      groupef.names = re_names,
      group_name    = inp$group_name
    ),
    result_class     = "rLMMindepNormalGamma_reg",
    parent_class     = "rLMMindepNormalGamma_reg",
    draw_engine      = "rLMMindepNormalGamma_reg_outer",
    sim_method_used  = "TWO_BLOCK_GIBBS",
    icm_info         = icm_info,
    offset              = inp$offset,
    weights             = if (is.null(inp$weights)) 1 else inp$weights
  )
}

#' @describeIn rLMM_reg Random \eqn{\sigma^2} (pooled or per-group ING
#'   observation prior, drawn jointly with group coefficients); all population
#'   components \code{dNormal}. Used by \code{\link{rlmerb}} when
#'   \code{group.dispersion} is \code{dGamma()} / \code{dGamma_list()} and
#'   population priors are all \code{dNormal}.
#' @export
rLMMindepNormalGamma_reg_known_vcov <- function(
    n,
    y,
    D,
    group,
    W,
    pfamily_list,
    dispprior_list,
    offset          = NULL,
    weights         = 1,
    icm_tol         = 1e-10,
    icm_maxit       = 200L,
    tv_tol          = 0.01,
    progbar         = TRUE,
    verbose         = FALSE
) {
  cl <- match.call()
  fn_name <- "rLMMindepNormalGamma_reg_known_vcov"

  group_name <- .lmebayes_resolve_group_name(
    group, substitute(group), fn_name = fn_name
  )

  inp <- .rLMM_validate_matrix_inputs(
    n, y, D, W, tv_tol,
    group_name, group
  )
  inp$offset <- offset
  inp$weights <- weights
  ing_prior_list <- .rLMM_validate_ing_measurement_prior_list(
    dispprior_list, length(inp$re_names), fn_name = fn_name,
    group_levels = inp$group_levels
  )

  pfamily_list <- .two_block_validate_pfamily_list(
    pfamily_list, inp$re_names, J = length(inp$group_levels)
  )
  P <- .rLMM_P_from_pfamily_list(pfamily_list, inp$re_names)
  ## The shared random-effects covariance every group's beta_j prior uses is
  ## always Psi = solve(P) -- P already derived above from pfamily_list --
  ## never a value read from 'dispprior_list' (see
  ## .rLMM_validate_ing_measurement_prior_list()'s header comment).
  ing_prior_list$Sigma <- solve(P)
  pf_summary <- .two_block_summarize_pfamily_list(pfamily_list)
  if (!pf_summary$all_dNormal) {
    stop(
      fn_name, "(): all population components must be dNormal(); ",
      "use rLMMindepNormalGamma_reg_estimated_vcov().",
      call. = FALSE
    )
  }

  .rLMMIngNormal_reg_run_with_pilot(
    inp                = inp,
    group              = group,
    P                  = P,
    ing_prior_list     = ing_prior_list,
    pfamily_list       = pfamily_list,
    pf_summary         = pf_summary,
    icm_tol            = icm_tol,
    icm_maxit          = icm_maxit,
    progbar            = progbar,
    verbose            = verbose,
    any_non_normal     = FALSE,
    random_measurement = TRUE,
    engine_label       = fn_name,
    result_class       = "rLMMindepNormalGamma_reg_known_vcov",
    cl                 = cl
  )
}

#' @describeIn rLMM_reg Random \eqn{\sigma^2} (pooled or per-group ING
#'   observation prior, drawn jointly with group coefficients); ING population
#'   components. Pilot/UB calibration via the sweep-outer sampler. Default path
#'   for \code{\link{rlmerb}} with \code{dGamma()} / \code{dGamma_list()}
#'   \code{group.dispersion} and ING population priors.
#' @export
rLMMindepNormalGamma_reg_estimated_vcov <- function(
    n,
    y,
    D,
    group,
    W,
    pfamily_list,
    dispprior_list,
    offset          = NULL,
    weights         = 1,
    icm_tol         = 1e-10,
    icm_maxit       = 200L,
    tv_tol          = 0.01,
    progbar         = TRUE,
    verbose         = FALSE,
    gap_tol         = 0.0196,
    mode_gap_max    = 1.0,
    diag_sweeps     = FALSE,
    stage_verbose   = FALSE
) {
  cl <- match.call()
  fn_name <- "rLMMindepNormalGamma_reg_estimated_vcov"

  group_name <- .lmebayes_resolve_group_name(
    group, substitute(group), fn_name = fn_name
  )

  inp <- .rLMM_validate_matrix_inputs(
    n, y, D, W, tv_tol,
    group_name, group
  )
  inp$offset <- offset
  inp$weights <- weights
  ing_prior_list <- .rLMM_validate_ing_measurement_prior_list(
    dispprior_list, length(inp$re_names), fn_name = fn_name,
    group_levels = inp$group_levels
  )

  pfamily_list <- .two_block_validate_pfamily_list(
    pfamily_list, inp$re_names, J = length(inp$group_levels)
  )
  P <- .rLMM_P_from_pfamily_list(pfamily_list, inp$re_names)
  ## The shared random-effects covariance every group's beta_j prior uses is
  ## always Psi = solve(P) -- P already derived above from pfamily_list --
  ## never a value read from 'dispprior_list' (see
  ## .rLMM_validate_ing_measurement_prior_list()'s header comment).
  ing_prior_list$Sigma <- solve(P)
  pf_summary <- .two_block_summarize_pfamily_list(pfamily_list)
  if (pf_summary$all_dNormal) {
    stop(
      fn_name, "(): at least one population component must not be dNormal(); ",
      "use rLMMindepNormalGamma_reg_known_vcov().",
      call. = FALSE
    )
  }

  .rLMMIngNormal_reg_run_with_pilot(
    inp                = inp,
    group              = group,
    P                  = P,
    ing_prior_list     = ing_prior_list,
    pfamily_list       = pfamily_list,
    pf_summary         = pf_summary,
    icm_tol            = icm_tol,
    icm_maxit          = icm_maxit,
    progbar            = progbar,
    verbose            = verbose,
    stage_verbose      = stage_verbose,
    gap_tol            = gap_tol,
    mode_gap_max       = mode_gap_max,
    diag_sweeps        = diag_sweeps,
    any_non_normal     = TRUE,
    random_measurement = TRUE,
    engine_label       = fn_name,
    result_class       = "rLMMindepNormalGamma_reg_estimated_vcov",
    cl                 = cl
  )
}

#' Format v2 batch output for staged \code{fixef.*} naming
#' @noRd
.rLMM_format_v2_out <- function(
    v2_out,
    n,
    re_names,
    group_levels,
    fixef_mode,
    fixef_init
) {
  x <- list(
    fixef_draws            = v2_out$fixef_draws,
    coefficients           = v2_out$coefficients,
    dispersion_fixef_draws = v2_out$dispersion_fixef_draws,
    iters_fixef_draws      = v2_out$iters_fixef_draws,
    iters_ranef_draws      = v2_out$iters_ranef_draws,
    n                      = n
  )
  .two_block_as_staged_names(
    x,
    fixef_mode = fixef_mode,
    fixef_init = fixef_init
  )
}

