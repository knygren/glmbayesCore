#' Model setup for generalized linear mixed models
#'
#' Wrapper around \code{\link[lme4]{lmer}} or \code{\link[lme4]{glmer}} for
#' models with exactly one grouping factor. Design matrices come from
#' \code{formula} (including cross-level RE moderation terms). Random-effect
#' and residual variance components (\code{Psi}, \code{dispersion}) come
#' from the same reference \code{lmer}/\code{glmer} fit on \code{formula}, so
#' printed summaries match \code{summary(lmer(...))} on that formula.
#'
#' \strong{Uncorrelated random effects (\code{||}).}
#' The sampler treats \code{Sigma_ranef} as diagonal (no off-diagonal
#' covariance). Multi-coefficient random terms must use \code{||}, e.g.
#' \code{(1 + x || group)} rather than \code{(1 + x | group)}. A single
#' random intercept may use \code{(1 | group)}; \code{(1 || group)} is not
#' supported by \code{lme4}.
#'
#' @details
#' \strong{Fixed-effect constraints.}
#' \code{model_setup} accepts the same formula language as
#' \code{\link[lme4]{lmer}}, subject to one structural rule: every fixed
#' effect that does \emph{not} correspond to a random-slope term must be a
#' \emph{group-constant} (level-2) covariate---a predictor whose value is the
#' same for every observation within a given group.  School-level attributes
#' such as \code{private_school} or \code{title1} satisfy this constraint.
#' Student-level covariates that vary \emph{within} groups may appear as fixed
#' main effects only when they \emph{also} appear as random slopes (they then
#' represent the population mean slope \eqn{\gamma_{10}}, e.g.,
#' \code{distracted_ppvt}).  Cross-level interactions of the form
#' \code{level2_var:random_slope} (e.g.,
#' \code{free_reduced_lunch:distracted_a1}) are additionally permitted; they
#' moderate the prior mean of the corresponding random slope across groups (see
#' \code{\link{extract_re_hyper_matrices}}).  Fixed terms that are none of
#' these three types---level-2 covariate, population mean slope, or cross-level
#' moderation interaction---are rejected with an informative error.
#'
#' \strong{Two-step identifiability assessment.}
#' After fitting \code{lmer}/\code{glmer}, \code{model_setup} calls
#' \code{\link{check_identifiability}} to assess whether the model is
#' empirically identified at both the within-group (\code{groupef.rank},
#' \code{groupef.estimable}, \code{groupef.glm_check}) and across-group
#' (\code{popef.rank}, \code{popef.deficient}, \code{popef.rank_ok}) levels,
#' and copies its fields onto the returned design object -- see
#' \code{\link{check_identifiability}} for the full algorithm description
#' (including the family-specific estimability checks and how
#' \code{popef.rank_ok} is derived). Non-estimable groups are flagged but
#' retained in the \code{lmer}/\code{glmer} fit;
#' \code{\link[lmebayesCore:Prior_Setup_GLMM]{Prior_Setup_GLMM}} excludes them when calibrating priors,
#' and requires \code{popef.rank_ok = TRUE} to derive default hyperpriors
#' automatically.
#'
#' The example uses \code{big_word_club} from \pkg{bayesrules}
#' (see \code{?bayesrules::big_word_club}) and the same
#' formula as the full \code{lmerb()} demo in lmebayes
#' (\code{demo("Ex_12_lmerb_BigWordClub", package = "lmebayes")}).
#'
#' @param formula A two-sided linear formula, as in \code{\link[lme4]{lmer}}:
#'   response on the left of \code{~}, fixed and random-effects terms
#'   (\code{(terms | group)}) on the right. Used for both design extraction
#'   and the reference \code{lmer}/\code{glmer} fit (fixed effects, hyper
#'   calibration, and variance components); see Details for the
#'   single-grouping-factor / \code{||} constraints.
#' @param data An optional data frame containing the variables named in
#'   \code{formula} (as in \code{\link[lme4]{lmer}}); if omitted, variables
#'   are taken from the environment of \code{formula}.
#' @param family As in \code{\link[stats]{glm}}: a character string naming a
#'   family function, a family function, or the result of a call to one.
#'   Defaults to \code{gaussian()}, in which case \code{\link[lme4]{lmer}} is
#'   used. Non-Gaussian families use \code{\link[lme4]{glmer}}.
#' @param REML Logical scalar, as in \code{lmer}'s own \code{REML} argument:
#'   should the reference fit's estimates optimize the REML criterion rather
#'   than the log-likelihood? Passed to \code{lmer} when
#'   \code{family = gaussian()}; ignored otherwise (\code{glmer} has no REML
#'   criterion).
#' @param control \code{\link[lme4]{lmerControl}} when
#'   \code{family = gaussian()}, otherwise \code{\link[lme4]{glmerControl}} --
#'   including the nonlinear optimizer and its settings (as in
#'   \code{lmer}/\code{glmer}'s own \code{control} argument); passed through
#'   to the reference fit when \code{fit_mer = TRUE}. \code{NULL} (default)
#'   resolves to \code{lmerControl()}'s defaults for Gaussian families, or
#'   \code{glmer}'s own default control otherwise.
#' @param fit_mer If \code{TRUE} (default), fit the reference
#'   \code{lmer}/\code{glmer} model and extract variance components --
#'   analogous to \code{glmmTMB}'s \code{doFit} argument. If \code{FALSE},
#'   skip the reference fit and return design matrices and rank diagnostics
#'   only (used by \code{glmerb()} in \pkg{lmebayes}, which fits its own
#'   \code{glmer} separately).
#' @param dispformula As in \code{glmmTMB}'s \code{dispformula} argument: a
#'   one-sided formula selecting the measurement-dispersion structure:
#'   \code{~1} (default, pooled) or \code{~<group_name>} (matching the
#'   random-effects grouping factor exactly, requesting per-group
#'   dispersion). Only affects \code{design$glmmTMB_fit} (see \code{@return});
#'   \code{design$lmer}/\code{design$glmer} are always the plain
#'   pooled-dispersion \code{lmer}/\code{glmer} fit, regardless of
#'   \code{dispformula}. \code{~<group_name>} additionally requires
#'   \code{family = gaussian()} (no observation-level dispersion parameter
#'   to model per group otherwise). Mirrors the \code{dispformula} argument
#'   on \code{\link[lmebayesCore:Prior_Setup_GLMM]{Prior_Setup_GLMM}} and on
#'   \code{lmerb()}/\code{glmerb()} in \pkg{lmebayes}; all three are
#'   independent arguments that must be kept consistent by the caller.
#' @param start Optional starting values, passed to the reference
#'   \code{lmer}/\code{glmer} fit's own \code{start} argument (numeric vector
#'   or named list with \code{par}/\code{theta}, and for \code{glmer} also
#'   \code{fixef}/\code{beta}; see their documentation). \code{NULL} (default)
#'   uses \code{lme4}'s own defaults.
#' @param verbose Integer scalar passed to the reference \code{lmer} (or
#'   \code{glmer} for non-Gaussian \code{family}); see their \code{verbose}
#'   argument (\code{> 0} prints optimizer iterations; \code{> 1} also prints
#'   inner PIRLS steps).
#' @param subset An optional expression indicating the subset of rows of
#'   \code{data} to use, as in \code{lmer}'s \code{subset} argument.  Passed
#'   to both \code{\link{extract_re_hyper_matrices}} (via
#'   \code{\link[lme4]{lFormula}}) and the reference \code{lmer}/\code{glmer}
#'   (and \code{glmmTMB} when used) so the sampler design matrices and the
#'   classical reference fit share the same rows.
#' @param weights An optional vector of prior weights, as in \code{lmer}'s
#'   \code{weights} argument.  Passed to design extraction
#'   (\code{lFormula}), the reference \code{lmer}/\code{glmer} fit, and
#'   \code{glmmTMB} when used, so dropped rows match \code{subset}/
#'   \code{na.action}.  Stored on the returned design as length-\code{n}
#'   \code{weights} (default all ones).  Mixed-model sampler engines that
#'   still hard-code unit weights do not consume this vector yet; see
#'   \code{\link[lmebayesCore:Prior_Setup_GLMM]{Prior_Setup_GLMM}}.
#' @param na.action A function indicating what should happen when \code{data}
#'   contain \code{NA}s (default \code{na.omit}, as in \code{lmer}).  Passed
#'   to both design extraction (\code{lFormula}) and the reference fit so
#'   dropped rows match.
#' @param offset An optional a priori known component to include in the
#'   linear predictor, as in \code{lmer}'s \code{offset} argument (plus any
#'   \code{offset()} term in \code{formula}).  Passed to design extraction
#'   and the reference fit; stored on the returned design as length-\code{n}
#'   \code{offset} (default all zeros).  Mixed-model sampler engines that
#'   still hard-code a zero offset do not consume this vector yet.
#' @param contrasts An optional list of contrasts (see the
#'   \code{contrasts.arg} of \code{model.matrix.default}, as in \code{lmer}).
#'   Passed to both design extraction (\code{lFormula}) and the reference
#'   fit so factor coding in \code{W}/\code{D} matches \code{fixef}/\code{vcov}.
#' @param devFunOnly Logical, as in \code{lmer}'s \code{devFunOnly} argument:
#'   if \code{TRUE}, return only the deviance evaluation function instead of
#'   a fitted model. Gaussian \code{lmer} fits only; ignored for \code{glmer}.
#' @param ... Additional arguments passed to \code{extract_re_hyper_matrices()}
#'   for design extraction and, when \code{fit_mer = TRUE}, to the reference
#'   \code{lmer}/\code{glmer} fit (as in \code{lm}/\code{glm}'s own
#'   \code{...}: further arguments passed to or from other methods).
#' @return Object of class \code{"model_setup"} with:
#'   \describe{
#'     \item{\code{y}}{Response vector, length \code{nrow(D)}.}
#'     \item{\code{weights}}{Numeric prior-weight vector of length
#'       \code{length(y)} (ones when unspecified), aligned with \code{y}/
#'       \code{D}/\code{group} after \code{subset}/\code{na.action}.}
#'     \item{\code{offset}}{Numeric offset vector of length \code{length(y)}
#'       (zeros when unspecified and no formula \code{offset()}), aligned
#'       with \code{y}/\code{D}/\code{group}.}
#'     \item{\code{D}}{Group-effects model matrix (\code{n_obs} x
#'       \code{p_re}): per-observation loadings on the within-group
#'       group-effect coefficient vector (columns \code{groupef.names}).}
#'     \item{\code{group}}{Factor of length \code{nrow(D)} for block
#'       subsetting; \code{levels(group)} fixes group order used elsewhere.
#'       Also carries \code{attr(group, "group_name")} (same value as
#'       \code{group_name}) so \code{group} alone can be passed to sampler
#'       engine functions without a separate \code{group_name} argument.}
#'     \item{\code{group_name}}{Grouping factor name.}
#'     \item{\code{groupef.names}}{Group-effect coefficient names from
#'       \code{lme4}; always identical to \code{names(W)} (and to
#'       \code{colnames(D)}).}
#'     \item{\code{W}}{Named list of group-level hyper-design matrices (one
#'       per \code{groupef.names} entry, one row per group level).}
#'     \item{\code{popef.moderation}}{Cross-level moderation metadata; see
#'       \code{\link{extract_re_hyper_matrices}}.}
#'     \item{\code{formula}}{The \code{formula} argument, as supplied.}
#'     \item{\code{family}}{The normalized \code{\link[stats]{family}}
#'       object.}
#'     \item{\code{call}}{The matched call.}
#'     \item{\code{terms}}{The \code{\link[stats]{terms}} object for
#'       \code{formula}, as in \code{lm}/\code{glm} (and
#'       \code{glmbayes::glmb}/\code{lmb}).}
#'     \item{\code{data}}{The \code{data} argument, as supplied (\code{NULL}
#'       when omitted), as in \code{lm}/\code{glm} (and
#'       \code{glmbayes::glmb}/\code{lmb}).}
#'     \item{\code{dispformula}}{The \code{dispformula} argument, as
#'       supplied.}
#'     \item{\code{lmer} / \code{glmer}}{The reference
#'       \code{\link[lme4]{lmer}} (Gaussian) or \code{\link[lme4]{glmer}}
#'       (otherwise) fit on \code{formula}; \code{NULL} when
#'       \code{fit_mer = FALSE}.}
#'     \item{\code{varcorr}}{The raw \code{\link[lme4]{VarCorr}} object from
#'       the reference fit.}
#'     \item{\code{Psi}}{Named numeric vector (names \code{groupef.names}) of
#'       group-effect variances.}
#'     \item{\code{dispersion}}{Residual (measurement) variance from the
#'       reference fit.}
#'     \item{\code{glmmTMB_fit}}{\code{NULL} when \code{dispformula = ~1} or
#'       \code{fit_mer = FALSE}; otherwise a \code{\link[glmmTMB]{glmmTMB}}
#'       fit on \code{formula} with the same \code{dispformula} (Gaussian
#'       models only), fit via the same \code{data}/\code{REML}/
#'       \code{control}-derived arguments as \code{lmer}. This is
#'       \strong{additive}: \code{lmer}/\code{glmer} are never replaced by a
#'       \code{glmmTMB} fit. \code{\link[lmebayesCore:Prior_Setup_GLMM]{Prior_Setup_GLMM}} and
#'       \code{lmerb()} (in \pkg{lmebayes}) reuse \code{glmmTMB_fit} as their
#'       per-group-dispersion calibration reference instead of fitting
#'       \code{glmmTMB} a second time.}
#'     \item{\code{groupef.rank}}{Named logical (names \code{levels(group)}):
#'       \code{TRUE} if that group's \eqn{D_j} is full column rank; see
#'       \code{\link{check_identifiability}}.}
#'     \item{\code{groupef.estimable}}{Named logical, same names as
#'       \code{groupef.rank}: \code{TRUE} if the group is additionally
#'       estimable (family-specific MLE/dispersion check); see
#'       \code{\link{check_identifiability}}.}
#'     \item{\code{groupef.glm_check}}{Per-group diagnostic data frame from
#'       \code{\link{check_identifiability}} (family-dependent; may be
#'       \code{NULL}).}
#'     \item{\code{popef.rank}}{Named logical (names \code{groupef.names}):
#'       \code{TRUE} if \code{W[[k]]} is full column rank when restricted to
#'       \code{groupef.estimable} groups.}
#'     \item{\code{popef.deficient}}{Negation of \code{popef.rank}.}
#'     \item{\code{popef.rank_ok}}{Scalar \code{TRUE} only when every
#'       \code{popef.rank} entry is \code{TRUE}; required by
#'       \code{\link[lmebayesCore:Prior_Setup_GLMM]{Prior_Setup_GLMM}} to derive default hyperpriors
#'       automatically.}
#'   }
#' @seealso \code{\link{check_identifiability}},
#'   \code{\link{extract_re_hyper_matrices}}
#' @example inst/examples/Ex_model_setup.R
#' @export
model_setup <- function(
    formula,
    data = NULL,
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
    devFunOnly = FALSE,
    fit_mer = TRUE,
    dispformula = ~1,
    ...
) {
  cl <- match.call()
  family <- .lmebayes_normalize_family(family)
  is_gaussian <- identical(family$family, "gaussian")
  if (is.null(control) && is_gaussian) {
    control <- lme4::lmerControl()
  }

  ## Shared model-frame args for design extraction and reference fits so
  ## design$y/D/group/W and lmer/glmer/glmmTMB see the same rows and
  ## factor coding (lFormula + lmer both accept these).
  ## Evaluate subset against data here (as model.frame does) so we can
  ## forward a concrete index through do.call(); passing the language
  ## object via do.call() would re-evaluate it in the wrong environment.
  ## weights/offset are force-evaluated (concrete vectors) and passed into
  ## lFormula so model.weights(fr)/model.offset(fr) match the design rows.
  frame_args <- list()
  if (!is.null(cl$subset)) {
    frame_args$subset <- eval(cl$subset, envir = data, enclos = parent.frame())
  }
  if (!missing(weights)) {
    frame_args$weights <- weights
  }
  if (!missing(offset)) {
    frame_args$offset <- offset
  }
  if (!missing(na.action)) {
    frame_args$na.action <- na.action
  }
  if (!is.null(contrasts)) {
    frame_args$contrasts <- contrasts
  }

  design <- do.call(
    extract_re_hyper_matrices,
    c(list(formula = formula, data = data), frame_args, list(...))
  )
  design$call    <- cl
  design$formula <- formula
  design$family  <- family
  design$terms   <- stats::terms(formula, data = data)
  design$data    <- data

  dispformula_kind <- .lmebayes_dispformula_kind(dispformula, design$group_name)
  if (identical(dispformula_kind, "group") && !is_gaussian) {
    stop(
      "'dispformula' must be ~1 for family = ", family$family,
      "() (no observation-level dispersion parameter for per-group ",
      "measurement-dispersion calibration).",
      call. = FALSE
    )
  }
  design$dispformula <- dispformula
  design$glmmTMB_fit <- NULL

  if (isTRUE(fit_mer)) {
    mer_args <- c(
      list(
        data = data,
        verbose = verbose
      ),
      if (!is.null(control)) list(control = control),
      .lmebayes_mer_optional_args(
        start = start,
        subset = if (!is.null(frame_args$subset)) frame_args$subset,
        weights = weights,
        na.action = na.action,
        offset = offset,
        contrasts = contrasts
      ),
      list(...)
    )

    if (is_gaussian) {
      fit_full <- do.call(
        lme4::lmer,
        c(list(formula = formula, REML = REML, devFunOnly = devFunOnly), mer_args)
      )
    } else {
      fit_full <- do.call(
        lme4::glmer,
        c(list(formula = formula, family = family), mer_args)
      )
    }

    if (lme4::isSingular(fit_full)) {
      message(
        if (is_gaussian) "lmer" else "glmer",
        " reference fit is singular -- check VarCorr; ",
        "RE variances may be on boundary."
      )
    }

    ## Guard: design extraction and the reference fit must agree on n.
    n_design <- length(design$y)
    n_fit <- length(stats::fitted(fit_full))
    if (n_design != n_fit) {
      stop(
        "Internal error: design extraction used ", n_design,
        " observation(s) but the reference ",
        if (is_gaussian) "lmer" else "glmer",
        " fit used ", n_fit,
        ". Pass the same subset/na.action/contrasts to both paths ",
        "(model_setup should do this automatically).",
        call. = FALSE
      )
    }

    vc <- extract_mer_variance_components(
      fit_full,
      design$groupef.names
    )
    if (is_gaussian) {
      design$lmer <- fit_full
    } else {
      design$glmer <- fit_full
    }
    design$varcorr <- vc$varcorr
    design$Psi <- vc$Psi
    design$dispersion <- vc$dispersion

    ## Additive: fit and store the glmmTMB per-group-dispersion reference
    ## alongside (never instead of) the lmer fit above, when 'dispformula'
    ## requests it. Prior_Setup_GLMM() and lmerb() (lmebayes) reuse this
    ## as their calibration reference instead of fitting glmmTMB a second
    ## time; see inst/DGAMMA_LIST_MARGINAL_AND_BOUNDS.md.
    if (is_gaussian && identical(dispformula_kind, "group")) {
      design$glmmTMB_fit <- do.call(
        .lmebayes_fit_glmmtmb_reference,
        c(
          list(
            formula     = formula,
            data        = data,
            family      = family,
            dispformula = dispformula,
            REML        = REML
          ),
          frame_args,
          if (!missing(weights)) list(weights = weights),
          if (!missing(offset)) list(offset = offset)
        )
      )
    }
  }

  # Two-step identifiability/estimability assessment (Level 1 rank +
  # estimability on D_j, Level 2 rank on W restricted to estimable
  # groups) -- see check_identifiability() for the full algorithm
  # description; model_setup() just copies its fields onto design.
  ident <- check_identifiability(
    y          = design$y,
    D          = design$D,
    group      = design$group,
    W          = design$W,
    family     = family,
    group_name = design$group_name
  )
  design$groupef.rank      <- ident$groupef.rank
  design$groupef.estimable <- ident$groupef.estimable
  design$groupef.glm_check <- ident$groupef.glm_check
  design$popef.rank        <- ident$popef.rank
  design$popef.deficient   <- ident$popef.deficient
  design$popef.rank_ok     <- ident$popef.rank_ok

  design
}

## Return character issue messages when an lme4 merMod fit failed checkConv
## or the inner optimizer (conv$opt != 0).  Empty character() = OK.
#' @noRd
.lmebayes_mer_convergence_issues <- function(fit, label = "reference fit") {
  if (is.null(fit)) {
    return(sprintf("%s: fit is NULL", label))
  }
  if (!inherits(fit, "merMod")) {
    return(sprintf("%s: not a merMod object", label))
  }
  issues <- character(0)
  conv   <- fit@optinfo$conv
  if (!is.null(conv$opt) && conv$opt != 0L) {
    issues <- c(
      issues,
      sprintf("%s: optimizer did not converge (conv$opt = %s)", label, conv$opt)
    )
  }
  lme4c <- conv$lme4
  if (!is.null(lme4c)) {
    code <- lme4c$code
    msgs <- lme4c$messages
    failed_code <- !is.null(code) && length(code) >= 1L &&
      !is.na(code[1L]) && code[1L] != 0L
    failed_msgs <- !is.null(msgs) && length(msgs) >= 1L
    if (failed_code || failed_msgs) {
      msg_txt <- if (failed_msgs) {
        paste(
          vapply(msgs, function(m) gsub("\\s+", " ", m), character(1L)),
          collapse = "; "
        )
      } else {
        sprintf("lme4 convergence code %s", code[1L])
      }
      issues <- c(issues, sprintf("%s: %s", label, msg_txt))
    }
  }
  issues
}

#' @noRd
.lmebayes_mer_optional_args <- function(
    start,
    subset,
    weights,
    na.action,
    offset,
    contrasts
) {
  args <- list()
  if (!missing(start) && !is.null(start)) {
    args$start <- start
  }
  if (!missing(subset)) {
    args$subset <- subset
  }
  if (!missing(weights)) {
    args$weights <- weights
  }
  if (!missing(na.action)) {
    args$na.action <- na.action
  }
  if (!missing(offset)) {
    args$offset <- offset
  }
  if (!missing(contrasts)) {
    args$contrasts <- contrasts
  }
  args
}

#' Phase~1 gate: mixed-model engines that still hard-code unit weights / zero
#' offset should call this before ignoring \code{design$weights}/
#' \code{design$offset}.
#' @noRd
.lmebayes_stop_if_nondefault_weights_offset <- function(
    weights,
    offset,
    where = "sampler"
) {
  n <- length(weights)
  if (n < 1L) {
    return(invisible(NULL))
  }
  if (!isTRUE(all.equal(as.numeric(weights), rep(1, n), tolerance = 0))) {
    stop(
      where, ": non-unit 'weights' are stored on the design but not yet ",
      "consumed by this mixed-model sampling path (Phase 1 storage only). ",
      "Use unit weights, or wait until sampler engines are wired.",
      call. = FALSE
    )
  }
  if (!isTRUE(all.equal(as.numeric(offset), rep(0, n), tolerance = 0))) {
    stop(
      where, ": non-zero 'offset' is stored on the design but not yet ",
      "consumed by this mixed-model sampling path (Phase 1 storage only). ",
      "Use a zero offset (or fold it into the formula once engines are wired).",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' @rdname model_setup
#' @method print model_setup
#' @param x A \code{model_setup} object.
#' @param ... Ignored.
#' @export
print.model_setup <- function(x, ...) {

  resp     <- deparse(x$formula[[2L]])
  re_names <- x$groupef.names
  grp      <- x$group_name
  n_obs    <- length(x$y)
  n_lev    <- nlevels(x$group)

  # ---- Call ------------------------------------------------------------------
  if (!is.null(x$call)) {
    cat("Call:\n  ", deparse1(x$call), "\n\n", sep = "")
  }

  # ---- Section 1: Level 1 (within-group) ------------------------------------
  cat("--- Level 1 (Within-Group) ---\n")
  cat(sprintf("  %s ~ %s\n\n", resp, paste(re_names, collapse = " + ")))
  cat(sprintf("  Observations             : %d\n", n_obs))
  cat(sprintf("  Group-effect coefficients: %d\n", length(re_names)))
  cat(sprintf("  Group                    : %s  [%d levels]\n", grp, n_lev))
  if (!is.null(x$groupef.rank)) {
    n_full <- sum(x$groupef.rank)
    cat(sprintf("  Full-rank D_j: %d of %d groups\n", n_full, n_lev))
    if (n_full < n_lev) {
      deficient <- names(x$groupef.rank)[!x$groupef.rank]
      shown     <- deficient[seq_len(min(10L, length(deficient)))]
      suffix    <- if (length(deficient) > 10L)
        sprintf(", ... (%d more)", length(deficient) - 10L) else ""
      cat(sprintf("    rank-deficient: %s%s\n",
                  paste(shown, collapse = ", "), suffix))
    }
  }
  if (!is.null(x$groupef.glm_check) && !is.null(x$groupef.estimable)) {
    n_est <- sum(x$groupef.estimable[x$groupef.rank])
    cat(sprintf(
      "  Full-rank & estimable: %d of %d full-rank (%d of %d total %s)\n",
      n_est, n_full, n_est, n_lev, grp
    ))
    if (n_est < n_full) {
      not_est <- names(x$groupef.estimable)[x$groupef.rank & !x$groupef.estimable]
      shown   <- not_est[seq_len(min(10L, length(not_est)))]
      suffix  <- if (length(not_est) > 10L)
        sprintf(", ... (%d more)", length(not_est) - 10L) else ""
      cat(sprintf("    full-rank but not estimable: %s%s\n",
                  paste(shown, collapse = ", "), suffix))
    }
  }
  cat("\n")

  # ---- Section 2: Level 2 (across-group) ------------------------------------
  cat("--- Level 2 (Across-Group) ---\n")

  w <- max(nchar(re_names))

  for (nm in re_names) {
    Xj    <- x$W[[nm]]
    other <- setdiff(colnames(Xj), "(Intercept)")

    hyper_rhs <- if (length(other) == 0L) "1" else paste(c("1", other), collapse = " + ")

    cat(sprintf("  %-*s ~ %s\n", w, nm, hyper_rhs))
  }
  cat("\n")

  # ---- Section 3: Hyper-design rank (estimable groups only) -----------------
  if (!is.null(x$popef.rank) && !is.null(x$groupef.estimable)) {
    n_est_groups <- sum(x$groupef.estimable)
    cat("--- Level 2 Rank ---\n")
    cat(sprintf("  (Restricted to %d estimable %s)\n\n", n_est_groups, grp))
    deficient_nms <- character(0)
    for (nm in re_names) {
      Xh      <- x$W[[nm]]
      p_hyper <- ncol(Xh)
      is_fr   <- if (nm %in% names(x$popef.rank)) x$popef.rank[[nm]] else NA
      status  <- if (isTRUE(is_fr)) "full-rank" else if (isFALSE(is_fr)) "RANK-DEFICIENT" else "unknown"
      cat(sprintf("  %-*s  groups=%-3d  predictors=%-2d  %s\n",
                  w, nm, n_est_groups, p_hyper, status))
      if (isFALSE(is_fr)) deficient_nms <- c(deficient_nms, nm)
    }
    # Per-RE deficient flags
    cat("\n")
    flag_strs <- ifelse(x$popef.deficient[re_names], "TRUE (deficient)", "FALSE")
    cat("  Rank-deficient flags:\n")
    for (nm in re_names) {
      cat(sprintf("    %-*s  %s\n", w, nm, flag_strs[nm]))
    }

    # Overall indicator
    ok_label <- if (isTRUE(x$popef.rank_ok)) "TRUE  -- model rank looks OK" else
                  "FALSE -- rank issues detected (see above)"
    cat(sprintf("\n  popef.rank_ok: %s\n", ok_label))

    if (length(deficient_nms) > 0L) {
      cat("\n")
      for (nm in deficient_nms) {
        cat(sprintf(
          "  NOTE: W for '%s' is rank-deficient after restricting to\n",
          nm))
        cat(sprintf(
          "  %d estimable %s. Consider removing predictors or merging\n",
          n_est_groups, grp))
        cat("  factor levels.\n")
      }
    }
    cat("\n")
  }

  invisible(x)
}
