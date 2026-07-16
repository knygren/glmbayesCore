#' Bayesian mixed model setup (single-factor \code{lmer}/\code{glmer} gate)
#'
#' Wrapper around \code{\link[lme4]{lmer}} or \code{\link[lme4]{glmer}} for
#' models with exactly one grouping factor. Design matrices come from
#' \code{formula} (including cross-level RE moderation terms). Random-effect
#' and residual variance components (\code{vcov_re}, \code{residual_var}) come
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
#' After fitting \code{lmer}, \code{model_setup} performs a two-step rank
#' check that assesses whether the model is empirically identified at both the
#' within-group and across-group levels:
#'
#' \enumerate{
#'   \item \emph{Level 1 (within-group):}  For each group \eqn{j}, the
#'     within-group random-effects design submatrix \eqn{\mathbf{Z}_j} is
#'     checked for full column rank (\code{re_rank}).  A rank-deficient group
#'     has too few distinct observations to estimate all random slopes
#'     independently; its BLUPs are identified through the prior rather than
#'     the data alone.  Such groups are flagged in \code{re_rank} but are
#'     retained in the \code{lmer} fit; \code{\link{Prior_Setup_lmebayes}}
#'     excludes them when calibrating priors.
#'
#'   \item \emph{Level 2 (across-group):}  Restricting to the full-rank groups
#'     from Step 1, each hyper-design matrix \code{X_hyper[[k]]} is checked
#'     for full column rank (\code{hyper_rank}).  Rank deficiency at this level
#'     means the level-2 hyperparameters \eqn{\boldsymbol{\mu}_k}---the prior
#'     means for random-effect coefficient \eqn{k} across groups---are not
#'     identified by the data, even as the number of full-rank groups grows.
#' }
#'
#' The scalar \code{rank_ok} is \code{TRUE} only when every
#' \code{X_hyper[[k]]} is full-rank after Step 2.  This is a necessary
#' condition for \code{\link{Prior_Setup_lmebayes}} to derive default priors
#' automatically; models with \code{rank_ok = FALSE} require user-supplied
#' hyperpriors.
#'
#' The example uses \code{big_word_club} from the Suggested package
#' \pkg{bayesrules} (see \code{?bayesrules::big_word_club}) and the same
#' formula as the full \code{lmerb()} demo in lmebayes
#' (\code{demo("Ex_12_lmerb_BigWordClub", package = "lmebayes")}).
#'
#' @param formula Mixed-model formula for design extraction and the reference
#'   \code{lmer}/\code{glmer} fit (fixed effects, hyper calibration, and
#'   variance components).
#' @param vcov_formula Ignored (deprecated). Variance components are taken from
#'   the full \code{formula} fit so \code{lmer} reference output is consistent.
#' @param data Optional data frame.
#' @param family A \code{\link[stats]{family}} object. Defaults to
#'   \code{gaussian()}, in which case \code{\link[lme4]{lmer}} is used.
#'   Non-Gaussian families use \code{\link[lme4]{glmer}}.
#' @param REML Logical; passed to \code{\link[lme4]{lmer}} when
#'   \code{family = gaussian()}.
#' @param control \code{\link[lme4]{lmerControl}} when \code{family = gaussian()},
#'   otherwise \code{\link[lme4]{glmerControl}}; passed through to the reference
#'   fit when \code{fit_mer = TRUE}.
#' @param fit_mer If \code{TRUE} (default), fit reference \code{lmer}/\code{glmer}
#'   models and extract variance components. If \code{FALSE}, return design
#'   matrices and rank diagnostics only (used by \code{glmerb()} in lmebayes).
#' @param start Optional starting values for the inner optimization.
#' @param verbose Passed to \code{\link[lme4]{lmer}}.
#' @param subset,weights,na.action,offset,contrasts Passed to
#'   \code{\link[lme4]{lmer}}.
#' @param devFunOnly If \code{TRUE}, return the deviance function only (Gaussian
#'   \code{lmer} fits only).
#' @param ... Passed to design extraction and, when \code{fit_mer = TRUE}, to the
#'   reference \code{lmer}/\code{glmer} fit.
#' @return Object of class \code{"model_setup"}: \code{y}, \code{Z},
#'   \code{groups}, \code{X_hyper}, \code{formula}, \code{family},
#'   \code{vcov_formula} (deprecated alias of \code{formula}),
#'   \code{lmer_fit} / \code{glmer_fit}, \code{lmer_vcov_fit} /
#'   \code{glmer_vcov_fit} (same object as the full-formula fit),
#'   \code{varcorr}, \code{vcov_re}, \code{residual_var}, and \code{re_rank} (named logical
#'   vector: \code{TRUE} if \code{Z_j} is full column rank for that group).
#' @seealso \code{\link{extract_re_hyper_matrices}},
#'   \code{\link{lmerb_default_vcov_formula}},
#'   \code{\link{extract_lmer_variance_components}}
#' @export
model_setup <- function(
    formula,
    data = NULL,
    vcov_formula = NULL,
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
    ...
) {
  cl <- match.call()
  family <- .lmebayes_normalize_family(family)
  is_gaussian <- identical(family$family, "gaussian")
  if (is.null(control) && is_gaussian) {
    control <- lme4::lmerControl()
  }

  design <- extract_re_hyper_matrices(formula = formula, data = data, ...)
  design$call    <- cl
  design$formula <- formula
  design$family  <- family

  if (!is.null(vcov_formula)) {
    warning(
      "'vcov_formula' is deprecated and ignored; variance components use ",
      "the full 'formula' reference fit.",
      call. = FALSE
    )
  }
  design$vcov_formula <- formula

  if (isTRUE(fit_mer)) {
    mer_args <- c(
      list(
        data = data,
        verbose = verbose
      ),
      if (!is.null(control)) list(control = control),
      .lmebayes_mer_optional_args(
        start = start,
        subset = subset,
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

    vc <- extract_mer_variance_components(
      fit_full,
      design$re_coef_names
    )
    if (is_gaussian) {
      design$lmer_fit <- fit_full
      design$lmer_vcov_fit <- fit_full
    } else {
      design$glmer_fit <- fit_full
      design$glmer_vcov_fit <- fit_full
    }
    design$varcorr <- vc$varcorr
    design$vcov_re <- vc$vcov_re
    design$residual_var <- vc$residual_var
  }

  # Per-group rank check: is Z_j full column rank for each factor level?
  p_re  <- ncol(design$Z)
  g_chr <- as.character(design$groups)
  design$re_rank <- vapply(
    levels(design$groups),
    function(lev) {
      rows <- which(g_chr == lev)
      Z_j  <- design$Z[rows, , drop = FALSE]
      nrow(Z_j) >= p_re &&
        Matrix::rankMatrix(Z_j, method = "qr")[1L] == p_re
    },
    logical(1L)
  )

  # Hyper-design rank check: for each RE coefficient, is the level-2 design
  # matrix X_hyper[[nm]] full column rank when restricted to the full-rank
  # groups?  Rank-deficient groups contribute a zero BLUP for the missing
  # slope and are excluded here so the check reflects only groups that
  # actually supply information about each RE.
  full_rank_levs <- names(design$re_rank)[design$re_rank]
  design$hyper_rank <- vapply(
    design$re_coef_names,
    function(nm) {
      Xh <- design$X_hyper[[nm]][full_rank_levs, , drop = FALSE]
      p  <- ncol(Xh)
      nrow(Xh) >= p && Matrix::rankMatrix(Xh, method = "qr")[1L] == p
    },
    logical(1L)
  )

  # Convenience summaries:
  #   hyper_deficient : named logical, TRUE = that RE's hyper-matrix is
  #                     rank-deficient (inverse of hyper_rank)
  #   rank_ok         : scalar TRUE only when every Z_j AND every hyper-matrix
  #                     is full-rank -- a quick go/no-go indicator
  design$hyper_deficient <- !design$hyper_rank

  # rank_ok reflects only the hyper-design matrices (level-2 estimability):
  # TRUE  = all X_hyper are full-rank after restricting to full-rank groups
  #         => the random-effects model can be estimated
  # FALSE = at least one X_hyper is rank-deficient => hyper parameters are
  #         not identified; Z_j rank deficiency is reported separately above
  design$rank_ok <- all(design$hyper_rank)

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
.lmebayes_normalize_family <- function(family) {
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  if (!inherits(family, "family") || is.null(family$family)) {
    stop("'family' must be a family object.", call. = FALSE)
  }
  family
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

#' @rdname model_setup
#' @method print model_setup
#' @param x A \code{model_setup} object.
#' @param ... Ignored.
#' @export
print.model_setup <- function(x, ...) {

  resp     <- deparse(x$formula[[2L]])
  re_names <- x$re_coef_names
  grp      <- x$group_name
  n_obs    <- length(x$y)
  n_lev    <- nlevels(x$groups)

  # ---- Call ------------------------------------------------------------------
  if (!is.null(x$call)) {
    cat("Call:\n  ", deparse1(x$call), "\n\n", sep = "")
  }

  # ---- Section 1: Measurement Model -----------------------------------------
  cat("--- Measurement Model ---\n")
  cat(sprintf("  %s ~ %s\n\n", resp, paste(re_names, collapse = " + ")))
  cat(sprintf("  Observations : %d\n", n_obs))
  cat(sprintf("  RE predictors: %d\n", length(re_names)))
  cat(sprintf("  Group        : %s  [%d levels]\n", grp, n_lev))
  if (!is.null(x$re_rank)) {
    n_full <- sum(x$re_rank)
    cat(sprintf("  Full-rank Z_j: %d of %d groups\n", n_full, n_lev))
    if (n_full < n_lev) {
      deficient <- names(x$re_rank)[!x$re_rank]
      shown     <- deficient[seq_len(min(10L, length(deficient)))]
      suffix    <- if (length(deficient) > 10L)
        sprintf(", ... (%d more)", length(deficient) - 10L) else ""
      cat(sprintf("    rank-deficient: %s%s\n",
                  paste(shown, collapse = ", "), suffix))
    }
  }
  cat("\n")

  # ---- Section 2: Random Effects Model --------------------------------------
  cat("--- Random Effects Model ---\n")

  w <- max(nchar(re_names))

  for (nm in re_names) {
    Xj    <- x$X_hyper[[nm]]
    other <- setdiff(colnames(Xj), "(Intercept)")

    hyper_rhs <- if (length(other) == 0L) "1" else paste(c("1", other), collapse = " + ")

    cat(sprintf("  %-*s ~ %s\n", w, nm, hyper_rhs))
  }
  cat("\n")

  # ---- Section 3: Hyper-design rank (full-rank groups only) -----------------
  if (!is.null(x$hyper_rank) && !is.null(x$re_rank)) {
    n_full_groups <- sum(x$re_rank)
    cat("--- Random Effects Model: Hyper-Design Rank ---\n")
    cat(sprintf("  (Restricted to %d full-rank %s)\n\n", n_full_groups, grp))
    deficient_nms <- character(0)
    for (nm in re_names) {
      Xh      <- x$X_hyper[[nm]]
      p_hyper <- ncol(Xh)
      is_fr   <- if (nm %in% names(x$hyper_rank)) x$hyper_rank[[nm]] else NA
      status  <- if (isTRUE(is_fr)) "full-rank" else if (isFALSE(is_fr)) "RANK-DEFICIENT" else "unknown"
      cat(sprintf("  %-*s  groups=%-3d  predictors=%-2d  %s\n",
                  w, nm, n_full_groups, p_hyper, status))
      if (isFALSE(is_fr)) deficient_nms <- c(deficient_nms, nm)
    }
    # Per-RE deficient flags
    cat("\n")
    flag_strs <- ifelse(x$hyper_deficient[re_names], "TRUE (deficient)", "FALSE")
    cat("  Rank-deficient flags:\n")
    for (nm in re_names) {
      cat(sprintf("    %-*s  %s\n", w, nm, flag_strs[nm]))
    }

    # Overall indicator
    ok_label <- if (isTRUE(x$rank_ok)) "TRUE  -- model rank looks OK" else
                  "FALSE -- rank issues detected (see above)"
    cat(sprintf("\n  rank_ok: %s\n", ok_label))

    if (length(deficient_nms) > 0L) {
      cat("\n")
      for (nm in deficient_nms) {
        cat(sprintf(
          "  NOTE: X_hyper for '%s' is rank-deficient after restricting to\n",
          nm))
        cat(sprintf(
          "  %d full-rank %s. Consider removing predictors or merging\n",
          n_full_groups, grp))
        cat("  factor levels.\n")
      }
    }
    cat("\n")
  }

  invisible(x)
}
