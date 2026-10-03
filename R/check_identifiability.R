## Family-specific "quick reject" pre-check for classical GLM MLE existence
## (or, for Gamma, dispersion estimability), evaluated *before* attempting
## the per-group glm() fit.  Returns a single note string, or character(0)
## if the group passes.  @noRd
.lmebayes_glm_estimable_precheck <- function(fam_name, y_j, n_j, p_j) {
  switch(
    fam_name,
    binomial = if (length(unique(y_j)) < 2L) "single outcome level",
    poisson  = if (all(y_j == 0)) {
      "all-zero response (MLE diverges under log link)"
    },
    Gamma = {
      if (any(y_j <= 0)) {
        "non-positive response (invalid for Gamma family)"
      } else if (n_j <= p_j) {
        "zero residual degrees of freedom for dispersion (n_j <= p_j)"
      }
    },
    NULL
  )
}

## Families for which a per-group glm() fit is attempted (mean-model MLE
## existence, plus family-specific pre-checks above).
.lmebayes_glm_fit_families <- c("binomial", "poisson", "Gamma")

#' Per-group classical-GLM MLE / dispersion-estimability check for Block-1
#'
#' For \code{binomial}, \code{poisson}, and \code{Gamma} families, fits
#' \code{glm(y ~ D_j - 1, family = family)} on each group with algebraically
#' full-rank \code{D_j} and marks the group estimable when a family-specific
#' pre-check passes (see below) \emph{and} all coefficients and \code{vcov}
#' entries from the fit are finite.  For \code{gaussian}, no fit is
#' attempted -- the coefficient MLE (OLS) always exists given full column
#' rank, but the group-level residual \emph{dispersion} additionally
#' requires strictly more observations than parameters (\code{n_j > p_j});
#' this is checked directly via arithmetic on \code{n_j}/\code{p_j}. Other
#' families return \code{groupef.estimable = groupef.rank} (no check).
#'
#' Family-specific pre-checks (before any fit is attempted):
#' \itemize{
#'   \item \code{binomial}: single outcome level, or a convergent
#'     \code{glm} fit whose fitted probabilities hit the 0/1 boundary
#'     (practical signal of (quasi-)complete separation).
#'   \item \code{poisson}: all-zero response (log-link intercept MLE
#'     diverges to \eqn{-\infty}).
#'   \item \code{Gamma}: non-positive response (invalid domain), or
#'     \code{n_j <= p_j} (no residual degrees of freedom for the
#'     dispersion/shape parameter).
#'   \item \code{gaussian}: \code{n_j <= p_j} (no residual degrees of
#'     freedom for the variance).
#' }
#' @noRd
.lmebayes_block_glm_estimable <- function(y, groups, Z, groupef.rank, family) {
  g_levs <- names(groupef.rank)
  if (is.null(g_levs)) {
    g_levs <- levels(groups)
    names(groupef.rank) <- g_levs
  }

  groupef.estimable <- stats::setNames(rep(FALSE, length(g_levs)), g_levs)
  fam_name <- family$family

  is_fit_family <- fam_name %in% .lmebayes_glm_fit_families
  is_gaussian   <- identical(fam_name, "gaussian")

  if (!is_fit_family && !is_gaussian) {
    groupef.estimable[groupef.rank] <- TRUE
    return(list(
      groupef.estimable = groupef.estimable,
      groupef.glm_check = NULL
    ))
  }

  y <- as.numeric(y)
  g_chr <- as.character(groups)
  fr_levs <- g_levs[groupef.rank]

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
    } else if (is_gaussian) {
      ## OLS coefficient MLE always exists given full column rank; only the
      ## residual degrees of freedom for the dispersion is at issue here --
      ## no fit is attempted.
      if (n_j <= p_j) {
        note <- "zero residual degrees of freedom for dispersion (n_j <= p_j)"
      } else {
        estimable <- TRUE
      }
    } else {
      note <- .lmebayes_glm_estimable_precheck(fam_name, y_j, n_j, p_j)
      if (length(note) == 0L) {
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
            } else if (identical(fam_name, "binomial") &&
                       any(fit$fitted.values < 1e-7 |
                             fit$fitted.values > 1 - 1e-7)) {
              ## glm() can still return finite coef/vcov under
              ## (quasi-)complete separation; boundary fitted
              ## probabilities are the practical MLE-existence signal.
              note <- "possible (quasi-)complete separation (fitted probabilities at boundary)"
            } else {
              estimable <- TRUE
            }
          }
        }
      }
    }

    groupef.estimable[[lev]] <- estimable
    rows_out[[ii]] <- data.frame(
      group     = lev,
      n         = n_j,
      p         = p_j,
      groupef.rank   = TRUE,
      estimable = estimable,
      note      = if (length(note)) paste(note, collapse = "; ") else "",
      stringsAsFactors = FALSE
    )
  }

  groupef.glm_check <- if (length(rows_out)) {
    do.call(rbind, rows_out)
  } else {
    data.frame(
      group = character(0),
      n = integer(0),
      p = integer(0),
      groupef.rank = logical(0),
      estimable = logical(0),
      note = character(0),
      stringsAsFactors = FALSE
    )
  }

  list(groupef.estimable = groupef.estimable, groupef.glm_check = groupef.glm_check)
}

#' Check identifiability and estimability of a single-factor mixed-model design
#'
#' @description
#' Runs the "two-step identifiability assessment" used internally by
#' \code{\link{model_setup}}: a Level~1 (within-group) rank/estimability
#' check on the level-1 design \code{D}, followed by a Level~2 (across-group)
#' rank check on the level-2 design matrices \code{W}, restricted to the
#' Level-1-estimable groups. Takes the same \code{D}/\code{group}/\code{W}
#' matrix inputs (and the same conventions -- see Details) as
#' \samp{rLMM_reg()}/\samp{rGLMM_reg()}, so it can be run directly
#' on hand-built matrices from those matrix-level workflows, or on
#' \code{model_setup()}'s \code{design$D}/\code{design$group}/
#' \code{design$W} directly, without requiring a full \code{lmer}/\code{glmer}
#' reference fit.
#'
#' @details
#' \strong{Level 1 (within-group).} For each group \eqn{j}, the within-group
#' submatrix \eqn{D_j} (\code{D[group == } \eqn{j}\code{, ]}) is checked for
#' full column rank (\code{groupef.rank}). A rank-deficient group has too few
#' distinct observations to estimate all random slopes independently.
#' Algebraically full-rank groups are then additionally checked for
#' estimability (\code{groupef.estimable}, with per-group detail in
#' \code{groupef.glm_check}), with the specific check depending on
#' \code{family}:
#' \itemize{
#'   \item \code{binomial()}, \code{poisson()}, \code{Gamma()}: a classical
#'     \code{glm(y ~ D_j - 1, family)} fit is attempted; a group can be full
#'     column rank yet have no finite MLE under complete or quasi-complete
#'     separation (\code{binomial}), an all-zero count response
#'     (\code{poisson}), or a non-positive response (\code{Gamma}) --
#'     \code{Gamma} groups are additionally required to have \code{n_j > p_j}
#'     so the dispersion/shape parameter is estimable.
#'   \item \code{gaussian()}: no fit is attempted (the OLS coefficient MLE
#'     always exists given full column rank), but the group must still have
#'     \code{n_j > p_j} so the residual variance/dispersion is estimable
#'     (otherwise the group achieves a perfect/saturated fit with zero
#'     residual degrees of freedom).
#'   \item other families: \code{groupef.estimable} is set equal to
#'     \code{groupef.rank} (no additional check yet).
#' }
#'
#' \strong{Level 2 (across-group).} Restricting to the \strong{estimable}
#' groups from Level 1 (\code{groupef.estimable}), each level-2 design matrix
#' \code{W[[k]]} is checked for full column rank (\code{popef.rank}), using
#' only the rows at the estimable positions (row \eqn{j} of \code{W[[k]]}
#' corresponds to \code{levels(group)[j]}, exactly as in
#' \samp{rLMM_reg()}). Rank deficiency at this level means the level-2
#' hyperparameters are not identified by the data, even as the number of
#' estimable groups grows. The scalar \code{popef.rank_ok} is \code{TRUE}
#' only when every \code{popef.rank} entry is \code{TRUE}.
#'
#' \strong{Inputs must already be consistently keyed.} This function does
#' \emph{not} infer or coerce missing names -- it validates them up front and
#' errors immediately on a mismatch, rather than risking a silently wrong
#' match deeper in the rank/estimability computation:
#' \itemize{
#'   \item \code{D} must have unique, non-empty \code{colnames(D)}: these are
#'     the random-effect coefficient names used to key \code{W} (there is no
#'     separate \code{groupef.names} argument to override them, matching
#'     \samp{rLMM_reg()}).
#'   \item \code{group} must be a \code{factor} (not coerced from a character
#'     vector); \code{levels(group)} fixes the row order assumed by
#'     \code{W} (there is no separate \code{group_levels} argument).
#'   \item \code{W} must be a named list, with names matching
#'     \code{colnames(D)} exactly (as sets; no missing, empty, or duplicate
#'     names), and each \code{W[[k]]} must have exactly
#'     \code{length(levels(group))} rows.
#'   \item \code{group_name}, used only for error/verbose messages, must be
#'     resolvable from either the \code{group_name} argument or
#'     \code{attr(group, "group_name")}.
#' }
#'
#' @param y Response vector of length \code{nrow(D)}.
#' @param D Level-1 design matrix (\code{nrow(D) x p_re}); must have unique,
#'   non-empty \code{colnames(D)} (random-effect coefficient names).
#' @param W Named list of group-level design matrices (\code{J x q_k}), one per
#'   column of \code{D}; names must match \code{colnames(D)} exactly.
#' @param group Grouping factor of length \code{nrow(D)} (must be a
#'   \code{factor}); \code{levels(group)} fixes the row order assumed by
#'   \code{W} and names \code{groupef.rank}/\code{groupef.estimable} in the
#'   return value -- there is no separate \code{group_levels} argument.
#'   Unlike \samp{rLMM_reg()}, \code{group_name} is a separate argument
#'   here (see below) rather than derived from \code{group} via
#'   \code{substitute()}.
#' @param family A \code{\link[stats]{family}} object (or a string/function
#'   accepted by \code{\link[stats]{family}}); normalized internally via the
#'   same logic \code{\link{model_setup}} uses. Defaults to \code{gaussian()}.
#' @param group_name Display name of the grouping factor (e.g.
#'   \code{"school_id"}), used only in error and \code{verbose} messages
#'   (not in any computation). Defaults to \code{attr(group, "group_name")}
#'   (the convention already used by \samp{rLMM_reg()}); must be
#'   resolvable to a non-empty string, either this way or by passing it
#'   explicitly. Unlike \samp{rLMM_reg()}, there is no fallback to
#'   \code{substitute(group)} -- this function may be called on matrices that
#'   are not bare variables.
#' @param verbose If \code{TRUE}, emit \code{message()}s summarizing the
#'   Level 1 and Level 2 outcomes. Default \code{FALSE} (silent), so calling
#'   this from \code{\link{model_setup}} introduces no new console output.
#' @return A list with:
#'   \describe{
#'     \item{\code{groupef.rank}}{Named logical, one entry per group (names
#'       \code{levels(group)}): \code{TRUE} if \eqn{D_j} is full column
#'       rank.}
#'     \item{\code{groupef.estimable}}{Named logical, same names as
#'       \code{groupef.rank}: \code{TRUE} if the group is additionally
#'       estimable (see Details).}
#'     \item{\code{groupef.glm_check}}{Per-group diagnostic data frame for
#'       \code{binomial()}/\code{poisson()}/\code{Gamma()}/\code{gaussian()},
#'       or \code{NULL} for other families.}
#'     \item{\code{popef.rank}}{Named logical, one entry per RE coefficient
#'       (names \code{colnames(D)}): \code{TRUE} if \code{W[[k]]} is full
#'       column rank when restricted to \code{groupef.estimable} groups.}
#'     \item{\code{popef.deficient}}{Negation of \code{popef.rank}.}
#'     \item{\code{popef.rank_ok}}{Scalar \code{TRUE} only when every
#'       \code{popef.rank} entry is \code{TRUE}.}
#'     \item{\code{group_levels}}{\code{levels(group)}, i.e. the exact
#'       names used on \code{groupef.rank}/\code{groupef.estimable}.}
#'     \item{\code{group_name}}{The resolved grouping-factor display name.}
#'   }
#' @seealso \code{\link{model_setup}},
#'   \samp{rLMM_reg()},
#'   \samp{rGLMM_reg()},
#'   \code{\link{extract_re_hyper_matrices}}
#' @examples
#' set.seed(1)
#' J <- 8L; n_j <- 6L
#' group <- factor(rep(paste0("g", seq_len(J)), each = n_j))
#' x1 <- rnorm(J * n_j)
#' D <- cbind("(Intercept)" = 1, x1 = x1)
#' y <- 1 + x1 + rnorm(J * n_j)
#' W <- list(
#'   "(Intercept)" = matrix(1, J, 1, dimnames = list(NULL, "(Intercept)")),
#'   x1            = matrix(1, J, 1, dimnames = list(NULL, "(Intercept)"))
#' )
#' ident <- check_identifiability(
#'   y = y, D = D, group = group, W = W,
#'   family = gaussian(), group_name = "group"
#' )
#' ident$popef.rank_ok
#' @export
check_identifiability <- function(
    y,
    D,
    group,
    W,
    family = gaussian(),
    group_name = attr(group, "group_name"),
    verbose = FALSE
) {
  family <- .lmebayes_normalize_family(family)

  if (!is.factor(group)) {
    stop(
      "'group' must be a factor (wrap with factor(group, levels = ...) to ",
      "control level order or supply a fixed superset of levels); there is ",
      "no 'group_levels' argument to override this.",
      call. = FALSE
    )
  }
  if (!is.character(group_name) || length(group_name) != 1L ||
      is.na(group_name) || !nzchar(group_name)) {
    stop(
      "'group_name' could not be resolved to a single non-empty string; ",
      "pass it explicitly or set attr(group, 'group_name').",
      call. = FALSE
    )
  }

  D <- as.matrix(D)
  n_obs <- nrow(D)
  if (length(y) != n_obs || length(group) != n_obs) {
    stop(
      sprintf(
        "'y' (length %d), 'D' (%d rows), and 'group' (length %d) must all agree.",
        length(y), n_obs, length(group)
      ),
      call. = FALSE
    )
  }

  re_names <- colnames(D)
  if (is.null(re_names) || length(re_names) != ncol(D) || anyNA(re_names) ||
      any(!nzchar(re_names)) || anyDuplicated(re_names)) {
    stop(
      "'D' must have unique, non-empty column names (colnames(D)); there is ",
      "no 'groupef.names' argument to override this.",
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

  g_levs <- levels(group)
  for (nm in re_names) {
    Wk <- W[[nm]]
    if (!is.matrix(Wk) && !inherits(Wk, "Matrix")) {
      stop(sprintf("W[['%s']] must be a matrix.", nm), call. = FALSE)
    }
    if (nrow(Wk) != length(g_levs)) {
      stop(
        sprintf(
          "W[['%s']] must have exactly one row per level of %s (%d), got %d.",
          nm, group_name, length(g_levs), nrow(Wk)
        ),
        call. = FALSE
      )
    }
  }

  ## Level 1: per-group rank check -- is D_j full column rank for each level?
  p_re  <- ncol(D)
  g_chr <- as.character(group)
  groupef.rank <- vapply(
    g_levs,
    function(lev) {
      rows <- which(g_chr == lev)
      D_j  <- D[rows, , drop = FALSE]
      nrow(D_j) >= p_re &&
        Matrix::rankMatrix(D_j, method = "qr")[1L] == p_re
    },
    logical(1L)
  )

  if (isTRUE(verbose)) {
    message(sprintf(
      "check_identifiability: Level 1 (rank) -- %d of %d %s full column rank.",
      sum(groupef.rank), length(g_levs), group_name
    ))
  }

  ## Level 1: estimability (classical-glm MLE / residual-df, per family) --
  ## must run before the Level 2 check below, since Level 2 restricts to
  ## *estimable* groups, not merely algebraically full-rank ones.
  glm_est <- .lmebayes_block_glm_estimable(
    y = y, groups = group, Z = D, groupef.rank = groupef.rank, family = family
  )
  groupef.estimable <- glm_est$groupef.estimable
  groupef.glm_check <- glm_est$groupef.glm_check

  if (!is.null(groupef.glm_check) && nrow(groupef.glm_check) > 0L &&
      !all(groupef.glm_check$group %in% g_levs)) {
    stop(
      "Internal error: groupef.glm_check$group contains values not in ",
      "levels(", group_name, ").",
      call. = FALSE
    )
  }

  if (isTRUE(verbose)) {
    message(sprintf(
      "check_identifiability: Level 1 (estimability) -- %d of %d full-rank %s are estimable.",
      sum(groupef.estimable[groupef.rank]), sum(groupef.rank), group_name
    ))
  }

  ## Level 2: hyper-design rank check, restricted to Level-1-estimable groups.
  ## Positional, not name-based: row j of W[[k]] corresponds to
  ## levels(group)[j] (the rLMM_reg/rGLMM_reg convention), and
  ## groupef.estimable preserves that same level order.
  estimable_idx <- which(groupef.estimable)
  popef.rank <- vapply(
    re_names,
    function(nm) {
      Wk <- W[[nm]][estimable_idx, , drop = FALSE]
      p  <- ncol(Wk)
      nrow(Wk) >= p && Matrix::rankMatrix(Wk, method = "qr")[1L] == p
    },
    logical(1L)
  )
  popef.deficient <- !popef.rank
  popef.rank_ok <- all(popef.rank)

  if (isTRUE(verbose)) {
    if (isTRUE(popef.rank_ok)) {
      message(sprintf(
        "check_identifiability: Level 2 -- all hyper-design matrices are full rank (%d estimable %s).",
        length(estimable_idx), group_name
      ))
    } else {
      message(sprintf(
        "check_identifiability: Level 2 -- rank-deficient hyper-design matrix(es): %s.",
        paste(names(popef.rank)[popef.deficient], collapse = ", ")
      ))
    }
  }

  list(
    groupef.rank      = groupef.rank,
    groupef.estimable = groupef.estimable,
    groupef.glm_check = groupef.glm_check,
    popef.rank        = popef.rank,
    popef.deficient   = popef.deficient,
    popef.rank_ok     = popef.rank_ok,
    group_levels      = g_levs,
    group_name        = group_name
  )
}
