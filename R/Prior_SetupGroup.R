#' Prior setup for row-block regressions
#'
#' Runs \code{\link{Prior_Setup}} on each block subset of the data.
#' Typical callers are \code{\link{rNormal_reg_group}} /
#' \code{\link{rNormalGLM_reg_group}} or the \code{lmbBlock} / \code{glmbBlock}
#' wrappers in \pkg{lmebayes}.
#'
#' @param formula A \code{\link{formula}} with a single response.
#' @param group Group partition: \code{factor} or vector of length \code{nrow(data)}
#'   (after \code{model.frame}), a column name in \code{data}, \code{l2_blocks}
#'   counts, or a list of row index vectors (see \code{\link{normalize_group}}).
#' @param pwt,n_prior,sd,mu,dispersion,k Each of these six calibration
#'   arguments (passed to \code{\link{Prior_Setup}} for every block) may be
#'   supplied in \strong{either} of two forms:
#'   \itemize{
#'     \item the ordinary shared form documented in \code{\link{Prior_Setup}}
#'       (a scalar, or a per-coefficient vector for \code{pwt}/\code{sd}/\code{mu})
#'       -- applied identically to every block (default); or
#'     \item a \strong{named list with one element per block}, keyed by the
#'       block/group level IDs (\code{names(Prior_SetupGroup(...))}, i.e.
#'       \code{normalize_group()$ids}) -- each element supplies that block's
#'       own value, in the same shape \code{Prior_Setup()} would otherwise
#'       accept. \code{names(x)} must exactly match the block IDs. An element
#'       may be \code{NULL}, meaning use \code{Prior_Setup()}'s default for
#'       that argument in that block.
#'   }
#'   The two forms may be mixed freely across the six arguments in one call.
#' @inheritParams Prior_Setup
#' @details
#' ## What it produces
#'
#' \code{Prior_SetupGroup()} calls \code{\link{Prior_Setup}} once per block,
#' on that block's rows only, and returns the results as a named list keyed by
#' block ID.
#'
#' ## Contrast with hierarchical mixed-model setup
#'
#' \code{\link[lmebayesCore:Prior_Setup_GLMM]{Prior_Setup_GLMM}} calibrates a
#' \strong{hierarchical} prior (one reference \code{lmer}/\code{glmer} fit,
#' Block~2 hyperpriors, cross-group shrinkage). \code{Prior_SetupGroup()}
#' calibrates \strong{independent} priors per row block with no shared
#' \eqn{\gamma} or \eqn{\tau^2}; row-block engines draw blocks iid.
#'
#' Names for per-block calibration lists must match \code{\link{normalize_group}}
#' block IDs exactly.
#'
#' ## Prior families
#'
#' After calibration, convert to sampler-ready \code{\link{pfamily}} objects with
#' \code{pfamily_list(x, ptypes = ...)} (see **Prior families via
#' \code{pfamily_list()}** below for \code{ptypes} and allowed families).
#'
#' @return A named list of class \code{"Prior_SetupGroup"}. Each element is a
#'   \code{\link{Prior_Setup}} result for one block.
#' @seealso \code{\link{Prior_Setup}}, \code{\link{multi_prior_setup}},
#'   \code{\link{normalize_group}}, \code{\link{pfamily_list}},
#'   \code{\link{rNormal_reg_group}}, \code{\link{rNormalGLM_reg_group}},
#'   \code{\link[lmebayesCore:Prior_Setup_GLMM]{Prior_Setup_GLMM}}
#' @example inst/examples/Ex_Prior_SetupGroup.R
#' @export
Prior_SetupGroup <- function(
    formula,
    group,
    family = gaussian(),
    data = NULL,
    weights = NULL,
    subset = NULL,
    na.action = na.fail,
    offset = NULL,
    contrasts = NULL,
    pwt = NULL,
    pwt_default_low = 0.01,
    pwt_default_high = 0.05,
    n_prior = NULL,
    sd = NULL,
    dispersion = NULL,
    intercept_source = c("null_model", "full_model"),
    effects_source = c("null_effects", "full_model"),
    mu = NULL,
    k = 1,
    ...
) {
  call <- match.call()
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  fam_ok <- family$family %in% c("gaussian", "poisson", "binomial")
  if (is.null(family$family) || !fam_ok) {
    stop(
      "Prior_SetupGroup() supports family = gaussian(), poisson(), or binomial() only.",
      call. = FALSE
    )
  }
  if (missing(data)) {
    data <- environment(formula)
  }

  meta <- .prior_setup_group_formula_block_meta(
    formula = formula,
    group = group,
    data = data,
    subset = if (!missing(subset)) subset else NULL,
    weights = if (!missing(weights)) weights else NULL,
    na.action = if (!missing(na.action)) na.action else NULL,
    offset = if (!missing(offset)) offset else NULL,
    contrasts = if (!missing(contrasts)) contrasts else NULL
  )
  block_ids <- meta$group_info$ids

  pwt_r        <- .prior_setup_group_resolve_calibration_arg(pwt, block_ids, "pwt")
  n_prior_r    <- .prior_setup_group_resolve_calibration_arg(n_prior, block_ids, "n_prior")
  sd_r         <- .prior_setup_group_resolve_calibration_arg(sd, block_ids, "sd")
  mu_r         <- .prior_setup_group_resolve_calibration_arg(mu, block_ids, "mu")
  dispersion_r <- .prior_setup_group_resolve_calibration_arg(dispersion, block_ids, "dispersion")
  k_r          <- .prior_setup_group_resolve_calibration_arg(k, block_ids, "k")

  ps_args <- list(
    family = family,
    data = data,
    weights = weights,
    na.action = na.action,
    offset = offset,
    contrasts = contrasts,
    pwt_default_low = pwt_default_low,
    pwt_default_high = pwt_default_high,
    intercept_source = intercept_source,
    effects_source = effects_source
  )

  setups <- vector("list", meta$group_info$k)
  for (b in seq_len(meta$group_info$k)) {
    block_id <- block_ids[[b]]
    rows_b <- .prior_setup_group_rows_to_data_subset(
      meta$group_info$rows[[b]], meta$mf, data
    )
    setups[[b]] <- tryCatch(
      do.call(
        Prior_Setup,
        c(
          list(formula = formula, subset = rows_b),
          ps_args,
          list(
            pwt        = pwt_r[[block_id]],
            n_prior    = n_prior_r[[block_id]],
            sd         = sd_r[[block_id]],
            dispersion = dispersion_r[[block_id]],
            mu         = mu_r[[block_id]],
            k          = k_r[[block_id]]
          ),
          list(...)
        )
      ),
      error = function(e) {
        stop(
          "Prior_SetupGroup(): block '", block_id, "': ", conditionMessage(e),
          call. = FALSE
        )
      }
    )
  }
  names(setups) <- block_ids

  attr(setups, "call") <- call
  attr(setups, "formula") <- formula
  attr(setups, "group") <- group
  attr(setups, "group_info") <- meta$group_info
  class(setups) <- c("Prior_SetupGroup", "list")
  setups
}

#' @rdname Prior_SetupGroup
#' @method print Prior_SetupGroup
#' @param x Object of class \code{"Prior_SetupGroup"}.
#' @param blocks \code{NULL} (header only), block IDs (\code{names(x)}), or
#'   integer indices. Selected blocks use \code{\link{print.PriorSetup}}.
#' @param ... Passed to \code{print.PriorSetup} when \code{blocks} is set.
#' @return \code{x} invisibly.
#' @export
print.Prior_SetupGroup <- function(x, blocks = NULL, ...) {
  cl <- attr(x, "call")
  cat("\nCall:\n")
  if (!is.null(cl)) {
    print(cl)
  } else {
    cat("Prior_SetupGroup()\n")
  }

  info <- attr(x, "group_info")
  k <- length(x)
  ids <- names(x)
  if (is.null(ids)) {
    ids <- if (!is.null(info$ids)) info$ids else as.character(seq_len(k))
  }

  cat("\n--- Row-block prior setup ---\n")
  form <- attr(x, "formula")
  if (!is.null(form)) {
    cat("  formula : ", paste(deparse(form), collapse = "\n"), "\n", sep = "")
  }
  cat(sprintf("  blocks  : %d (%s)\n", k, paste(ids, collapse = ", ")))
  cat("  Each element is a Prior_Setup() result; build pfamilies with\n")
  cat("  pfamily_list(). Use print(x, blocks = ...) for details.\n")

  if (!is.null(blocks)) {
    sel <- .prior_setup_group_select_named_list_keys(
      x, blocks, arg = "blocks", what = "group"
    )
    for (nm in sel) {
      cat(sprintf("\n[[%s]]\n", nm))
      print(x[[nm]], ...)
    }
  }

  invisible(x)
}

#' @rdname Prior_SetupGroup
#' @section Prior families via \code{pfamily_list()}:
#' Each row block (\code{names(object)}) maps to one \code{\link{pfamily}}:
#' \itemize{
#'   \item \code{"dNormal"} (default): \code{dNormal(mu, Sigma, dispersion)}
#'     from that block's \code{\link{Prior_Setup}} result.
#'   \item \code{"dNormal_Gamma"}: \code{dNormal_Gamma(mu, Sigma_0, shape, rate)}
#'     (Gaussian only).
#'   \item \code{"dIndependent_Normal_Gamma"}:
#'     \code{dIndependent_Normal_Gamma(mu, Sigma, shape_ING, rate)}
#'     (Gaussian only).
#' }
#' @param object A \code{Prior_SetupGroup} object.
#' @param ptypes Character: one string recycled to every block, or a named
#'   vector / list with one string per block ID (\code{names(object)}).
#'   \code{NULL} (default) uses \code{"dNormal"} for all blocks.
#' @param ... Not used.
#' @return An object of class \code{"pfamily_list"}: a named list of
#'   \code{pfamily} objects with attribute \code{"ptypes"}.
#' @export
#' @method pfamily_list Prior_SetupGroup
pfamily_list.Prior_SetupGroup <- function(object, ptypes = NULL, ...) {
  allowed <- c("dNormal", "dNormal_Gamma", "dIndependent_Normal_Gamma")
  block_ids <- names(object)
  if (is.null(block_ids) || any(!nzchar(block_ids))) {
    stop("'object' must be a named Prior_SetupGroup list.", call. = FALSE)
  }
  k <- length(block_ids)

  if (is.null(ptypes)) {
    ptypes <- "dNormal"
  }

  if (is.list(ptypes)) {
    ok <- vapply(
      ptypes,
      function(p) is.character(p) && length(p) == 1L && !is.na(p),
      logical(1L)
    )
    if (!all(ok)) {
      stop("'ptypes' list elements must each be a single string.",
           call. = FALSE)
    }
    nms <- names(ptypes)
    ptypes <- vapply(ptypes, identity, character(1L))
    names(ptypes) <- nms
  }
  if (!is.character(ptypes) || length(ptypes) < 1L || anyNA(ptypes)) {
    stop("'ptypes' must be a character vector or list of strings.",
         call. = FALSE)
  }
  bad <- setdiff(unique(ptypes), allowed)
  if (length(bad) > 0L) {
    stop(
      "Invalid 'ptypes' value(s): ", paste(bad, collapse = ", "),
      ". Allowed: ", paste(allowed, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (length(ptypes) == 1L) {
    ptypes <- stats::setNames(rep(unname(ptypes), k), block_ids)
  } else {
    if (length(ptypes) != k) {
      stop(
        sprintf(
          "'ptypes' has length %d but Prior_SetupGroup has %d block(s): %s.",
          length(ptypes), k, paste(block_ids, collapse = ", ")
        ),
        call. = FALSE
      )
    }
    if (!is.null(names(ptypes)) && any(nzchar(names(ptypes)))) {
      if (!setequal(names(ptypes), block_ids)) {
        stop(
          "Names of 'ptypes' must match the block IDs: ",
          paste(block_ids, collapse = ", "), ".",
          call. = FALSE
        )
      }
      ptypes <- ptypes[block_ids]
    } else {
      names(ptypes) <- block_ids
    }
  }

  out <- stats::setNames(vector("list", k), block_ids)
  for (id in block_ids) {
    ps <- object[[id]]
    out[[id]] <- switch(
      ptypes[[id]],
      dNormal = {
        disp <- ps$dispersion
        if (is.null(disp)) {
          disp <- 1
        }
        dNormal(
          mu         = ps$mu,
          Sigma      = ps$Sigma,
          dispersion = disp
        )
      },
      dNormal_Gamma = {
        if (is.null(ps$Sigma_0) || is.null(ps$shape) || is.null(ps$rate)) {
          stop(
            "Block '", id, "': ptypes = \"dNormal_Gamma\" requires ",
            "Sigma_0, shape, and rate on the Prior_Setup result ",
            "(family = gaussian()).",
            call. = FALSE
          )
        }
        dNormal_Gamma(
          mu      = ps$mu,
          Sigma_0 = ps$Sigma_0,
          shape   = ps$shape,
          rate    = ps$rate
        )
      },
      dIndependent_Normal_Gamma = {
        if (is.null(ps$Sigma) || is.null(ps$shape_ING) || is.null(ps$rate)) {
          stop(
            "Block '", id, "': ptypes = \"dIndependent_Normal_Gamma\" requires ",
            "Sigma, shape_ING, and rate on the Prior_Setup result ",
            "(family = gaussian()).",
            call. = FALSE
          )
        }
        dIndependent_Normal_Gamma(
          mu    = ps$mu,
          Sigma = ps$Sigma,
          shape = ps$shape_ING,
          rate  = ps$rate
        )
      }
    )
  }

  attr(out, "ptypes") <- ptypes
  class(out) <- c("pfamily_list", "list")
  out
}

#' @noRd
.prior_setup_group_resolve_calibration_arg <- function(x, block_ids, arg_name) {
  if (!is.list(x)) {
    return(stats::setNames(rep(list(x), length(block_ids)), block_ids))
  }
  nm <- names(x)
  if (is.null(nm) || any(!nzchar(nm))) {
    stop(
      "Prior_SetupGroup(): '", arg_name, "' is a list but is not fully named; ",
      "supply one named element per block/group level (",
      paste(block_ids, collapse = ", "), ").",
      call. = FALSE
    )
  }
  missing_ids <- setdiff(block_ids, nm)
  extra_ids   <- setdiff(nm, block_ids)
  if (length(missing_ids) || length(extra_ids)) {
    stop(
      "Prior_SetupGroup(): names('", arg_name, "') must exactly match the ",
      "block/group level IDs.",
      if (length(missing_ids)) {
        paste0("\n  Missing: ", paste(missing_ids, collapse = ", "))
      } else {
        ""
      },
      if (length(extra_ids)) {
        paste0("\n  Unexpected: ", paste(extra_ids, collapse = ", "))
      } else {
        ""
      },
      call. = FALSE
    )
  }
  x[block_ids]
}

#' @noRd
.prior_setup_group_formula_block_meta <- function(
    formula,
    group,
    data,
    subset = NULL,
    weights = NULL,
    na.action = NULL,
    offset = NULL,
    contrasts = NULL
) {
  mf_args <- list(
    formula = formula,
    data = data,
    drop.unused.levels = TRUE
  )
  if (!is.null(subset)) mf_args$subset <- subset
  if (!is.null(weights)) mf_args$weights <- weights
  if (!is.null(na.action)) mf_args$na.action <- na.action
  if (!is.null(offset)) mf_args$offset <- offset
  if (!is.null(contrasts)) mf_args$contrasts <- contrasts
  mf <- do.call(stats::model.frame, mf_args)

  l2 <- nrow(mf)
  block_vec <- .prior_setup_group_resolve_block(group, data, mf, l2)
  group_info <- normalize_group(block_vec, l2)

  mt <- attr(mf, "terms")
  x_mat <- stats::model.matrix(mt, mf, contrasts)
  p <- ncol(x_mat)
  pred_names <- colnames(x_mat)
  if (is.null(pred_names) || length(pred_names) != p) {
    pred_names <- paste0("X", seq_len(p))
  }

  list(
    mf = mf,
    group_info = group_info,
    p = p,
    pred_names = pred_names
  )
}

#' @noRd
.prior_setup_group_rows_to_data_subset <- function(rows_mf, mf, data) {
  rn <- rownames(mf)
  if (is.null(rn)) {
    return(rows_mf)
  }
  if (!is.null(rownames(data))) {
    out <- match(rn[rows_mf], rownames(data))
    if (anyNA(out)) {
      stop("Could not map model.frame rows to rownames(data).", call. = FALSE)
    }
    return(out)
  }
  as.integer(rn[rows_mf])
}

#' @noRd
.prior_setup_group_resolve_block <- function(block, data, mf, l2) {
  if (is.list(block)) {
    return(block)
  }
  if (is.character(block) && length(block) == 1L && block %in% names(data)) {
    rn_mf <- rownames(mf)
    col <- data[[block]]
    if (length(col) == l2 && is.null(rn_mf)) {
      return(col)
    }
    if (!is.null(rn_mf)) {
      if (!is.null(rownames(data))) {
        return(col[match(rn_mf, rownames(data))])
      }
      return(col[as.integer(rn_mf)])
    }
  }
  block <- as.vector(block)
  if (length(block) == l2) {
    return(block)
  }
  stop(
    "'block' must have length nrow(model.frame), be a list of row indices, ",
    "l2_blocks counts, or a single column name in 'data'.",
    call. = FALSE
  )
}

#' @noRd
.prior_setup_group_select_named_list_keys <- function(x, keys, arg = "groups",
                                                      what = "name") {
  all_names <- names(x)
  if (is.null(all_names)) {
    all_names <- as.character(seq_along(x))
  }
  if (is.null(keys)) {
    return(all_names)
  }
  if (is.numeric(keys)) {
    if (any(is.na(keys)) || any(keys < 1) || any(keys > length(x)) ||
        any(keys != as.integer(keys))) {
      stop(
        "'", arg, "' numeric indices must be integers in 1:", length(x), ".",
        call. = FALSE
      )
    }
    return(all_names[as.integer(keys)])
  }
  if (!is.character(keys) || length(keys) < 1L || anyNA(keys)) {
    stop(
      "'", arg, "' must be NULL, a character vector of ", what, "s, ",
      "or integer indices.",
      call. = FALSE
    )
  }
  keys <- as.character(keys)
  missing <- setdiff(keys, all_names)
  if (length(missing) > 0L) {
    stop(
      "'", arg, "' contains unknown ", what, "(s): ",
      paste(missing, collapse = ", "),
      ".\nAvailable: ", paste(all_names, collapse = ", "), ".",
      call. = FALSE
    )
  }
  keys
}
