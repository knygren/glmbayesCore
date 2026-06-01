#' Multi-response Bayesian linear models (\code{lmb})
#'
#' @description
#' Fits one \code{\link{lmb}} model per column of a multivariate response
#' (formula left-hand side with \code{cbind(...)}), sharing the same predictors
#' on the right-hand side. Returns a named list of \code{"lmb"} objects with
#' class \code{"mlmb"}.
#'
#' @details
#' This is the formula / \code{data} interface counterpart to
#' \code{\link{multi_rlmb}} (matrix \code{y}, \code{x}). Each response column
#' uses its own \code{pfamily_list[[j]]}. Use \code{\link{multi_prior_setup}}
#' to build aligned priors, then \code{\link{summary.mlmb}} or
#' \code{\link{print.mlmb}} for output styled like \code{\link[stats]{summary.mlm}}
#' (cf.\ \code{\link[stats]{lm}} with a matrix response via \code{cbind(...)}).
#'
#' @param formula A \code{\link{formula}} with a matrix response on the left-hand
#'   side (typically \code{cbind(...)}).
#' @param pfamily_list Named or unnamed list of length equal to the number of
#'   response columns; each element is a \code{\link{pfamily}} object for
#'   \code{\link{lmb}}.
#' @inheritParams lmb
#' @inheritParams multi_rlmb
#' @return A named list of class \code{"mlmb"}. Element \code{j} is an
#'   \code{"lmb"} fit for response \code{j}. Attributes include \code{call},
#'   \code{formula}, \code{coef_names}, \code{pred_names}, and
#'   \code{pfamily_lists}.
#' @seealso \code{\link{lmb}}, \code{\link{multi_rlmb}}, \code{\link{multi_prior_setup}},
#'   \code{\link{summary.mlmb}}, \code{\link{print.mlmb}},
#'   \code{\link[stats]{lm}} with \code{cbind} responses.
#' @family modelfuns
#' @example inst/examples/Ex_multi_lmb.R
#' @export
multi_lmb <- function(
    formula,
    pfamily_list,
    n = 1000,
    data,
    subset,
    weights,
    na.action,
    method = "qr",
    model = TRUE,
    x = TRUE,
    y = TRUE,
    qr = TRUE,
    singular.ok = TRUE,
    contrasts = NULL,
    offset,
    Gridtype = 2,
    n_envopt = NULL,
    use_parallel = TRUE,
    use_opencl = FALSE,
    verbose = FALSE,
    ...
) {
  mc <- match.call(expand.dots = FALSE)
  call <- match.call()
  if (missing(data)) {
    data <- environment(formula)
  }

  meta_args <- list(
    formula = formula,
    data = data,
    contrasts = contrasts,
    has_subset = !missing(subset),
    has_weights = !missing(weights),
    has_na.action = !missing(na.action),
    has_offset = !missing(offset)
  )
  if (!missing(subset)) meta_args$subset <- subset
  if (!missing(weights)) meta_args$weights <- weights
  if (!missing(na.action)) meta_args$na.action <- na.action
  if (!missing(offset)) meta_args$offset <- offset
  meta <- do.call(.mlmb_formula_meta, meta_args)

  pfamily_lists <- .mrglmb_normalize_pfamily_lists(
    pfamily_list,
    meta$l1,
    meta$p,
    .validate_pfamily_for_rlmb
  )

  lmb_args <- list(
    n = n,
    data = data,
    method = method,
    model = model,
    x = x,
    y = y,
    qr = qr,
    singular.ok = singular.ok,
    contrasts = contrasts,
    Gridtype = Gridtype,
    n_envopt = n_envopt,
    use_parallel = use_parallel,
    use_opencl = use_opencl,
    verbose = verbose
  )
  if (!missing(subset)) lmb_args$subset <- subset
  if (!missing(weights)) lmb_args$weights <- weights
  if (!missing(na.action)) lmb_args$na.action <- na.action
  if (!missing(offset)) lmb_args$offset <- offset
  if (length(list(...))) {
    lmb_args <- c(lmb_args, list(...))
  }

  block_results <- vector("list", meta$l1)
  for (j in seq_len(meta$l1)) {
    f_j <- stats::reformulate(meta$termlabels, response = meta$coef_names[j])
    fit_j <- do.call(
      lmb,
      c(list(formula = f_j, pfamily = pfamily_lists[[j]]), lmb_args)
    )
    fit_j$call <- .mlmb_lmb_display_call(mc, f_j)
    block_results[[j]] <- fit_j
  }

  .mlmb_assemble(
    block_results = block_results,
    coef_names = meta$coef_names,
    call = call,
    formula = formula,
    l1 = meta$l1,
    p = meta$p,
    pred_names = meta$pred_names,
    pfamily_lists = pfamily_lists
  )
}

#' @keywords internal
#' Build a short \code{lmb()} call for printing (omits \code{pfamily}).
.mlmb_lmb_display_call <- function(mc_multi, formula_j) {
  cl <- call("lmb", formula = formula_j)
  pass <- c(
    "n", "data", "subset", "weights", "na.action", "offset",
    "use_parallel", "use_opencl", "verbose", "Gridtype", "n_envopt"
  )
  for (nm in pass) {
    if (nm %in% names(mc_multi)) {
      cl[[nm]] <- mc_multi[[nm]]
    }
  }
  cl
}

#' @keywords internal
.mlmb_formula_meta <- function(
    formula,
    data,
    subset = NULL,
    weights = NULL,
    na.action = NULL,
    offset = NULL,
    contrasts = NULL,
    has_subset = FALSE,
    has_weights = FALSE,
    has_na.action = FALSE,
    has_offset = FALSE
) {
  mf_args <- list(
    formula = formula,
    data = data,
    drop.unused.levels = TRUE
  )
  if (has_subset) mf_args$subset <- subset
  if (has_weights) mf_args$weights <- weights
  if (has_na.action) mf_args$na.action <- na.action
  if (has_offset) mf_args$offset <- offset
  if (!is.null(contrasts)) mf_args$contrasts <- contrasts
  mf <- do.call(stats::model.frame, mf_args)
  mt <- attr(mf, "terms")
  Y <- as.matrix(stats::model.response(mf, "any"))
  l1 <- ncol(Y)
  if (l1 < 1L) {
    stop(
      "formula must specify at least one response column (e.g. cbind(...)).",
      call. = FALSE
    )
  }
  coef_names <- colnames(Y)
  if (is.null(coef_names) || length(coef_names) != l1) {
    coef_names <- paste0("Y", seq_len(l1))
  }

  x_mat <- stats::model.matrix(mt, mf, contrasts)
  p <- ncol(x_mat)
  pred_names <- colnames(x_mat)
  if (is.null(pred_names) || length(pred_names) != p) {
    pred_names <- paste0("X", seq_len(p))
  }

  list(
    coef_names = coef_names,
    termlabels = attr(mt, "term.labels"),
    l1 = l1,
    p = p,
    pred_names = pred_names
  )
}

#' @keywords internal
.mlmb_assemble <- function(
    block_results,
    coef_names,
    call,
    formula,
    l1,
    p,
    pred_names,
    pfamily_lists = NULL
) {
  outlist <- setNames(block_results, coef_names)
  attr(outlist, "call") <- call
  attr(outlist, "formula") <- formula
  attr(outlist, "l1") <- l1
  attr(outlist, "p") <- p
  attr(outlist, "coef_names") <- coef_names
  attr(outlist, "pred_names") <- pred_names
  if (!is.null(pfamily_lists)) {
    attr(outlist, "pfamily_lists") <- pfamily_lists
  }
  class(outlist) <- "mlmb"
  outlist
}

#' @keywords internal
.mlmb_coef_means_matrix <- function(object) {
  nm <- names(object)
  if (length(nm) < 1L) {
    return(NULL)
  }
  cm <- do.call(cbind, lapply(object, function(fit) fit$coef.means))
  rn <- names(object[[1L]]$coef.means)
  if (is.null(rn)) {
    rn <- colnames(object[[1L]]$x)
  }
  if (!is.null(rn) && nrow(cm) == length(rn)) {
    rownames(cm) <- rn
  }
  colnames(cm) <- nm
  cm
}

#' @keywords internal
.mlmb_dic_table <- function(object) {
  nm <- names(object)
  if (length(nm) < 1L) {
    return(NULL)
  }
  pD <- vapply(object, function(fit) fit$pD, numeric(1))
  dic <- vapply(object, function(fit) fit$DIC, numeric(1))
  cbind(pD = pD, DIC = dic)
}
