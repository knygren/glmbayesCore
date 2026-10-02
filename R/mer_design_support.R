## Minimal helpers for extract_re_hyper_matrices() (Step 1 MER migration).
## Full model_setup() helpers remain in lmebayesCore until later steps.

#' @noRd
.lmebayes_normalize_weights <- function(weights, n, arg = "weights") {
  if (is.null(weights)) {
    return(rep(1, n))
  }
  if (!is.numeric(weights)) {
    stop(sprintf("'%s' must be numeric.", arg), call. = FALSE)
  }
  wt <- as.numeric(weights)
  if (length(wt) == 1L) {
    wt <- rep(wt, n)
  }
  if (length(wt) != n) {
    stop(
      sprintf(
        "'%s' must be scalar or length %d (number of observations).",
        arg, n
      ),
      call. = FALSE
    )
  }
  if (anyNA(wt) || any(wt < 0)) {
    stop(sprintf("'%s' must be nonnegative and non-missing.", arg), call. = FALSE)
  }
  wt
}

#' @noRd
.lmebayes_normalize_offset <- function(offset, n, arg = "offset") {
  if (is.null(offset)) {
    return(rep(0, n))
  }
  if (!is.numeric(offset)) {
    stop(sprintf("'%s' must be numeric.", arg), call. = FALSE)
  }
  off <- as.numeric(offset)
  if (length(off) == 1L) {
    off <- rep(off, n)
  }
  if (length(off) != n) {
    stop(
      sprintf(
        "'%s' must be scalar or length %d (number of observations).",
        arg, n
      ),
      call. = FALSE
    )
  }
  if (anyNA(off)) {
    stop(sprintf("'%s' must be non-missing.", arg), call. = FALSE)
  }
  off
}

#' @noRd
.lmebayes_weights_offset_from_frame <- function(fr, n) {
  list(
    weights = .lmebayes_normalize_weights(stats::model.weights(fr), n),
    offset  = .lmebayes_normalize_offset(stats::model.offset(fr), n)
  )
}
