# Row-group partition helpers for conditionally independent groupwise simulators.

#' Normalize a row-group partition
#'
#' Parse a group specification into a common partition layout used by
#' \code{\link{rNormal_reg_group}}, \code{\link{Prior_SetupGroup}}, and
#' row-block drivers in \pkg{lmebayes}.
#'
#' @param group Block partition: \code{factor} or integer vector of length
#'   \code{l2}, \code{l2_blocks} counts summing to \code{l2}, or a list of
#'   disjoint row-index vectors covering \code{1:l2}.
#' @param l2 Number of observations (rows) after \code{model.frame}.
#' @return A list with:
#'   \describe{
#'     \item{\code{k}}{Number of blocks.}
#'     \item{\code{ids}}{Character block labels (factor levels, list names, or
#'       \code{"block1"}, \ldots).}
#'     \item{\code{l2_blocks}}{Integer vector of length \code{k}: observations
#'       per block.}
#'     \item{\code{starts}}{Integer start index of each block in a contiguous
#'       stacking of block rows (1-based).}
#'     \item{\code{rows}}{List of length \code{k}: integer row indices in
#'       \code{1:l2} for each block.}
#'   }
#' @export
#' @example inst/examples/Ex_normalize_group.R
normalize_group <- function(group, l2) {
  l2 <- as.integer(l2)
  if (length(l2) != 1L || l2 < 1L) {
    stop("'l2' must be a positive integer (length of y).", call. = FALSE)
  }

  if (is.list(group)) {
    if (length(group) < 1L) {
      stop("'group' list must have at least one element.", call. = FALSE)
    }
    rows <- lapply(group, function(idx) {
      idx <- as.integer(idx)
      if (anyNA(idx) || any(idx < 1L) || any(idx > l2)) {
        stop("Row indices in 'group' must be integers in 1:l2.", call. = FALSE)
      }
      unique(idx)
    })
    all_idx <- unlist(rows, use.names = FALSE)
    if (any(duplicated(all_idx))) {
      stop("Row indices in 'group' list must be disjoint.", call. = FALSE)
    }
    if (!identical(sort(all_idx), seq_len(l2))) {
      stop("Row indices in 'group' list must cover exactly 1:l2.", call. = FALSE)
    }
    k <- length(rows)
    ids <- names(group)
    if (is.null(ids) || any(ids == "")) {
      ids <- paste0("group", seq_len(k))
    }
    l2_blocks <- vapply(rows, length, integer(1L))
    starts <- c(1L, cumsum(l2_blocks)[-k] + 1L)
    return(list(
      k = k,
      ids = ids,
      l2_blocks = l2_blocks,
      starts = starts,
      rows = rows
    ))
  }

  # If group is already a factor, preserve its level order before any coercion.
  # as.vector() on a factor returns integer indices, and re-calling factor() on
  # those integers sorts the labels lexicographically, silently destroying the
  # caller-supplied level order.
  if (is.factor(group) && length(group) == l2) {
    blk <- group
    k <- nlevels(blk)
    rows <- split(seq_len(l2), blk)
    ids <- levels(blk)
    l2_blocks <- vapply(rows, length, integer(1L))
    starts <- c(1L, cumsum(l2_blocks)[-k] + 1L)
    return(list(
      k = k,
      ids = ids,
      l2_blocks = l2_blocks,
      starts = starts,
      rows = rows
    ))
  }

  group <- as.vector(group)

  if (length(group) == l2) {
    blk <- factor(group)
    k <- nlevels(blk)
    rows <- split(seq_len(l2), blk)
    ids <- levels(blk)
    l2_blocks <- vapply(rows, length, integer(1L))
    starts <- c(1L, cumsum(l2_blocks)[-k] + 1L)
    return(list(
      k = k,
      ids = ids,
      l2_blocks = l2_blocks,
      starts = starts,
      rows = rows
    ))
  }

  if (length(group) >= 1L && length(group) < l2 && all(group >= 1L)) {
    l2_blocks <- as.integer(group)
    if (sum(l2_blocks) != l2) {
      stop(
        "When 'group' has length k < l2, it is treated as l2_blocks; ",
        "sum(group) must equal length(y) (", l2, ").",
        call. = FALSE
      )
    }
    k <- length(l2_blocks)
    ends <- cumsum(l2_blocks)
    starts <- c(1L, ends[-k] + 1L)
    rows <- lapply(seq_len(k), function(j) seq.int(starts[j], ends[j]))
    ids <- paste0("group", seq_len(k))
    return(list(
      k = k,
      ids = ids,
      l2_blocks = l2_blocks,
      starts = starts,
      rows = rows
    ))
  }

  stop(
    "'group' must be a factor or integer vector of length l2, ",
    "a list of row indices, or an integer vector of l2_blocks counts.",
    call. = FALSE
  )
}
