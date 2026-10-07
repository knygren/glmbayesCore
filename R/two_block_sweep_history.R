#' Build structured two-block sweep history from per-sweep fixef snapshots
#' @param stage_label Character stage label (e.g. \code{"main"}).
#' @param sweep_stats Nested list from \code{.two_block_snapshot_fixef_stats}.
#' @param fixef_mode Named list of ICM mode vectors per RE component.
#' @param re_names Character vector of RE component names.
#' @param sweep_cov Optional list (length \code{length(sweep_stats)}) of
#'   \code{.two_block_snapshot_fixef_cov()} results, one per inner sweep. When
#'   supplied, stored as \code{cov_by_sweep} (list of covariance matrices) and
#'   \code{coef_index} (data frame of stacking order) on the returned object,
#'   enabling \code{\link{plot_var_convergence}}'s whitened mode.
#' @param disp_stats Optional list (length \code{length(sweep_stats)}) of
#'   \code{\link{.two_block_snapshot_disp_stats}} results, one per inner
#'   sweep. When supplied (and non-empty), stored as \code{disp_table} --
#'   same 5-column shape as \code{table} (\code{covariate} is always
#'   \code{"precision"}) -- enabling \code{\link{plot_var_convergence}}'s
#'   \code{component = "precision"} mode. \code{NULL}/empty (known-vcov
#'   engines, or zero ING Block~2 components) yields a 0-row \code{disp_table}.
#' @return Object of class \code{"two_block_sweep_history"}.
#' @noRd
.two_block_build_sweep_history <- function(
    stage_label,
    sweep_stats,
    fixef_mode,
    re_names,
    sweep_cov = NULL,
    disp_stats = NULL
) {
  stage_label <- as.character(stage_label)[1L]
  if (!nzchar(stage_label)) {
    stage_label <- "stage"
  }
  n_sweep <- length(sweep_stats)
  rows <- list()

  for (k in re_names) {
    cn <- if (n_sweep >= 1L) {
      names(sweep_stats[[1L]][[k]]$mean)
    } else {
      character(0L)
    }
    if (is.null(cn) || !length(cn)) {
      cn <- names(fixef_mode[[k]])
    }
    if (is.null(cn)) {
      cn <- character(0L)
    }
    for (nm in cn) {
      mode_v <- .two_block_fixef_mode_at(fixef_mode, k, nm, cn)
      rows[[length(rows) + 1L]] <- data.frame(
        re_component = k,
        covariate    = nm,
        sweep        = 0L,
        mean         = mode_v,
        sd           = NA_real_,
        stringsAsFactors = FALSE
      )
      for (m in seq_len(n_sweep)) {
        rows[[length(rows) + 1L]] <- data.frame(
          re_component = k,
          covariate    = nm,
          sweep        = as.integer(m),
          mean         = sweep_stats[[m]][[k]]$mean[[nm]],
          sd           = sweep_stats[[m]][[k]]$sd[[nm]],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  table <- if (length(rows)) {
    do.call(rbind, rows)
  } else {
    data.frame(
      re_component = character(),
      covariate    = character(),
      sweep        = integer(),
      mean         = numeric(),
      sd           = numeric(),
      stringsAsFactors = FALSE
    )
  }

  cov_by_sweep <- NULL
  coef_index   <- NULL
  if (!is.null(sweep_cov) && length(sweep_cov)) {
    cov_by_sweep <- lapply(sweep_cov, function(sc) sc$cov)
    coef_index   <- sweep_cov[[1L]]$coef_index
  }

  disp_ing_names <- if (!is.null(disp_stats) && length(disp_stats)) {
    names(disp_stats[[1L]])
  } else {
    character(0L)
  }
  disp_rows <- list()
  for (k in disp_ing_names) {
    for (m in seq_len(n_sweep)) {
      disp_rows[[length(disp_rows) + 1L]] <- data.frame(
        re_component = k,
        covariate    = "precision",
        sweep        = as.integer(m),
        mean         = disp_stats[[m]][[k]]$mean,
        sd           = disp_stats[[m]][[k]]$sd,
        stringsAsFactors = FALSE
      )
    }
  }
  disp_table <- if (length(disp_rows)) {
    do.call(rbind, disp_rows)
  } else {
    data.frame(
      re_component = character(),
      covariate    = character(),
      sweep        = integer(),
      mean         = numeric(),
      sd           = numeric(),
      stringsAsFactors = FALSE
    )
  }

  structure(
    list(
      stage        = stage_label,
      n_sweeps     = as.integer(n_sweep),
      re_names     = re_names,
      fixef_mode   = fixef_mode,
      table        = table,
      cov_by_sweep = cov_by_sweep,
      coef_index   = coef_index,
      disp_table   = disp_table
    ),
    class = c("two_block_sweep_history", "list")
  )
}

#' Filter sweep-history table rows for printing
#' @noRd
.two_block_filter_sweep_history_table <- function(
    table,
    max_sweeps = Inf,
    sweeps = NULL,
    components = NULL,
    covariate = NULL
) {
  if (is.null(table) || !nrow(table)) {
    return(table)
  }
  out <- table
  if (!is.null(components)) {
    components <- as.character(components)
    out <- out[out$re_component %in% components, , drop = FALSE]
  }
  if (!is.null(covariate)) {
    covariate <- as.character(covariate)
    out <- out[out$covariate %in% covariate, , drop = FALSE]
  }
  if (!is.null(sweeps)) {
    sweeps <- as.integer(sweeps)
    out <- out[out$sweep %in% c(0L, sweeps), , drop = FALSE]
  } else if (is.finite(max_sweeps) && max_sweeps >= 0L) {
    max_s <- max(out$sweep, na.rm = TRUE)
    keep <- out$sweep == 0L |
      out$sweep >= max(1L, max_s - as.integer(max_sweeps) + 1L)
    out <- out[keep, , drop = FALSE]
  }
  out
}

#' Sweep count shown in the print header
#' @noRd
.two_block_sweep_history_header_n <- function(x, sweeps = NULL, tab = NULL) {
  if (!is.null(sweeps) && length(sweeps) == 1L) {
    return(as.integer(sweeps[1L]))
  }
  if (!is.null(tab) && nrow(tab)) {
    sweep_rows <- tab$sweep[tab$sweep > 0L]
    if (length(sweep_rows)) {
      return(max(sweep_rows, na.rm = TRUE))
    }
  }
  x$n_sweeps
}

#' Print one Block~2 sweep-history table (mode + optional sweep rows)
#' @noRd
.two_block_print_sweep_history_body <- function(
    x,
    tab,
    sweeps = NULL,
    digits = 4L
) {
  if (!nrow(tab)) {
    cat(sprintf(
      "\n--- two-block [%s stage summary: fixef by sweep (0 sweeps)] ---\n\n",
      x$stage
    ))
    return(invisible(NULL))
  }

  n_hdr <- .two_block_sweep_history_header_n(x, sweeps = sweeps, tab = tab)
  cat(sprintf(
    "\n--- two-block [%s stage summary: fixef by sweep (%d sweeps)] ---\n",
    x$stage, n_hdr
  ))
  cat("  Block 2 fixed effects (mode and chain stats after each sweep):\n")
  hdr <- sprintf(
    "  %-18s  %-30s  %12s  %12s  %12s",
    "Random effect", "Covariate", "mode/sweep", "mean", "sd"
  )
  sep <- paste0("  ", strrep("-", nchar(hdr) - 2L))
  cat(hdr, "\n")
  cat(sep, "\n")

  row_keys <- unique(tab[, c("re_component", "covariate"), drop = FALSE])
  for (i in seq_len(nrow(row_keys))) {
    k  <- row_keys$re_component[i]
    nm <- row_keys$covariate[i]
    sub <- tab[tab$re_component == k & tab$covariate == nm, , drop = FALSE]
    sub <- sub[order(sub$sweep), , drop = FALSE]
    for (j in seq_len(nrow(sub))) {
      sweep_label <- if (sub$sweep[j] == 0L) {
        "mode"
      } else {
        paste0("sweep ", sub$sweep[j])
      }
      sd_str <- if (is.na(sub$sd[j])) {
        ""
      } else {
        sprintf("%12.*f", digits, sub$sd[j])
      }
      cat(sprintf(
        "  %-18s  %-30s  %12s  %12.*f  %s\n",
        k, nm, sweep_label, digits, sub$mean[j], sd_str
      ))
    }
  }
  cat("\n")
  invisible(NULL)
}

#' Print one Block~2 RE-precision (\code{1/tau^2}) sweep-history table
#'
#' Companion to \code{\link{.two_block_print_sweep_history_body}} for
#' \code{x$disp_table} (see \code{\link{.two_block_snapshot_disp_stats}});
#' a no-op when \code{tab} has no rows (known-vcov stages, or an
#' estimated-vcov stage with zero ING Block~2 components).
#' @noRd
.two_block_print_disp_history_body <- function(x, tab, digits = 4L) {
  if (is.null(tab) || !nrow(tab)) {
    return(invisible(NULL))
  }
  cat("  Block 2 RE precision (1/tau^2; chain mean and SD after each sweep):\n")
  hdr <- sprintf("  %-18s  %12s  %12s  %12s", "Random effect", "sweep", "mean", "sd")
  sep <- paste0("  ", strrep("-", nchar(hdr) - 2L))
  cat(hdr, "\n")
  cat(sep, "\n")

  row_keys <- unique(tab$re_component)
  for (k in row_keys) {
    sub <- tab[tab$re_component == k, , drop = FALSE]
    sub <- sub[order(sub$sweep), , drop = FALSE]
    for (j in seq_len(nrow(sub))) {
      cat(sprintf(
        "  %-18s  %12s  %12.*f  %12.*f\n",
        k, paste0("sweep ", sub$sweep[j]), digits, sub$mean[j], digits, sub$sd[j]
      ))
    }
  }
  cat("\n")
  invisible(NULL)
}

#' Two-block sweep history object
#' @name two_block_sweep_history
#' @keywords internal
NULL

#' Print two-block sweep history (fixef mode + per-sweep chain stats)
#'
#' @param x Object of class \code{"two_block_sweep_history"}.
#' @param max_sweeps When \code{sweeps} is \code{NULL}, show mode rows plus
#'   the last \code{max_sweeps} inner sweeps.
#' @param sweeps Optional integer vector of sweep indices to include (mode
#'   rows are always retained).
#' @param components Optional character vector of RE components to include.
#' @param covariate Optional character vector of covariate names to include
#'   (use with \code{components} to trace one coefficient at a time).
#' @param by_sweep When \code{TRUE}, print a separate table for each inner
#'   sweep (mode row plus that sweep only), matching live \code{glmerb()} /
#'   \code{diag_sweeps} output. When \code{FALSE} (default), all selected
#'   sweeps appear in one table per coefficient.
#' @param digits Number of decimal places for numeric columns.
#' @param ... Ignored.
#' @importFrom utils tail
#' @return \code{x} invisibly.
#' @rdname two_block_sweep_history
#' @method print two_block_sweep_history
#' @export
print.two_block_sweep_history <- function(
    x,
    max_sweeps = Inf,
    sweeps = NULL,
    components = NULL,
    covariate = NULL,
    by_sweep = FALSE,
    digits = 4L,
    ...
) {
  if (isTRUE(by_sweep)) {
    all_sweeps <- seq_len(x$n_sweeps)
    if (!is.null(sweeps)) {
      use <- as.integer(sweeps)
    } else if (is.finite(max_sweeps)) {
      use <- tail(all_sweeps, as.integer(max_sweeps))
    } else {
      use <- all_sweeps
    }
    use <- use[use >= 1L & use <= x$n_sweeps]
    for (m in use) {
      tab_m <- .two_block_filter_sweep_history_table(
        x$table,
        sweeps = m,
        components = components,
        covariate = covariate
      )
      .two_block_print_sweep_history_body(x, tab_m, sweeps = m, digits = digits)
      disp_tab_m <- .two_block_filter_sweep_history_table(
        x$disp_table,
        sweeps = m,
        components = components
      )
      .two_block_print_disp_history_body(x, disp_tab_m, digits = digits)
    }
    return(invisible(x))
  }

  tab <- .two_block_filter_sweep_history_table(
    x$table,
    max_sweeps = max_sweeps,
    sweeps = sweeps,
    components = components,
    covariate = covariate
  )
  .two_block_print_sweep_history_body(x, tab, sweeps = sweeps, digits = digits)
  disp_tab <- .two_block_filter_sweep_history_table(
    x$disp_table,
    max_sweeps = max_sweeps,
    sweeps = sweeps,
    components = components
  )
  .two_block_print_disp_history_body(x, disp_tab, digits = digits)
  invisible(x)
}

