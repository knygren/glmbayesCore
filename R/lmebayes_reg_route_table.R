## Routing table for the Gaussian LMM matrix engines rlmerb() calls.
## The non-Gaussian glmm_* routes stay in lmebayesCore with rglmerb().

#' Route key for REG_ROUTE_TABLE lookup
#' @noRd
.lmebayes_reg_route_key <- function(family, disp_mode, any_non_normal) {
  is_gaussian <- is.null(family) || identical(family$family, "gaussian")
  if (!is_gaussian) {
    stop(
      "rlmerb() only fits gaussian() responses. Non-Gaussian families ",
      "use rglmerb() in lmebayesCore.",
      call. = FALSE
    )
  }
  est_suffix <- if (isTRUE(any_non_normal)) "estimated" else "known"
  if (disp_mode %in% c("gamma", "gamma_list")) {
    paste0("lmm_gamma_", est_suffix)
  } else {
    paste0("lmm_fixed_", est_suffix)
  }
}

#' Named route metadata (export symbol + pilot / draw labels).
#' @noRd
REG_ROUTE_TABLE <- list(
  lmm_fixed_known = list(
    export            = "rLMMNormal_reg_known_vcov",
    needs_pilot       = FALSE,
    draw_engine_label = "rGLMM_sweep"
  ),
  lmm_fixed_estimated = list(
    export            = "rLMMNormal_reg_estimated_vcov",
    needs_pilot       = TRUE,
    draw_engine_label = "rGLMM_sweep"
  ),
  lmm_gamma_known = list(
    export            = "rLMMindepNormalGamma_reg_known_vcov",
    needs_pilot       = FALSE,
    draw_engine_label = "rGLMM_sweep"
  ),
  lmm_gamma_estimated = list(
    export            = "rLMMindepNormalGamma_reg_estimated_vcov",
    needs_pilot       = TRUE,
    draw_engine_label = "rGLMM_sweep"
  )
)

#' Resolve a REG_ROUTE_TABLE entry to a callable export
#' @noRd
.lmebayes_reg_route_fn <- function(route_key) {
  entry <- REG_ROUTE_TABLE[[route_key]]
  if (is.null(entry)) {
    stop("Unknown reg route key: ", route_key, call. = FALSE)
  }
  fn <- get(
    entry$export,
    mode     = "function",
    envir    = asNamespace("glmbayesCore"),
    inherits = FALSE
  )
  list(
    export_fn         = fn,
    needs_pilot       = entry$needs_pilot,
    draw_engine_label = entry$draw_engine_label
  )
}
