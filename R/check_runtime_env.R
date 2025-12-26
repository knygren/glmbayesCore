#' Check GPU Runtime Environment
#'
#' Compares the directories detected by \code{detect_compute_runtimes()} against
#' the current PATH and LD_LIBRARY_PATH. Reports what is present and what is missing,
#' and suggests what should be added to ensure GPU runtimes are found correctly.
#'
#' This function is environment-agnostic: it works across Linux, WSL, and Windows/MSYS2
#' because it consumes the structured output of \code{detect_compute_runtimes()}.
#'
#' @param runtime_info The structured list returned by \code{detect_compute_runtimes()}.
#' @return A structured list with diagnostics for each runtime, including:
#'   \itemize{
#'     \item installed: whether the runtime was detected
#'     \item found_path_dirs: directories already present in PATH
#'     \item missing_path_dirs: directories that should be added to PATH
#'     \item found_lib_dirs: directories already present in LD_LIBRARY_PATH
#'     \item missing_lib_dirs: directories that should be added to LD_LIBRARY_PATH
#'     \item include_dirs: include directories detected for headers
#'   }
#' @examples
#' info <- detect_compute_runtimes(list(environment="wsl"))
#' diag <- check_runtime_env(info)
#' str(diag)
#' @export
check_runtime_env <- function(runtime_info) {
  # Helper: normalize paths for comparison
  normalize_for_compare <- function(p) {
    if (length(p) == 0) return(character(0))
    vapply(p, function(x) {
      tryCatch(
        tolower(normalizePath(x, winslash = "/", mustWork = FALSE)),
        error = function(e) tolower(gsub("\\\\", "/", x))
      )
    }, character(1))
  }
  
  env <- runtime_info$environment
  
  # Current environment variables
  current_path <- strsplit(Sys.getenv("PATH"),
                           ifelse(.Platform$OS.type == "windows", ";", ":"))[[1]]
  current_path <- normalize_for_compare(current_path)
  
  current_lib <- strsplit(Sys.getenv("LD_LIBRARY_PATH", unset=""), ":")[[1]]
  current_lib <- normalize_for_compare(current_lib)
  
  results <- list(environment = env, diagnostics = list())
  
  for (vendor in names(runtime_info$runtimes)) {
    vendor_info <- runtime_info$runtimes[[vendor]]
    results$diagnostics[[vendor]] <- list()
    
    for (runtime in names(vendor_info)) {
      rt <- vendor_info[[runtime]]
      
      expected_bin <- normalize_for_compare(rt$bin_dirs)
      expected_lib <- normalize_for_compare(rt$lib_dirs)
      
      # PATH comparison
      found_path <- intersect(expected_bin, current_path)
      missing_path <- setdiff(expected_bin, current_path)
      
      # Library comparison depends on environment
      if (env %in% c("windows", "msys2")) {
        # On Windows/MSYS2, LD_LIBRARY_PATH is not used; check existence
        found_lib <- Filter(dir.exists, rt$lib_dirs)
        missing_lib <- setdiff(rt$lib_dirs, found_lib)
      } else {
        # On Linux/WSL/macOS, compare against LD_LIBRARY_PATH
        found_lib <- intersect(expected_lib, current_lib)
        missing_lib <- setdiff(expected_lib, current_lib)
      }
      
      rt_result <- list(
        installed = rt$installed,
        found_path_dirs   = found_path,
        missing_path_dirs = missing_path,
        found_lib_dirs    = found_lib,
        missing_lib_dirs  = missing_lib,
        include_dirs      = rt$include_dirs
      )
      
      results$diagnostics[[vendor]][[runtime]] <- rt_result
    }
  }
  
  return(results)
}