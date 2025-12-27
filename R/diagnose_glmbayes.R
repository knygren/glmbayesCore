#' @title Diagnose GPU and OpenCL Environment
#'
#' @description
#' Runs a full diagnostic pipeline to assess GPU availability, driver
#' installation, compute runtimes (CUDA/OpenCL), environment variable
#' configuration, and whether a usable OpenCL runtime platform exists.
#' A detailed human‑readable diagnostic report is printed to the console.
#'
#' The printed report includes:
#' \itemize{
#'   \item Operating system (Windows, Linux, WSL, macOS)
#'   \item GPU vendor selection (NVIDIA → AMD → Intel → CPU fallback)
#'   \item Driver installation status
#'   \item OpenCL header availability (\code{CL/cl.h})
#'   \item OpenCL runtime availability (ICD loader)
#'   \item Whether OpenCL is fully available (headers + runtime)
#'   \item PATH and library directory validation
#'   \item OpenCL runtime probe (Linux/WSL only)
#'   \item Whether \code{glmbayes} was compiled with OpenCL support
#' }
#'
#' The runtime probe attempts to create an OpenCL platform, device, context,
#' command queue, and compile a minimal kernel. On WSL, this probe typically
#' fails even when headers and the ICD loader are present; this is expected
#' because WSL does not expose OpenCL GPU devices.
#'
#' @return
#' Invisibly returns a named list containing:
#' \describe{
#'   \item{\code{environment_info}}{Output of \code{detect_environment_and_gpus()}}
#'   \item{\code{driver_status}}{Output of \code{detect_or_install_gpu_drivers()}}
#'   \item{\code{runtime_status}}{Output of \code{detect_compute_runtimes()}}
#'   \item{\code{env_diag}}{Output of \code{check_runtime_env()}}
#'   \item{\code{opencl_runtime_probe}}{TRUE/FALSE/NA depending on probe outcome}
#'   \item{\code{opencl_enabled}}{Whether \code{glmbayes} was compiled with OpenCL}
#' }
#'
#' @export
diagnose_glmbayes <- function() {
  cat("=== glmbayes Diagnostic Report ===\n")
  
  # Step 1: Environment + GPU detection
  info     <- detect_environment_and_gpus()
  drivers  <- detect_or_install_gpu_drivers(info)
  runtimes <- detect_compute_runtimes(info)
  env_diag <- check_runtime_env(runtimes)
  
  cat("Environment:", info$environment, "\n\n")
  
  # Step 2: Preference order (NVIDIA > AMD > Intel)
  gpu_vendor <- if (info$nvidia$present) "nvidia"
  else if (info$amd$present) "amd"
  else if (info$intel$present) "intel"
  else NULL
  
  if (!is.null(gpu_vendor)) {
    cat("GPU:", toupper(gpu_vendor), "\n")
    drv  <- drivers$drivers[[gpu_vendor]]
    rt   <- runtimes$runtimes[[gpu_vendor]]
    diag <- env_diag$diagnostics[[gpu_vendor]]
    
    # 1. Driver check
    if (drv$installed) {
      cat("  [OK] Driver installed\n")
    } else {
      cat("  [FAIL] Driver not installed\n")
      if (length(drv$issues) > 0)
        cat("    Issues:", paste(drv$issues, collapse=", "), "\n")
    }
    
    # 2. OpenCL header/runtime presence (NEW LOGIC)
    hdr <- rt$opencl$headers_present
    rtm <- rt$opencl$runtime_present
    inst <- rt$opencl$installed   # AND of both
    
    if (hdr) {
      cat("  [OK] OpenCL headers found (CL/cl.h)\n")
    } else {
      cat("  [FAIL] OpenCL headers not found (CL/cl.h missing)\n")
    }
    
    if (rtm) {
      cat("  [OK] OpenCL runtime found (OpenCL.dll / ICD)\n")
    } else {
      cat("  [FAIL] OpenCL runtime not found\n")
    }
    
    if (inst) {
      cat("  [OK] OpenCL fully available (headers + runtime)\n")
    } else {
      cat("  [FAIL] OpenCL incomplete (missing headers or runtime)\n")
    }
    
    # 3. PATH/lib environment validation
    paths_ok <- (length(diag$opencl$missing_path_dirs) == 0 &&
                   length(diag$opencl$missing_lib_dirs) == 0)
    
    if (paths_ok) {
      cat("  [OK] Required PATH and library dirs present\n")
    } else {
      if (length(diag$opencl$missing_path_dirs) > 0)
        cat("  [WARN] Missing PATH entries:",
            paste(diag$opencl$missing_path_dirs, collapse=", "), "\n")
      if (length(diag$opencl$missing_lib_dirs) > 0)
        cat("  [WARN] Missing library dirs:",
            paste(diag$opencl$missing_lib_dirs, collapse=", "), "\n")
    }
    
    # 4. Runtime probe (Linux/WSL only)
    runtime_ok <- NA
    if (paths_ok && tolower(info$environment) %in% c("linux", "wsl")) {
      runtime_ok <- verify_opencl_runtime(rt$opencl$lib_dirs)
      if (runtime_ok) {
        cat("  [OK] OpenCL runtime probe succeeded (platform available)\n")
      } else {
        cat("  [FAIL] OpenCL runtime probe failed (no usable platform)\n")
      }
    } else if (!paths_ok) {
      cat("  [SKIP] Runtime probe skipped (missing PATH/lib dirs)\n")
    } else {
      cat("  [SKIP] Runtime probe skipped on Windows\n")
    }
    
  } else {
    cat("[FAIL] No supported GPU detected. glmbayes will run in CPU-only mode.\n")
  }
  
  # Step 3: Report compile-time OpenCL status
  opencl_enabled <- NA
  if (exists("has_opencl") && is.function(has_opencl)) {
    opencl_enabled <- has_opencl()
    if (opencl_enabled) {
      cat("\n[OK] glmbayes was compiled with OpenCL support.\n")
    } else {
      cat("\n[FAIL] glmbayes was compiled without OpenCL support.\n")
    }
  }
  
  # Step 4: Interactive PATH/lib fixes
  missing_items <- (length(diag$opencl$missing_path_dirs) > 0 ||
                      length(diag$opencl$missing_lib_dirs) > 0)
  
  if (missing_items && !isTRUE(opencl_enabled)) {
    cat("\n[INFO] Missing PATH/lib entries detected and OpenCL is not enabled.\n")
    
    if (length(diag$opencl$missing_path_dirs) > 0) {
      cat("  Missing PATH entries:\n")
      cat("   -", paste(diag$opencl$missing_path_dirs, collapse="\n   - "), "\n")
      ans <- readline("Would you like to permanently add missing PATH dirs? [y/N]: ")
      if (tolower(ans) == "y") {
        if (tolower(info$environment) == "windows") {
          add_to_path_windows(diag$opencl$missing_path_dirs)
        } else {
          add_to_path_linux(diag$opencl$missing_path_dirs)
        }
      }
    }
    
    if (length(diag$opencl$missing_lib_dirs) > 0 &&
        tolower(info$environment) %in% c("linux", "wsl")) {
      cat("  Missing library dirs:\n")
      cat("   -", paste(diag$opencl$missing_lib_dirs, collapse="\n   - "), "\n")
      ans_lib <- readline("Would you like to permanently add missing library dirs to LD_LIBRARY_PATH? [y/N]: ")
      if (tolower(ans_lib) == "y") {
        add_to_libpath_linux(diag$opencl$missing_lib_dirs)
      }
    }
  }
  
  cat("\n=== End of Diagnostic Report ===\n")
  
  invisible(list(
    environment_info      = info,
    driver_status         = drivers,
    runtime_status        = runtimes,
    env_diag              = env_diag,
    opencl_runtime_probe  = runtime_ok,
    opencl_enabled        = opencl_enabled
  ))
}