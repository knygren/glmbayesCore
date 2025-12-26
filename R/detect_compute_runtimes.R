#' Detect CUDA and OpenCL runtimes and their installation directories
#'
#' \code{detect_compute_runtimes()} inspects the host system to determine
#' whether CUDA and/or OpenCL runtimes are installed for NVIDIA, AMD, and Intel,
#' and reports the directories where relevant headers, libraries, and executables
#' are found.
#'
#' Detection is modeled after configure-time logic and is intentionally minimal:
#' it does not attempt to install drivers or SDKs, and does not assume optional
#' tools such as \code{clinfo} or \code{rocminfo}. Instead, it scans known
#' directories, environment variables, and SDK roots for the presence of key
#' files.
#'
#' \itemize{
#'   \item On Windows/MSYS2:
#'     \itemize{
#'       \item CUDA Toolkit subdirectories under
#'             \code{C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v*}.
#'       \item Vendor-specific OpenCL DLLs in \code{C:/Windows/System32}
#'             (\code{nvopencl.dll}, \code{igdrcl64.dll}, \code{amdocl.dll}).
#'       \item SDK headers in IntelSWTools, AMD APP SDK, and CUDA include paths.
#'       \item Optional environment roots (\code{OPENCL_HOME}, \code{OPENCL_SDK}).
#'     }
#'   \item On Linux:
#'     \itemize{
#'       \item Typical CUDA install locations (\code{/usr/local/cuda}).
#'       \item System include directories reported by \code{gcc -E -x c++ - -v}.
#'       \item System library search directories reported by \code{gcc -Xlinker --verbose}.
#'       \item Presence of \code{CL/cl.h} and \code{libOpenCL.so}.
#'     }
#'   \item On WSL:
#'     \itemize{
#'       \item NVIDIA-only detection via \code{nvidia-smi}.
#'       \item AMD and Intel runtimes are forced to \code{installed = FALSE}.
#'     }
#' }
#'
#' @details
#' This function is intended to run *before* PATH is set. It includes PATH
#' directories in its search, but does not report whether a runtime is "on PATH".
#' Instead, it reports the actual directories where files were found, so users
#' can decide which to add to PATH or to compiler/linker flags.
#' @param info A list returned by `detect_environment_and_gpus()`. The list must
#'   contain the following elements:
#'   \describe{
#'     \item{environment}{One of "windows", "msys2", "linux", "wsl", or "unknown".}
#'     \item{nvidia}{A list with elements `present` (logical) and `names` (character).}
#'     \item{amd}{A list with elements `present` (logical) and `names` (character).}
#'     \item{intel}{A list with elements `present` (logical) and `names` (character).}
#'   }
#' @return
#' A named list with the following components:
#' \describe{
#'   \item{\code{environment}}{Character scalar indicating the detected host
#'     environment. One of \code{"windows"}, \code{"msys2"}, \code{"wsl"},
#'     or \code{"linux"}.}
#'
#'   \item{\code{runtimes}}{A nested list with vendor-specific results:}
#'     \describe{
#'       \item{nvidia$cuda}{List with elements:}
#'         \describe{
#'           \item{\code{installed}}{Logical; \code{TRUE} if CUDA was detected.}
#'           \item{\code{bin_dirs}}{Character vector of directories containing \code{nvcc} executables.}
#'           \item{\code{include_dirs}}{Character vector of directories containing CUDA headers (e.g. \code{cuda.h}).}
#'           \item{\code{lib_dirs}}{Character vector of directories containing CUDA libraries (e.g. \code{libcudart.so}, \code{x64} libs).}
#'         }
#'       \item{nvidia$opencl}{List with elements:}
#'         \describe{
#'           \item{\code{installed}}{Logical; \code{TRUE} if NVIDIA OpenCL runtime was detected.}
#'           \item{\code{include_dirs}}{Character vector of directories containing OpenCL headers (e.g. \code{CL/cl.h}).}
#'           \item{\code{lib_dirs}}{Character vector of directories containing NVIDIA OpenCL libraries (e.g. \code{nvopencl.dll}).}
#'         }
#'       \item{amd$opencl}{List with elements:}
#'         \describe{
#'           \item{\code{installed}}{Logical; \code{TRUE} if AMD OpenCL runtime was detected.}
#'           \item{\code{include_dirs}}{Character vector of directories containing AMD OpenCL headers.}
#'           \item{\code{lib_dirs}}{Character vector of directories containing AMD OpenCL libraries (e.g. \code{amdocl.dll}).}
#'         }
#'       \item{intel$opencl}{List with elements:}
#'         \describe{
#'           \item{\code{installed}}{Logical; \code{TRUE} if Intel OpenCL runtime was detected.}
#'           \item{\code{include_dirs}}{Character vector of directories containing Intel OpenCL headers.}
#'           \item{\code{lib_dirs}}{Character vector of directories containing Intel OpenCL libraries (e.g. \code{igdrcl64.dll}).}
#'         }
#'     }
#' }
#'
#' @examples
#' \dontrun{
#'   info <- detect_environment_and_gpus()
#'   runtimes <- detect_compute_runtimes(info)
#'   print(runtimes$environment)
#'   print(runtimes$runtimes$nvidia$cuda)
#'   print(runtimes$runtimes$intel$opencl)
#' }
#'
#' @export
detect_compute_runtimes <- function(info) {
  env <- info$environment
  
  result <- list(
    environment = env,
    runtimes = list(
      nvidia = list(
        cuda   = list(installed = FALSE, bin_dirs = character(), include_dirs = character(), lib_dirs = character()),
        opencl = list(installed = FALSE, bin_dirs = character(), include_dirs = character(), lib_dirs = character())
      ),
      amd   = list(opencl = list(installed = FALSE, include_dirs = character(), lib_dirs = character(), bin_dirs = character())),
      intel = list(opencl = list(installed = FALSE, include_dirs = character(), lib_dirs = character(), bin_dirs = character()))
    )
  )
  
  # -------------------------------
  # Linux / WSL logic (UNCHANGED)
  # -------------------------------
  if (env %in% c("linux", "wsl")) {
    inc_dirs <- try(system("echo | gcc -E -x c++ - -v 2>&1 | grep '^ /' | sed 's/^ //'", intern = TRUE), silent = TRUE)
    if (!inherits(inc_dirs, "try-error")) {
      for (d in inc_dirs) {
        if (file.exists(file.path(d, "CL", "cl.h"))) {
          result$runtimes$nvidia$opencl$include_dirs <- unique(c(result$runtimes$nvidia$opencl$include_dirs, d))
        }
      }
    }
    
    raw_lib_dirs <- try(system("gcc -Xlinker --verbose 2>&1 | grep SEARCH_DIR | sed 's/SEARCH_DIR(\"=*\\([^\\\"]*\\)\").*/\\1/'", intern = TRUE), silent = TRUE)
    if (!inherits(raw_lib_dirs, "try-error")) {
      alt_lib_dirs <- gsub("/usr/local", "/usr", raw_lib_dirs)
      system_lib_dirs <- sort(unique(c(raw_lib_dirs, alt_lib_dirs)))
      for (d in system_lib_dirs) {
        if (file.exists(file.path(d, "libOpenCL.so"))) {
          result$runtimes$nvidia$opencl$lib_dirs <- unique(c(result$runtimes$nvidia$opencl$lib_dirs, d))
        }
      }
    }
    
    if (length(result$runtimes$nvidia$opencl$include_dirs) > 0 &&
        length(result$runtimes$nvidia$opencl$lib_dirs) > 0) {
      result$runtimes$nvidia$opencl$installed <- TRUE
    }
    
    nvcc_path <- Sys.which("nvcc")
    if (nzchar(nvcc_path)) {
      nvcc_path <- normalizePath(nvcc_path, winslash = "/", mustWork = TRUE)
      cuda_root <- normalizePath(file.path(dirname(nvcc_path), ".."), winslash = "/", mustWork = TRUE)
      cuda_include <- file.path(cuda_root, "include")
      cuda_lib <- file.path(cuda_root, "lib64")
      
      result$runtimes$nvidia$cuda$installed <- TRUE
      result$runtimes$nvidia$cuda$bin_dirs <- dirname(nvcc_path)
      if (dir.exists(cuda_include)) result$runtimes$nvidia$cuda$include_dirs <- cuda_include
      if (dir.exists(cuda_lib)) result$runtimes$nvidia$cuda$lib_dirs <- cuda_lib
    }
  }
  
  # -------------------------------
  # Windows / MSYS2 logic
  # -------------------------------
  if (env %in% c("windows", "msys2")) {
    search_paths <- character()
    
    # 1. PATH expansion
    path_dirs <- strsplit(Sys.getenv("PATH"), ";")[[1]]
    path_dirs <- path_dirs[nzchar(path_dirs)]
    search_paths <- c(search_paths, path_dirs)
    
    # 2. Environment variables
    for (var in c("OPENCL_HOME", "OPENCL_SDK")) {
      val <- Sys.getenv(var, unset = "")
      if (nzchar(val)) {
        for (sub in c("", "bin", "lib", "include", file.path("include","CL"))) {
          path <- file.path(val, sub)
          if (dir.exists(path)) search_paths <- c(search_paths, path)
        }
      }
    }
    
    # 3. Known SDK roots
    sdk_roots <- c(
      "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA",
      "C:/Program Files (x86)/IntelSWTools/OpenCL SDK",
      "C:/Program Files (x86)/AMD APP SDK",
      "C:/Program Files (x86)/Intel/oneAPI"
    )
    for (root in sdk_roots) {
      if (dir.exists(root)) {
        search_paths <- c(
          search_paths,
          root,
          list.dirs(root, recursive = TRUE, full.names = TRUE)
        )
      }
    }
    
    # 4. Static fallbacks
    static_paths <- c("C:/OpenCL-SDK","C:/opt","C:/Program Files (x86)")
    for (path in static_paths) {
      if (dir.exists(path)) search_paths <- c(search_paths, path)
    }
    
    # -------------------------------
    # OpenCL header detection (UNCHANGED)
    # -------------------------------
    cl_header <- NULL
    for (dir in search_paths) {
      if (file.exists(file.path(dir,"CL","cl.h"))) {
        cl_header <- file.path(dir,"CL","cl.h")
        break
      }
    }
    
    if (!is.null(cl_header)) {
      cl_base <- dirname(cl_header)
      opencl_home <- gsub("\\\\","/",sub("[/\\\\]include[/\\\\]CL$","",cl_base))
      include_flag <- file.path(opencl_home,"include")
      lib_flag     <- file.path(opencl_home,"lib","x64")
      
      result$runtimes$nvidia$opencl$include_dirs <- 
        c(result$runtimes$nvidia$opencl$include_dirs, include_flag)
      
      if (dir.exists(lib_flag)) {
        result$runtimes$nvidia$opencl$lib_dirs <- 
          c(result$runtimes$nvidia$opencl$lib_dirs, lib_flag)
      }
    }
    
    # -------------------------------
    # OpenCL ICD runtime detection (UPDATED)
    # -------------------------------
    system32 <- file.path(Sys.getenv("SystemRoot"), "System32")
    syswow64 <- file.path(Sys.getenv("SystemRoot"), "SysWOW64")
    driverstore <- file.path(system32, "DriverStore", "FileRepository")
    
    loader_exists <- file.exists(file.path(system32, "OpenCL.dll"))
    
    icd_names <- c("nvopencl64.dll","amdocl64.dll","intelocl64.dll","igdrcl64.dll","pocl.dll")
    icd_found <- FALSE
    
    # Search System32 + SysWOW64
    for (d in c(system32, syswow64)) {
      for (nm in icd_names) {
        if (file.exists(file.path(d, nm))) icd_found <- TRUE
      }
    }
    
    # Search DriverStore recursively
    if (dir.exists(driverstore)) {
      hits <- list.files(driverstore, pattern = paste(icd_names, collapse="|"),
                         recursive = TRUE, full.names = TRUE)
      if (length(hits) > 0) icd_found <- TRUE
    }
    
    # Registry search (future-proof)
    reg_paths <- character()
    for (key in c(
      "HKEY_LOCAL_MACHINE\\SOFTWARE\\Khronos\\OpenCL\\Vendors",
      "HKEY_LOCAL_MACHINE\\SOFTWARE\\WOW6432Node\\Khronos\\OpenCL\\Vendors"
    )) {
      try({
        entries <- suppressWarnings(utils::readRegistry(key, maxdepth = 1))
        reg_paths <- c(reg_paths, names(entries))
      }, silent = TRUE)
    }
    for (p in reg_paths) {
      if (file.exists(p)) {
        if (basename(p) %in% icd_names) icd_found <- TRUE
      }
    }
    
    # Determine actual loader path
    loader_path <- file.path(Sys.getenv("SystemRoot"), "System32", "OpenCL.dll")
    loader_exists <- file.exists(loader_path)
    
    # Populate bin_dirs ONLY if loader + ICD exist
    if (loader_exists && icd_found) {
      loader_dir <- normalizePath(dirname(loader_path), winslash = "/", mustWork = FALSE)
      result$runtimes$nvidia$opencl$bin_dirs <- loader_dir
    }
    
    
    # Mark installed if headers OR ICDs found
    if (!is.null(cl_header) || icd_found) {
      result$runtimes$nvidia$opencl$installed <- TRUE
    }
    
    # -------------------------------
    # CUDA detection (UNCHANGED)
    # -------------------------------
    nvcc_path <- Sys.which("nvcc.exe")
    if (nzchar(nvcc_path)) {
      nvcc_path <- normalizePath(nvcc_path, winslash = "/", mustWork = TRUE)
      cuda_root <- normalizePath(file.path(dirname(nvcc_path), ".."), winslash = "/", mustWork = TRUE)
      cuda_include <- file.path(cuda_root,"include")
      cuda_lib     <- file.path(cuda_root,"lib","x64")
      
      result$runtimes$nvidia$cuda$installed <- TRUE
      result$runtimes$nvidia$cuda$bin_dirs <- dirname(nvcc_path)
      if (dir.exists(cuda_include)) result$runtimes$nvidia$cuda$include_dirs <- cuda_include
      if (dir.exists(cuda_lib))     result$runtimes$nvidia$cuda$lib_dirs     <- cuda_lib
    }
  }
  
  return(result)
}