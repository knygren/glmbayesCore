#' Detect host environment and available GPUs
#'
#' \code{detect_environment_and_gpus()} inspects the host system to determine
#' whether R is running on Windows, Linux, WSL, or MSYS2, and then performs
#' lightweight, vendor-agnostic checks for NVIDIA, AMD, and Intel GPUs.
#'
#' Detection is intentionally minimal and relies only on tools that are expected
#' to be present when a GPU stack is installed:
#' \itemize{
#'   \item On Windows/MSYS2: \code{where.exe}, \code{nvidia-smi}, and
#'         PowerShell WMI (\code{Get-CimInstance Win32_VideoController}).
#'   \item On Linux/WSL: \code{nvidia-smi} (if NVIDIA drivers are installed)
#'         and \code{lspci} (for AMD/Intel detection).
#' }
#'
#' This function does not attempt to install drivers, CUDA, or OpenCL, and does
#' not assume any optional tools such as \code{clinfo} or \code{rocminfo}.
#' It is intended as a basic diagnostic to inform subsequent setup steps
#' (e.g., OpenCL or CUDA configuration) rather than a full configuration
#' manager.
#'
#' @details
#' Environment detection follows a simple heuristic:
#' \itemize{
#'   \item If \code{uname -s} contains \code{"MINGW"}, the environment is
#'         reported as \code{"msys2"}.
#'   \item Else if \code{Sys.info()[["sysname"]]} is \code{"Windows"}, the
#'         environment is reported as \code{"windows"}.
#'   \item Else if \code{/proc/version} exists and contains \code{"Microsoft"},
#'         the environment is reported as \code{"wsl"}.
#'   \item Otherwise, the environment is reported as \code{"linux"}.
#' }
#'
#' GPU detection logic:
#' \describe{
#'   \item{NVIDIA}{
#'     \itemize{
#'       \item On Windows/MSYS2, \code{where.exe nvidia-smi} is used to test
#'             for the presence of \code{nvidia-smi.exe}. If found, GPU names
#'             are obtained via
#'             \code{nvidia-smi --query-gpu=name --format=csv,noheader}.
#'       \item On Linux/WSL, \code{command -v nvidia-smi} is used to test
#'             for the presence of \code{nvidia-smi}. If found, GPU names
#'             are obtained via the same query.
#'     }
#'   }
#'
#'   \item{AMD and Intel on Windows/MSYS2}{
#'     PowerShell WMI is used:
#'     \code{Get-CimInstance Win32_VideoController | Select-Object -ExpandProperty Name}.
#'     GPU names containing \code{"AMD"} or \code{"Radeon"} are classified as
#'     AMD; names containing \code{"Intel"} are classified as Intel.
#'   }
#'
#'   \item{AMD and Intel on Linux/WSL}{
#'     If \code{lspci} is available, its output is scanned:
#'     \itemize{
#'       \item Lines containing \code{"AMD"} or \code{"ATI"} are classified
#'             as AMD GPUs.
#'       \item Lines containing \code{"Intel"} and \code{"VGA"} are classified
#'             as Intel GPUs.
#'     }
#'   }
#' }
#'
#' This function performs best-effort detection. For example, on systems without
#' \code{lspci}, AMD/Intel detection on Linux may return \code{present = FALSE}
#' even if a GPU is physically present. Similarly, if \code{nvidia-smi} is not
#' on \code{PATH}, NVIDIA GPUs will not be reported even if drivers are
#' installed.
#'
#' @return
#' A named list with the following components:
#' \describe{
#'   \item{\code{environment}}{Character scalar indicating the detected host
#'     environment. One of \code{"windows"}, \code{"msys2"}, \code{"wsl"},
#'     \code{"linux"}, or \code{"unknown"}.}
#'
#'   \item{\code{nvidia}}{A list with elements:
#'     \describe{
#'       \item{\code{present}}{Logical; \code{TRUE} if NVIDIA GPUs were
#'         detected via \code{nvidia-smi}, \code{FALSE} otherwise.}
#'       \item{\code{names}}{Character vector of GPU names as reported by
#'         \code{nvidia-smi}, or \code{NULL} if not available.}
#'     }
#'   }
#'
#'   \item{\code{amd}}{A list with elements:
#'     \describe{
#'       \item{\code{present}}{Logical; \code{TRUE} if AMD GPUs were detected
#'         via WMI (Windows/MSYS2) or \code{lspci} (Linux/WSL).}
#'       \item{\code{names}}{Character vector of matching lines from WMI or
#'         \code{lspci}, or \code{NULL}/empty if none were found.}
#'     }
#'   }
#'
#'   \item{\code{intel}}{A list with elements:
#'     \describe{
#'       \item{\code{present}}{Logical; \code{TRUE} if Intel GPUs were detected
#'         via WMI (Windows/MSYS2) or \code{lspci} (Linux/WSL).}
#'       \item{\code{names}}{Character vector of matching lines from WMI or
#'         \code{lspci}, or \code{NULL}/empty if none were found.}
#'     }
#'   }
#' }
#'
#' @examples
#' \dontrun{
#'   info <- detect_environment_and_gpus()
#'   print(info$environment)
#'   print(info$nvidia)
#'   print(info$amd)
#'   print(info$intel)
#' }
#'
#' @export
detect_environment_and_gpus <- function() {
  
  # -------------------------------
  # 1. Detect environment
  # -------------------------------
  sysname <- Sys.info()[["sysname"]]
  
  uname_s <- try(system("uname -s", intern = TRUE), silent = TRUE)
  if (inherits(uname_s, "try-error")) {
    uname_s <- ""
  }
  
  if (grepl("MINGW", uname_s)) {
    env <- "msys2"
  } else if (identical(sysname, "Windows")) {
    env <- "windows"
  } else {
    if (file.exists("/proc/version")) {
      v <- readLines("/proc/version", warn = FALSE)
      if (any(grepl("Microsoft", v, ignore.case = TRUE))) {
        env <- "wsl"
      } else {
        env <- "linux"
      }
    } else {
      env <- "unknown"
    }
  }
  
  # -------------------------------
  # 2. Detect NVIDIA GPU
  # -------------------------------
  has_nvidia   <- FALSE
  nvidia_names <- NULL
  
  if (env %in% c("windows", "msys2")) {
    
    nvidia_path <- suppressWarnings(
      try(system("where.exe nvidia-smi", intern = TRUE), silent = TRUE)
    )
    
    has_nvidia <- !inherits(nvidia_path, "try-error") &&
      length(nvidia_path) > 0L &&
      any(nzchar(nvidia_path))
    
    if (has_nvidia) {
      nvidia_names <- try(
        system("nvidia-smi --query-gpu=name --format=csv,noheader", intern = TRUE),
        silent = TRUE
      )
      if (inherits(nvidia_names, "try-error")) {
        nvidia_names <- NULL
      }
    }
    
  } else {
    
    has_nvidia <- (system("command -v nvidia-smi >/dev/null 2>&1",
                          intern = FALSE) == 0L)
    
    if (has_nvidia) {
      nvidia_names <- try(
        system("env LD_LIBRARY_PATH=/usr/lib/wsl/lib:$LD_LIBRARY_PATH /usr/bin/nvidia-smi --query-gpu=name --format=csv,noheader", intern = TRUE),
        silent = TRUE
      )
      if (inherits(nvidia_names, "try-error")) {
        nvidia_names <- NULL
      }
    }
  }
  
  # -------------------------------
  # 3. Detect AMD + Intel on Windows/MSYS2
  # -------------------------------
  amd_names   <- NULL
  intel_names <- NULL
  has_amd     <- FALSE
  has_intel   <- FALSE
  
  if (env %in% c("windows", "msys2")) {
    
    gpu_list <- suppressWarnings(
      try(
        system(
          'powershell -Command "Get-CimInstance Win32_VideoController | Select-Object -ExpandProperty Name"',
          intern = TRUE
        ),
        silent = TRUE
      )
    )
    
    if (!inherits(gpu_list, "try-error") && length(gpu_list) > 0L) {
      
      amd_names   <- grep("AMD|Radeon", gpu_list, value = TRUE, ignore.case = TRUE)
      intel_names <- grep("Intel", gpu_list, value = TRUE, ignore.case = TRUE)
      
      has_amd   <- length(amd_names)   > 0L
      has_intel <- length(intel_names) > 0L
    }
  }
  
  # -------------------------------
  # 4. Detect AMD + Intel on Linux (NOT WSL)
  # -------------------------------
  if (env == "linux") {
    
    if (system("command -v lspci >/dev/null 2>&1", intern = FALSE) == 0L) {
      
      pci <- try(system("lspci", intern = TRUE), silent = TRUE)
      
      if (!inherits(pci, "try-error") && length(pci) > 0L) {
        
        controller_lines <- grep("VGA compatible controller|3D controller",
                                 pci, value = TRUE, ignore.case = TRUE)
        
        amd_names_linux <- grep("(AMD|ATI)", controller_lines,
                                value = TRUE, ignore.case = TRUE)
        
        intel_names_linux <- grep("Intel", controller_lines,
                                  value = TRUE, ignore.case = TRUE)
        
        if (length(amd_names_linux) > 0L) {
          amd_names <- amd_names_linux
          has_amd <- TRUE
        }
        
        if (length(intel_names_linux) > 0L) {
          intel_names <- intel_names_linux
          has_intel <- TRUE
        }
      }
    }
  }
  
  # -------------------------------
  # 5. WSL: force AMD/Intel to FALSE
  # -------------------------------
  if (env == "wsl") {
    has_amd <- FALSE
    has_intel <- FALSE
    amd_names <- NULL
    intel_names <- NULL
  }
  
  # -------------------------------
  # Return structured result
  # -------------------------------
  list(
    environment = env,
    nvidia = list(
      present = isTRUE(has_nvidia),
      names   = nvidia_names
    ),
    amd = list(
      present = isTRUE(has_amd),
      names   = amd_names
    ),
    intel = list(
      present = isTRUE(has_intel),
      names   = intel_names
    )
  )
}