#' Detect GPU Drivers and Provide Installation Instructions
#'
#' This function analyzes the output of `detect_environment_and_gpus()` and
#' determines whether the appropriate GPU drivers are installed for the detected
#' hardware and environment. It does not install drivers automatically. Instead,
#' it returns a structured list describing which drivers are present, which are
#' missing, and what actions the user should take.
#'
#' @param info A list returned by `detect_environment_and_gpus()`. The list must
#'   contain the following elements:
#'   \describe{
#'     \item{environment}{One of "windows", "msys2", "linux", "wsl", or "unknown".}
#'     \item{nvidia}{A list with elements `present` (logical) and `names` (character).}
#'     \item{amd}{A list with elements `present` (logical) and `names` (character).}
#'     \item{intel}{A list with elements `present` (logical) and `names` (character).}
#'   }
#'
#' @return A structured list describing driver status for each GPU vendor. Each
#'   vendor entry contains:
#'   \describe{
#'     \item{installed}{Logical indicating whether the driver appears installed.}
#'     \item{issues}{Character vector of detected problems.}
#'     \item{actions}{Character vector of recommended user actions.}
#'   }
#'
#' @examples
#' \dontrun{
#'   info <- detect_environment_and_gpus()
#'   detect_or_install_gpu_drivers(info)
#' }
#'
#' @export
detect_or_install_gpu_drivers <- function(info) {
  env <- info$environment
  
  result <- list(
    environment = env,
    drivers = list(
      nvidia = list(installed = FALSE, issues = character(), actions = character()),
      amd    = list(installed = FALSE, issues = character(), actions = character()),
      intel  = list(installed = FALSE, issues = character(), actions = character())
    )
  )
  
  # ------------------------------------------------------------------
  # NVIDIA DRIVER CHECK
  # ------------------------------------------------------------------
  if (isTRUE(info$nvidia$present)) {
    
    if (env %in% c("windows", "msys2")) {
      
      if (nzchar(Sys.which("nvidia-smi"))) {
        result$drivers$nvidia$installed <- TRUE
      } else {
        result$drivers$nvidia$issues <- c("NVIDIA driver not detected on Windows")
        result$drivers$nvidia$actions <- c(
          "Install NVIDIA Studio Driver from https://www.nvidia.com/Download"
        )
      }
      
    } else if (env == "linux") {
      
      if (nzchar(Sys.which("nvidia-smi"))) {
        result$drivers$nvidia$installed <- TRUE
      } else {
        result$drivers$nvidia$issues <- c("NVIDIA driver not detected on Linux")
        result$drivers$nvidia$actions <- c(
          "Install NVIDIA driver: sudo ubuntu-drivers autoinstall",
          "Or manually: sudo apt install nvidia-driver-550"
        )
      }
      
    } else if (env == "wsl") {
      
      if (!file.exists("/dev/dxg")) {
        result$drivers$nvidia$issues <- c(
          "WSL GPU virtualization not available (/dev/dxg missing)"
        )
        result$drivers$nvidia$actions <- c(
          "Your GPU cannot be exposed to WSL. Use Windows-native GPU compute instead."
        )
      } else if (!nzchar(Sys.which("nvidia-smi"))) {
        result$drivers$nvidia$issues <- c(
          "NVIDIA driver not detected inside WSL"
        )
        result$drivers$nvidia$actions <- c(
          "Install NVIDIA Studio Driver in Windows. WSL does not use Linux NVIDIA drivers."
        )
      } else {
        result$drivers$nvidia$installed <- TRUE
      }
    }
  }
  
  # ------------------------------------------------------------------
  # AMD DRIVER CHECK
  # ------------------------------------------------------------------
  if (isTRUE(info$amd$present)) {
    
    if (env %in% c("windows", "msys2")) {
      
      result$drivers$amd$installed <- TRUE
      
    } else if (env == "linux") {
      
      icd_files <- character()
      if (dir.exists("/etc/OpenCL/vendors")) {
        icd_files <- list.files("/etc/OpenCL/vendors", full.names = TRUE)
      }
      
      if (any(grepl("amdocl", icd_files))) {
        result$drivers$amd$installed <- TRUE
      } else {
        result$drivers$amd$issues <- c("AMD OpenCL ICD not detected")
        result$drivers$amd$actions <- c(
          "Install AMD OpenCL or ROCm packages (distribution specific)"
        )
      }
      
    } else if (env == "wsl") {
      
      result$drivers$amd$issues <- c(
        "AMD GPUs cannot be exposed to WSL. WSL GPU compute requires NVIDIA."
      )
      result$drivers$amd$actions <- c(
        "Use Windows-native GPU compute instead."
      )
    }
  }
  
  # ------------------------------------------------------------------
  # INTEL DRIVER CHECK
  # ------------------------------------------------------------------
  if (isTRUE(info$intel$present)) {
    
    if (env %in% c("windows", "msys2")) {
      
      result$drivers$intel$installed <- TRUE
      
    } else if (env == "linux") {
      
      icd_files <- character()
      if (dir.exists("/etc/OpenCL/vendors")) {
        icd_files <- list.files("/etc/OpenCL/vendors", full.names = TRUE)
      }
      
      if (any(grepl("intel", icd_files, ignore.case = TRUE))) {
        result$drivers$intel$installed <- TRUE
      } else {
        result$drivers$intel$issues <- c("Intel OpenCL ICD not detected")
        result$drivers$intel$actions <- c(
          "Install Intel OpenCL runtime (package name varies by distribution)"
        )
      }
      
    } else if (env == "wsl") {
      
      result$drivers$intel$issues <- c(
        "Intel GPUs cannot be exposed to WSL. WSL GPU compute requires NVIDIA."
      )
      result$drivers$intel$actions <- c(
        "Use Windows-native GPU compute instead."
      )
    }
  }
  
  return(result)
}