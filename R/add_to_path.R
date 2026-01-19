#' @title Add Directories to PATH or LD_LIBRARY_PATH
#'
#' @description
#' These helper functions allow you to add missing directories to the PATH
#' or library search environment variables in a permanent way, minimizing
#' manual editing.
#'
#' @details
#' - On **Windows**, updates the user-level PATH via PowerShell.
#' - On **Linux/WSL**, appends export lines to ~/.bashrc for PATH or LD_LIBRARY_PATH.
#'
#' @param dirs Character vector of directories to add.
#'
#' @return No return value; called for side effects.
#'
#' @examples
#' \dontrun{
#' add_to_path_windows(c("C:/OpenCL/bin"))
#' add_to_path_linux(c("/usr/local/cuda/bin"))
#' add_to_libpath_linux(c("/usr/local/cuda/lib64"))
#' }
#'
#' @seealso [Sys.getenv], [Sys.setenv]
#'
#' @name add_to_path
#' @rdname add_to_path
NULL

#' @export
#' @rdname add_to_path
add_to_path_windows <- function(dirs) {
  for (d in dirs) {
    cmd <- sprintf(
      '[System.Environment]::SetEnvironmentVariable("PATH", "%s;%s", "User")',
      d, Sys.getenv("PATH")
    )
    system2("powershell", c("-Command", cmd))
  }
  message("[ACTION] Permanently added missing PATH dirs to user environment.")
  message("         Restart your shell/session to see changes.")
}

#' @export
#' @rdname add_to_path
add_to_path_linux <- function(dirs) {
  bashrc <- file.path(Sys.getenv("HOME"), ".bashrc")
  for (d in dirs) {
    line <- sprintf('export PATH="%s:$PATH"', d)
    write(line, bashrc, append = TRUE)
  }
  message("[ACTION] Permanently added missing PATH dirs to ~/.bashrc.")
  message("         Restart your shell or run `source ~/.bashrc` to apply changes.")
}

#' @export
#' @rdname add_to_path
add_to_libpath_linux <- function(dirs) {
  bashrc <- file.path(Sys.getenv("HOME"), ".bashrc")
  for (d in dirs) {
    line <- sprintf('export LD_LIBRARY_PATH="%s:$LD_LIBRARY_PATH"', d)
    write(line, bashrc, append = TRUE)
  }
  message("[ACTION] Permanently added missing library dirs to ~/.bashrc (LD_LIBRARY_PATH).")
  message("         Restart your shell or run `source ~/.bashrc` to apply changes.")
}