## Compare glmbayes vs glmbayesCore (R / src / inst/cl) and write assessment README.
##
## Usage (from glmbayesCore package root):
##   source("data-raw/make_compare_glmbayes_glmbayesCore.R")
##
## Optional:
##   Sys.setenv(GLMBAYES_COMPARE_ROOT = "C:/path/to/glmbayes")
##
## Outputs:
##   data-raw/COMPARE_glmbayes_vs_glmbayesCore.md
##   data-raw/compare_glmbayes_glmbayesCore.json

root <- if (basename(getwd()) == "data-raw") {
  normalizePath("..")
} else {
  normalizePath(".")
}

py <- file.path(root, "data-raw", "make_compare_glmbayes_glmbayesCore.py")
if (!file.exists(py)) {
  stop("Missing regenerator script: ", py)
}

# Prefer `python`, then Windows `py -3`
runners <- list(
  c("python", py),
  c("py", "-3", py)
)

ok <- FALSE
last_status <- NULL
for (cmd in runners) {
  status <- tryCatch(
    system2(cmd[[1]], args = cmd[-1], stdout = TRUE, stderr = TRUE),
    error = function(e) e
  )
  if (inherits(status, "error")) {
    last_status <- conditionMessage(status)
    next
  }
  code <- attr(status, "status")
  if (is.null(code) || identical(as.integer(code), 0L)) {
    writeLines(as.character(status))
    ok <- TRUE
    break
  }
  last_status <- paste(as.character(status), collapse = "\n")
}

if (!ok) {
  stop(
    "Failed to run Python regenerator. Install Python 3 and ensure ",
    "`python` or `py -3` is on PATH.\n",
    last_status
  )
}

message(
  "Wrote:\n  ",
  file.path(root, "data-raw", "COMPARE_glmbayes_vs_glmbayesCore.md"),
  "\n  ",
  file.path(root, "data-raw", "compare_glmbayes_glmbayesCore.json")
)
