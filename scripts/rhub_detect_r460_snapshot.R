# Used by .github/workflows/rhub.yaml: "Detect R 4.6.0 low-svn snapshot" on all OSes.
# Windows: `Rscript scripts/rhub_detect_r460_snapshot.R` (avoids long `Rscript -e` in Git Bash).
# Excluded from R CMD build (see .Rbuildignore); version-controlled via .gitignore exception.
# Writes two lines to stdout: (1) true/false, (2) short note (env MATRIX_LABEL).
lab <- Sys.getenv("MATRIX_LABEL", "")
rv <- getRversion()
v <- R.version
raw <- if (!is.null(v[["svn.rev"]])) v[["svn.rev"]] else v[["svn rev"]]
svn <- suppressWarnings(as.integer(as.character(raw)))
is460 <- identical(as.character(rv), "4.6.0")
relax <- !identical(lab, "rchk") &&
  isTRUE(is460) && length(svn) == 1L && !is.na(svn) && svn < 89746L
line1 <- if (relax) "true" else "false"
line2 <- if (relax) {
  sprintf("R 4.6.0 svn %d (<89746): DESCRIPTION relax + Ensure Rcpp", svn)
} else {
  "none"
}
writeLines(c(line1, line2), stdout())
invisible()
