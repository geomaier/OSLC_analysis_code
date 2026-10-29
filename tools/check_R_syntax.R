#!/usr/bin/env Rscript
# Parse source only. Do not source or execute analysis scripts.
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else "."
root <- normalizePath(root, winslash = "/", mustWork = TRUE)
files <- list.files(root, pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE,
                    all.files = FALSE)
if (!length(files)) stop("No R files found under: ", root)
failed <- character()
for (path in sort(files)) {
  error <- tryCatch({parse(file = path, keep.source = FALSE); NULL},
                    error = function(e) conditionMessage(e))
  relative <- substring(path, nchar(root) + 2L)
  if (is.null(error)) {
    cat("PASS ", relative, "\n", sep = "")
  } else {
    cat("FAIL ", relative, ": ", error, "\n", sep = "")
    failed <- c(failed, relative)
  }
}
cat(length(files), "R files parsed;", length(failed), "failure(s).\n")
if (length(failed)) quit(save = "no", status = 1L)
