#!/usr/bin/env Rscript
# Regression checks for reporting and metadata guards; no study data or packages.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[1] else ".", winslash = "/", mustWork = TRUE)

local({
# Focused source-expression checks; no study data or external R packages.
repo_root <- root
read_source <- function(path) parse(file.path(repo_root, path), keep.source = FALSE)
lhs_text <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = "")
get_assignment <- function(expressions, lhs) {
  keep <- vapply(expressions, function(x)
    is.call(x) && identical(x[[1L]], as.name("<-")) &&
      identical(lhs_text(x[[2L]]), lhs), logical(1))
  if (sum(keep) != 1L) stop("Expected one assignment to ", lhs)
  expressions[[which(keep)]]
}
get_function <- function(expressions, name) {
  eval(get_assignment(expressions, name)[[3L]], new.env(parent = globalenv()))
}
body_expressions <- function(fun) as.list(body(fun))[-1L]
checked <- character()
check <- function(name, expr) {
  force(expr)
  checked <<- c(checked, name)
}
throws <- function(expr) inherits(tryCatch({force(expr); NULL}, error = identity), "error")

main_expr <- read_source("16S/scripts/160_deseq2_volcano_all_comparisons.R")
matched_expr <- read_source("16S/scripts/175_matched_analysis.R")
report_expr <- read_source("16S/R/reporting_corrections.R")

volcano_rules <- function(fun, object_name) {
  pipeline <- get_assignment(body_expressions(fun), object_name)[[3L]]
  filter_call <- pipeline[[2L]][[3L]]
  mutate_call <- pipeline[[3L]]
  stopifnot(identical(filter_call[[1L]], as.name("filter")),
            identical(mutate_call[[1L]], as.name("mutate")))
  list(predicates = as.list(filter_call)[-1L],
       ordinate = as.list(mutate_call)[["FDR_plot"]])
}
main_rules <- volcano_rules(get_function(main_expr, "run_deseq2_comparison"), "volcano_df")
matched_rules <- volcano_rules(get_function(matched_expr, "make_aligned_volcano"), "plot_df")
final_main_rules <- volcano_rules(get_function(report_expr, "replot_main_volcano"), "volcano_df")
final_matched_rules <- volcano_rules(get_function(report_expr, "make_aligned_volcano"), "plot_df")

check("Fresh volcano rules equal the established final reporting rules", {
  stopifnot(identical(main_rules, final_main_rules),
            identical(matched_rules, final_matched_rules))
})
check("Missing and invalid FDRs are omitted without replacing raw values", {
  fixture <- data.frame(FDR = c(0.01, 0.5, NA, NaN, Inf, -0.1, 1.1, 0, 0.1),
                        log2FoldChange = c(1, -1, 1, 1, 1, 1, 1, 0, Inf))
  before <- fixture
  for (rules in list(main_rules, matched_rules)) {
    include <- Reduce(`&`, lapply(rules$predicates, eval, envir = fixture))
    stopifnot(identical(which(include), c(1L, 2L, 8L)))
    displayed <- eval(rules$ordinate, envir = fixture[include, , drop = FALSE])
    stopifnot(identical(displayed, c(0.01, 0.5, .Machine$double.xmin)),
              all(is.finite(-log10(displayed))), identical(fixture, before))
  }
})

ancom_function <- get_function(matched_expr, "run_ancombc2_rank")
ancom_body <- body_expressions(ancom_function)
check("ANCOM unavailable and zero-information tests remain NA in reporting fields", {
  env <- new.env(parent = globalenv())
  env$tidy <- data.frame(
    Residual_df = c(5, 0, 5, 5, 5),
    Contrast_estimable = c(TRUE, TRUE, FALSE, TRUE, TRUE),
    lfc_natural_log = c(1, 0, 2, 3, Inf), se = c(1, 0, 1, 1, 1),
    W = c(1, 0, 2, 3, 1), p_value = c(.01, 1, 1, NaN, .1),
    q_value_BH = c(.02, 1, 1, NaN, .2), q_value_Holm = c(.04, 1, 1, NaN, .4),
    log2_abundance_difference = c(1, 0, 2, 3, Inf) / log(2)
  )
  before <- env$tidy
  for (field in c("Test_available", "Reporting_P", "Reporting_BH_FDR",
                  "Reporting_Holm_P", "Reporting_log2_difference")) {
    eval(get_assignment(ancom_body, paste0("tidy$", field)), env)
  }
  stopifnot(identical(env$tidy$Test_available, c(TRUE, FALSE, FALSE, FALSE, FALSE)),
            identical(env$tidy[, names(before)], before),
            identical(env$tidy$Reporting_P[1L], .01))
  reporting <- env$tidy[, c("Reporting_P", "Reporting_BH_FDR", "Reporting_Holm_P",
                           "Reporting_log2_difference")]
  stopifnot(all(is.na(reporting[-1L, , drop = FALSE])))
})
cross_function <- get_function(matched_expr, "make_cross_method_table")
final_cross_function <- get_function(report_expr, "make_cross_method_table")
check("Fresh cross-method table uses exactly the established reporting function", {
  stopifnot(identical(body(cross_function), body(final_cross_function)),
            identical(formals(cross_function), formals(final_cross_function)))
})
check("Unavailable tests and structural-zero-only rows have no direction agreement", {
  pipeline <- get_assignment(body_expressions(cross_function), "out")[[3L]]
  mutate_call <- pipeline[[2L]][[3L]]
  cases <- as.list(mutate_call)[["Same_effect_direction"]]
  unavailable_condition <- cases[[2L]][[2L]]
  fixture <- data.frame(
    DESeq2_FDR = c(.01, .01, NA, .01, .01),
    ANCOMBC2_test_available = c(TRUE, FALSE, TRUE, NA, TRUE),
    DESeq2_log2FoldChange = c(1, 1, 1, 1, Inf),
    ANCOMBC2_log2_difference = c(1, 0, 1, NA, 1)
  )
  stopifnot(identical(eval(unavailable_condition, fixture),
                      c(FALSE, TRUE, TRUE, TRUE, TRUE)))
})

metadata_expr <- read_source("16S/scripts/80_fix_metadata_names.R")
begin <- which(vapply(metadata_expr, function(x)
  is.call(x) && identical(x[[1L]], as.name("<-")) &&
    identical(lhs_text(x[[2L]]), "sample_ids"), logical(1)))
finish <- which(vapply(metadata_expr, function(x)
  is.call(x) && identical(x[[1L]], as.name("<-")) &&
    identical(lhs_text(x[[2L]]), "meta[[METADATA_SAMPLE_COL]]"), logical(1)))
stopifnot(length(begin) == 1L, length(finish) == 1L)
metadata_segment <- metadata_expr[seq.int(begin, finish)]
old_id <- "B3-HP_TRVS_8979"
new_id <- "B3-SSL_TRVS_8979"
alias_fixture <- data.frame(from = old_id, to = new_id)
run_metadata <- function(ids, supplied_ids = ids, aliases = alias_fixture) {
  env <- new.env(parent = globalenv())
  env$meta <- data.frame(Sample_clean = supplied_ids,
                        Sample_Type = rep("SSL", length(ids)),
                        row.names = ids, stringsAsFactors = FALSE)
  env$METADATA_SAMPLE_COL <- "Sample_clean"
  env$METADATA_PATCHES <- aliases
  for (expr in metadata_segment) eval(expr, env)
  list(meta = env$meta, patches = env$applied_patches)
}
check("Configured aliases harmonise identifiers and preserve tissue metadata", {
  answer <- run_metadata(c(old_id, "control"))
  stopifnot(identical(rownames(answer$meta), c(new_id, "control")),
            identical(answer$meta$Sample_clean, c(new_id, "control")),
            identical(answer$meta$Sample_Type, c("SSL", "SSL")),
            identical(answer$patches, alias_fixture))
  mixed <- run_metadata(c(old_id, "control"), c(new_id, "control"))
  stopifnot(identical(mixed$meta, answer$meta))
})
check("Alias application is idempotent for already corrected input", {
  answer <- run_metadata(c(new_id, "control"))
  stopifnot(nrow(answer$patches) == 0L,
            identical(rownames(answer$meta), c(new_id, "control")))
})
check("Aliases reject collisions, unrelated mismatches and ambiguous mappings", {
  stopifnot(throws(run_metadata(c(old_id, new_id))),
            throws(run_metadata(c("one", "two"), c("different", "two"))),
            throws(run_metadata(c(old_id, "control"),
                                aliases = rbind(alias_fixture, alias_fixture))))
})

oslc_16s_review_results <- data.frame(Check = checked, Passed = TRUE)
print(oslc_16s_review_results, row.names = FALSE)
cat("Passed", nrow(oslc_16s_review_results), "16S source-expression checks.\n")

})

local({
  script_path <- file.path(root, "Organoid_Bulk_RNAseq/reviewed/03_MUC5AC_qPCR.R")
program <- parse(file = script_path)
hits <- which(vapply(program, function(node) {
  is.call(node) && length(node) >= 3L &&
    identical(node[[1]], as.name("<-")) &&
    identical(node[[2]], as.name("replicate_qc"))
}, logical(1)))
stopifnot(length(hits) == 1L)
assignment <- program[[hits]]
mutation <- assignment[[3]][[3]]
stopifnot(is.call(mutation), identical(mutation[[1]], as.name("mutate")))
mutation_arguments <- as.list(mutation)
removed_expression <- mutation_arguments[["Automatically_excluded_as_outlier"]]
used_expression <- mutation_arguments[["Used_in_primary_sample_mean"]]
stopifnot(!is.null(removed_expression), !is.null(used_expression))
fixture <- list2env(list(
  Replicate_number = c(1L, 2L, 3L, 1L, 2L, 3L, 1L, 2L, 3L),
  Excluded_replicate_number = c(NA_integer_, NA_integer_, NA_integer_,
                               2L, 2L, 2L, NA_integer_, NA_integer_, NA_integer_),
  Finite_measurement = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
  Primary_analysis_included = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)
), parent = baseenv())
fixture$Automatically_excluded_as_outlier <- eval(removed_expression, fixture)
used <- eval(used_expression, fixture)
stopifnot(
  identical(fixture$Automatically_excluded_as_outlier,
            c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)),
  identical(used,
            c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE)),
  !anyNA(fixture$Automatically_excluded_as_outlier),
  !anyNA(used)
)
cat('PASS: actual qPCR source expressions correctly flag no-outlier, excluded-well, missing-well and excluded-sample cases.\n')

})
