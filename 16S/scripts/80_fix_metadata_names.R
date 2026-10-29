#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

source("config/config.R")

stopifnot(file.exists(METADATA_CSV))
meta <- read.csv(
  METADATA_CSV,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

meta <- meta[, names(meta) != "", drop = FALSE]
stopifnot(METADATA_SAMPLE_COL == "Sample_clean")
sample_ids <- rownames(meta)
stopifnot(!anyNA(sample_ids), all(nzchar(sample_ids)), !anyDuplicated(sample_ids))
stopifnot(METADATA_SAMPLE_COL %in% names(meta))
existing_ids <- as.character(meta[[METADATA_SAMPLE_COL]])

# Apply only the explicit identifier aliases in config; retain tissue metadata.
applied_patches <- data.frame(from = character(), to = character())
if (exists("METADATA_PATCHES") && nrow(METADATA_PATCHES) > 0L) {
  stopifnot(is.data.frame(METADATA_PATCHES),
            all(c("from", "to") %in% names(METADATA_PATCHES)))
  patches <- METADATA_PATCHES[, c("from", "to"), drop = FALSE]
  patches[] <- lapply(patches, as.character)
  stopifnot(!anyNA(patches), all(nzchar(patches$from)), all(nzchar(patches$to)),
            !anyDuplicated(patches$from), !anyDuplicated(patches$to))
  used <- patches$from %in% c(sample_ids, existing_ids)
  applied_patches <- patches[used, , drop = FALSE]
  apply_aliases <- function(ids) {
    matched <- match(ids, patches$from)
    ids[!is.na(matched)] <- patches$to[matched[!is.na(matched)]]
    ids
  }
  sample_ids <- apply_aliases(sample_ids)
  existing_ids <- apply_aliases(existing_ids)
}
has_id <- !is.na(existing_ids) & nzchar(trimws(existing_ids))
stopifnot(all(existing_ids[has_id] == sample_ids[has_id]), !anyDuplicated(sample_ids))
rownames(meta) <- sample_ids
meta[[METADATA_SAMPLE_COL]] <- sample_ids

out_csv <- file.path("outputs", "metadata_fixed.csv")
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
fwrite(meta, out_csv)
write.csv(applied_patches, file.path("outputs", "metadata_name_corrections_applied.csv"), row.names = FALSE)
cat("Saved: ", out_csv, "\n", sep = "")
