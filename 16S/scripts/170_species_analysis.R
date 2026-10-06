#!/usr/bin/env Rscript
local({
patch_ancombc_bias_em_if_needed <- function() {
  namespace <- asNamespace("ANCOMBC")

  if (!exists(".bias_em", envir = namespace, inherits = FALSE)) {
    stop(
      "The installed ANCOMBC namespace does not contain .bias_em. ",
      "Do not run the analysis with an unknown implementation.",
      call. = FALSE
    )
  }

  installed_version <- packageVersion("ANCOMBC")
  bias_em_function <- get(
    ".bias_em",
    envir = namespace,
    inherits = FALSE
  )

  source_before <- paste(
    deparse(body(bias_em_function)),
    collapse = "\n"
  )

  if (grepl("neither_na", source_before, fixed = TRUE)) {
    return(
      data.frame(
        ANCOMBC_version = as.character(installed_version),
        Bias_EM_patch_required = FALSE,
        Bias_EM_patch_applied_in_memory = FALSE,
        Patch_status = "Installed ANCOMBC already contains the corrected joint NA mask",
        stringsAsFactors = FALSE
      )
    )
  }

  if (installed_version != package_version("2.4.0")) {
    stop(
      "ANCOMBC ",
      as.character(installed_version),
      " does not contain the corrected joint NA mask. ",
      "This script patches only the audited ANCOMBC 2.4.0 implementation. ",
      "Install a current official ANCOMBC release before running.",
      call. = FALSE
    )
  }

  body_parts <- as.list(body(bias_em_function))

  normalise_expression <- function(expression_object) {
    output <- paste(deparse(expression_object), collapse = "")
    output <- gsub("\\s+", "", output)
    gsub("<-", "=", output, fixed = TRUE)
  }

  if (length(body_parts) < 5L) {
    stop("Unexpected ANCOMBC 2.4.0 .bias_em structure.", call. = FALSE)
  }

  observed_opening <- vapply(
    body_parts[2:4],
    normalise_expression,
    character(1)
  )

  expected_opening <- c(
    "beta=beta[!is.na(beta)]",
    "nu0=var_hat",
    "nu0=nu0[!is.na(nu0)]"
  )

  if (!identical(observed_opening, expected_opening)) {
    stop(
      "The ANCOMBC 2.4.0 .bias_em function does not match the audited ",
      "implementation. The namespace was not modified.",
      call. = FALSE
    )
  }

  corrected_body_parts <- c(
    body_parts[1],
    list(
      quote(nu0 <- var_hat),
      quote(neither_na <- !(is.na(beta) | is.na(nu0))),
      quote(beta <- beta[neither_na]),
      quote(nu0 <- nu0[neither_na])
    ),
    body_parts[5:length(body_parts)]
  )

  body(bias_em_function) <- as.call(corrected_body_parts)

  assignInNamespace(
    x = ".bias_em",
    value = bias_em_function,
    ns = "ANCOMBC"
  )

  source_after <- paste(
    deparse(
      body(
        get(
          ".bias_em",
          envir = asNamespace("ANCOMBC"),
          inherits = FALSE
        )
      )
    ),
    collapse = "\n"
  )

  if (!grepl("neither_na", source_after, fixed = TRUE)) {
    stop("The in-memory ANCOMBC patch could not be verified.", call. = FALSE)
  }

  data.frame(
    ANCOMBC_version = as.character(installed_version),
    Bias_EM_patch_required = TRUE,
    Bias_EM_patch_applied_in_memory = TRUE,
    Patch_status = paste0(
      "Applied the corrected joint NA mask in memory; ",
      "installed package files were not modified"
    ),
    stringsAsFactors = FALSE
  )
}

clean_character <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}

sample_id_key <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\.", "-", x)
  toupper(x)
}

clean_plain_genus <- function(x) {
  x <- clean_character(x)
  x <- gsub("^g__", "", x, ignore.case = TRUE)
  x <- gsub("^Genus:", "", x, ignore.case = TRUE)
  x <- trimws(x)

  bad <- is.na(x) |
    toupper(x) %in% c(
      "NA", "N/A", "UNKNOWN", "UNCLASSIFIED",
      "UNCULTURED", "UNIDENTIFIED", "UNASSIGNED"
    )

  x[bad] <- NA_character_

  key <- tolower(gsub("[^a-z0-9]", "", x))
  x[!is.na(key) & key == "escherichiashigella"] <- "Escherichia-Shigella"
  x[!is.na(key) & key == "serratia"] <- "Serratia"
  x[!is.na(key) & key == "fusobacterium"] <- "Fusobacterium"

  x
}

canonical_site <- function(x) {
  x <- clean_character(x)
  if (length(x) == 0L) return(character(0))

  z <- toupper(x)
  z <- gsub("[^A-Z0-9]", "", z)

  out <- rep(NA_character_, length(z))

  out[z %in% c("CAEC", "CAECUM", "CECUM", "CECAL")] <- "CAEC"
  out[z %in% c("ASC", "ASCENDING", "ASCENDINGCOLON")] <- "ASC"
  out[z %in% c("TRVS", "TRANSVERSE", "TRANSVERSECOLON")] <- "TRVS"
  out[z %in% c("DESC", "DESCENDING", "DESCENDINGCOLON")] <- "DESC"
  out[z %in% c("SIG", "SIGMOID", "SIGMOIDCOLON")] <- "SIG"
  out[z %in% c("RECT", "RECTUM", "RECTAL")] <- "RECT"

  out
}

extract_patient_from_sample_id <- function(sample_id) {
  out <- sub(".*_([0-9]+)$", "\\1", sample_id)
  out[out == sample_id] <- NA_character_
  clean_character(out)
}

extract_site_from_sample_id <- function(sample_id) {
  site_codes <- c("CAEC", "ASC", "TRVS", "DESC", "SIG", "RECT")

  vapply(
    sample_id,
    function(id) {
      hits <- site_codes[
        vapply(
          site_codes,
          function(site) grepl(paste0("_", site, "_"), id, fixed = TRUE),
          logical(1)
        )
      ]

      if (length(hits) == 1L) hits else NA_character_
    },
    character(1)
  )
}

clean_taxonomy_value <- function(x) {
  x <- clean_character(x)

  bad <- is.na(x) |
    grepl(
      "^(NA|N/A|UNKNOWN|UNCLASSIFIED|UNCULTURED|UNIDENTIFIED|UNASSIGNED)$",
      toupper(x)
    )

  x[bad] <- NA_character_
  x
}

make_lineage_aware_genus <- function(tax_df) {
  tax_df <- as.data.frame(
    tax_df,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  required_ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

  for (rank in required_ranks) {
    if (!rank %in% colnames(tax_df)) tax_df[[rank]] <- NA_character_
    tax_df[[rank]] <- clean_taxonomy_value(tax_df[[rank]])
  }

  fallback_parent <- rep("Bacteria", nrow(tax_df))

  for (rank in c("Family", "Order", "Class", "Phylum", "Kingdom")) {
    use <- fallback_parent == "Bacteria" & !is.na(tax_df[[rank]])
    fallback_parent[use] <- tax_df[[rank]][use]
  }

  genus_base <- tax_df$Genus
  missing_genus <- is.na(genus_base)
  genus_base[missing_genus] <- paste0(
    "Unclassified_",
    gsub("[^A-Za-z0-9]+", "_", fallback_parent[missing_genus])
  )

  lineage_key <- apply(
    tax_df[, c("Phylum", "Class", "Order", "Family"), drop = FALSE],
    1,
    function(z) paste(ifelse(is.na(z), "NA", z), collapse = "|")
  )

  genus_label <- genus_base

  for (this_genus in unique(genus_base)) {
    idx <- which(genus_base == this_genus)
    distinct_lineages <- sort(unique(lineage_key[idx]))

    if (length(distinct_lineages) > 1L) {
      lineage_number <- match(lineage_key[idx], distinct_lineages)
      genus_label[idx] <- paste0(this_genus, "__lineage", lineage_number)
    }
  }

  mapping <- data.frame(
    Original_taxon = rownames(tax_df),
    Original_genus = tax_df$Genus,
    Fallback_parent = fallback_parent,
    Lineage_key = lineage_key,
    Genus_for_analysis = genus_label,
    stringsAsFactors = FALSE
  )

  list(genus = genus_label, mapping = mapping)
}
starts_with_genus <- function(species, genus) {
  if (is.na(species) || is.na(genus)) return(FALSE)
  species_lower <- tolower(species)
  genus_lower <- tolower(genus)
  if (!startsWith(species_lower, genus_lower)) return(FALSE)
  if (nchar(species_lower) == nchar(genus_lower)) return(TRUE)
  substr(species_lower, nchar(genus_lower) + 1L, nchar(genus_lower) + 1L) == " "
}

classify_species_label <- function(genus, species) {
  if (is.na(genus)) {
    return(c(
      Species_Key = NA_character_,
      Resolution_status = "Missing genus"
    ))
  }

  if (is.na(species)) {
    return(c(
      Species_Key = NA_character_,
      Resolution_status = "Missing species"
    ))
  }

  if (!grepl("^[A-Za-z][A-Za-z.-]+$", genus)) {
    return(c(
      Species_Key = NA_character_,
      Resolution_status = "Genus label is not a conventional single-word genus"
    ))
  }

  genus_canonical <- paste0(
    toupper(substr(genus, 1L, 1L)),
    tolower(substr(genus, 2L, nchar(genus)))
  )
  genus_lower <- tolower(genus)
  species_lower <- tolower(species)

  if (grepl("escherichia|shigella", paste(genus_lower, species_lower))) {
    return(c(
      Species_Key = NA_character_,
      Resolution_status = paste0(
        "Escherichia/Shigella lineage; not accepted for species-level 16S claims"
      )
    ))
  }

  if (grepl("[/;|]", species) || grepl("\\bor\\b", species_lower)) {
    return(c(
      Species_Key = NA_character_,
      Resolution_status = "Multiple species or ambiguous exact-match label"
    ))
  }

  ambiguous_pattern <- paste0(
    "\\b(uncultured|unclassified|unidentified|unknown|metagenome|",
    "environmental|microbiome|clone|bacterium|bacteria|candidatus)\\b|",
    "\\bcf\\.?\\b|\\baff\\.?\\b|",
    "(^|[[:space:]])sp\\.?$|(^|[[:space:]])spp\\.?$|",
    "[[:space:]]group$|[[:space:]]complex$"
  )

  if (grepl(ambiguous_pattern, species_lower, perl = TRUE)) {
    return(c(
      Species_Key = NA_character_,
      Resolution_status = "Ambiguous or non-species label"
    ))
  }

  if (identical(species_lower, genus_lower)) {
    return(c(
      Species_Key = NA_character_,
      Resolution_status = "Species field repeats genus only"
    ))
  }

  species_tokens_raw <- strsplit(
    gsub("\\s+", " ", trimws(species)),
    "\\s+"
  )[[1]]
  if (
    length(species_tokens_raw) >= 2L &&
    !starts_with_genus(species, genus)
  ) {
    return(c(
      Species_Key = NA_character_,
      Resolution_status = "Species binomial was inconsistent with the Genus column"
    ))
  }

  full_label <- if (starts_with_genus(species, genus)) {
    species
  } else {
    paste(genus_canonical, species)
  }

  full_label <- gsub("\\s+", " ", trimws(full_label))
  tokens <- strsplit(full_label, "\\s+")[[1]]

  if (length(tokens) < 2L) {
    return(c(
      Species_Key = NA_character_,
      Resolution_status = "Label is not a binomial species name"
    ))
  }

  tokens[1] <- genus_canonical
  tokens[2] <- paste0(
    tolower(substr(tokens[2], 1L, 1L)),
    substr(tokens[2], 2L, nchar(tokens[2]))
  )
  species_key <- paste(tokens[1:2], collapse = " ")

  if (!grepl("^[A-Z][A-Za-z.-]+ [a-z][A-Za-z0-9.-]+$", species_key)) {
    return(c(
      Species_Key = NA_character_,
      Resolution_status = "Binomial format did not pass the conservative audit"
    ))
  }

  c(
    Species_Key = species_key,
    Resolution_status = "Resolved provisional species label"
  )
}

fit_deseq_with_fallback <- function(dds, comparison) {
  fit_attempt_messages <- character(0)

  fit_types <- c("parametric", "local", "mean")
  for (fit_type in fit_types) {
    fitted_dds <- tryCatch(
      DESeq(
        dds,
        test = "Wald",
        fitType = fit_type,
        quiet = FALSE
      ),
      error = function(e) {
        fit_attempt_messages <<- c(
          fit_attempt_messages,
          paste0(
            "fitType = ", fit_type, " failed: ",
            conditionMessage(e)
          )
        )
        NULL
      }
    )

    if (!is.null(fitted_dds)) {
      return(list(
        dds = fitted_dds,
        strategy = paste0("DESeq fitType = ", fit_type),
        attempt_messages = fit_attempt_messages
      ))
    }
  }

  gene_wise_dds <- tryCatch(
    {
      candidate <- estimateDispersionsGeneEst(dds, quiet = FALSE)
      dispersions(candidate) <- S4Vectors::mcols(candidate)$dispGeneEst
      nbinomWaldTest(candidate, quiet = FALSE)
    },
    error = function(e) {
      fit_attempt_messages <<- c(
        fit_attempt_messages,
        paste0(
          "Gene-wise dispersion fallback failed: ",
          conditionMessage(e)
        )
      )
      NULL
    }
  )

  if (!is.null(gene_wise_dds)) {
    return(list(
      dds = gene_wise_dds,
      strategy = "Gene-wise dispersion estimates + nbinomWaldTest",
      attempt_messages = fit_attempt_messages
    ))
  }

  stop(
    "All DESeq2 dispersion/model-fitting strategies failed for ",
    comparison,
    ". Attempts: ",
    paste(fit_attempt_messages, collapse = " | "),
    call. = FALSE
  )
}

check_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
}

new_run <- function(work_dir, name, inputs, script_name) {
  if (!dir.exists(work_dir)) stop("Working directory does not exist: ", work_dir)
  absent <- inputs[!file.exists(inputs)]
  if (length(absent)) stop("Missing input: ", paste(absent, collapse = ", "))
  parent <- file.path(work_dir, name)
  dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(parent, paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid()))
  if (dir.exists(out) || !dir.create(out)) stop("Cannot create a new output folder: ", out)
  write.csv(data.frame(File = normalizePath(inputs), MD5 = unname(tools::md5sum(inputs))),
            file.path(out, "input_files.csv"), row.names = FALSE)
  script <- file.path(work_dir, script_name)
  if (file.exists(script)) {
    file.copy(script, file.path(out, script_name), overwrite = FALSE)
    write.csv(data.frame(File = script_name, MD5 = unname(tools::md5sum(script))),
              file.path(out, "script_checksum.csv"), row.names = FALSE)
  }
  out
}

csv <- function(x, path) {
  write.csv(as.data.frame(x, check.names = FALSE), path, row.names = FALSE, na = "")
}

matrix_csv <- function(x, path) {
  csv(data.frame(Feature = rownames(x), x, check.names = FALSE), path)
}

capture_fit <- function(expr) {
  warnings <- character()
  error <- NULL
  result <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, paste0(conditionMessage(w), " | ",
                                 paste(deparse(conditionCall(w)), collapse = " ")))
    invokeRestart("muffleWarning")
  }), error = function(e) {
    error <<- conditionMessage(e)
    NULL
  })
  list(value = result, warnings = unique(warnings), error = error)
}

write_messages <- function(x, path) {
  writeLines(if (length(x)) x else "None", path)
}

find_column <- function(meta, candidates) {
  found <- candidates[candidates %in% names(meta)]
  if (length(found)) found[1] else NULL
}

collapse_groups <- function(x) {
  x <- clean_character(x)
  x[x %in% c("SSL", "HP", "SSL_HP")] <- "SSL_HP"
  x[x %in% c("H_ADJ_SSL", "H_ADJ_HP", "H_ADJ_SSL_HP")] <- "H_ADJ_SSL_HP"
  x[x %in% c("H_COL", "H_COL_D", "H_COL_ND")] <- "H_COL"
  x[!x %in% c("SSL_HP", "H_ADJ_SSL_HP", "H_COL")] <- NA_character_
  x
}

load_input <- function(path, out) {
  phy <- readRDS(path)
  if (!inherits(phy, "phyloseq")) stop("Input is not a phyloseq object.")
  x <- as(phyloseq::otu_table(phy), "matrix")
  if (!phyloseq::taxa_are_rows(phy)) x <- t(x)
  meta <- data.frame(as(phyloseq::sample_data(phy), "data.frame"), check.names = FALSE)
  meta$SampleID <- rownames(meta)
  if (anyDuplicated(meta$SampleID) || anyDuplicated(colnames(x))) stop("Duplicate sample IDs.")
  meta <- meta[match(colnames(x), meta$SampleID), , drop = FALSE]
  if (anyNA(meta$SampleID)) stop("Sample metadata cannot be aligned to counts.")
  gc <- find_column(meta, c("Sample_Type", "SampleType", "sample_type", "Group", "New_Sample_Type"))
  if (is.null(gc)) stop("No tissue-group column found.")
  meta$Group <- collapse_groups(meta[[gc]])
  pc <- find_column(meta, c("Patient", "patient", "Patient_ID", "patient_id", "PatientID"))
  pm <- if (is.null(pc)) rep(NA_character_, nrow(meta)) else clean_character(meta[[pc]])
  pi <- extract_patient_from_sample_id(meta$SampleID)
  conflict <- !is.na(pm) & !is.na(pi) & pm != pi
  csv(meta[conflict, , drop = FALSE], file.path(out, "patient_id_conflicts.csv"))
  if (any(conflict)) stop("Patient metadata conflicts with sample IDs; inspect patient_id_conflicts.csv.")
  meta$Patient <- ifelse(is.na(pm), pi, pm)
  sc <- find_column(meta, c("Location_Colon", "LocationColon", "Colon_Location", "Site", "Location", "location"))
  sm <- if (is.null(sc)) rep(NA_character_, nrow(meta)) else canonical_site(meta[[sc]])
  si <- extract_site_from_sample_id(meta$SampleID)
  conflict <- !is.na(sm) & !is.na(si) & sm != si
  csv(meta[conflict, , drop = FALSE], file.path(out, "site_conflicts.csv"))
  if (any(conflict)) stop("Site metadata conflicts with sample IDs; inspect site_conflicts.csv.")
  meta$Site <- ifelse(is.na(sm), si, sm)
  keep <- !is.na(meta$Group)
  meta <- meta[keep, , drop = FALSE]
  x <- x[, keep, drop = FALSE]
  if (anyNA(meta$Patient)) stop("Patient IDs are missing.")
  if (any(meta$Patient == "7699")) stop("This workflow expects the final object excluding patient 7699.")
  if (any(!is.finite(x)) || any(x < 0) || any(abs(x - round(x)) > 1e-8)) {
    stop("Input must contain finite, non-negative integer counts.")
  }
  if (any(colSums(x) == 0)) stop("Input contains an empty biopsy.")
  expected <- c(SSL_HP = 38L, H_ADJ_SSL_HP = 12L, H_COL = 14L)
  observed <- table(factor(meta$Group, levels = names(expected)))
  if (!identical(as.integer(observed), unname(expected))) stop("Input cohort differs from 38 SSL_HP, 12 adjacent and 14 H_COL biopsies.")
  if (length(unique(meta$Patient)) != 31L) stop("Expected 31 patients in the final cohort.")
  tt <- as(phyloseq::tax_table(phy), "matrix")
  tt <- tt[match(rownames(x), rownames(tt)), , drop = FALSE]
  if (!identical(rownames(tt), rownames(x))) stop("Taxonomy and counts do not align.")
  rownames(meta) <- meta$SampleID
  csv(meta, file.path(out, "biopsy_metadata.csv"))
  matrix_csv(x, file.path(out, "asv_counts.csv"))
  csv(data.frame(ASV = rownames(tt), tt, check.names = FALSE), file.path(out, "taxonomy.csv"))
  list(counts = round(x), metadata = meta, taxonomy = tt)
}

exact_pairs <- function(meta, out) {
  m <- meta[meta$Group %in% c("SSL_HP", "H_ADJ_SSL_HP"), , drop = FALSE]
  csv(m[is.na(m$Site), , drop = FALSE], file.path(out, "missing_sites.csv"))
  m <- m[!is.na(m$Site), , drop = FALSE]
  m$PairID <- paste(m$Patient, m$Site, sep = "__")
  t <- table(m$PairID, factor(m$Group, levels = c("H_ADJ_SSL_HP", "SSL_HP")))
  audit <- data.frame(PairID = rownames(t), H_ADJ_SSL_HP = t[, 1], SSL_HP = t[, 2], row.names = NULL)
  audit$Retained <- audit$H_ADJ_SSL_HP == 1L & audit$SSL_HP == 1L
  csv(audit, file.path(out, "pairing_audit.csv"))
  if (any(audit$H_ADJ_SSL_HP > 1L & audit$SSL_HP > 0L) ||
      any(audit$SSL_HP > 1L & audit$H_ADJ_SSL_HP > 0L)) stop("Ambiguous exact-site pairing; inspect pairing_audit.csv.")
  m <- m[m$PairID %in% audit$PairID[audit$Retained], , drop = FALSE]
  m <- m[order(m$Patient, m$Site, m$Group), , drop = FALSE]
  if (nrow(m) != 24L || length(unique(m$PairID)) != 12L || length(unique(m$Patient)) != 10L) {
    stop("Exact matched cohort differs from 24 biopsies, 12 site pairs and 10 patients.")
  }
  csv(m, file.path(out, "matched_biopsies.csv"))
  m
}

patient_metadata <- function(m, reference, test, paired) {
  m <- m[m$Group %in% c(reference, test), , drop = FALSE]
  if (paired) {
    shared <- intersect(m$Patient[m$Group == reference], m$Patient[m$Group == test])
    m <- m[m$Patient %in% shared, , drop = FALSE]
  } else if (length(intersect(m$Patient[m$Group == reference], m$Patient[m$Group == test]))) {
    stop("Shared patients detected in a comparison declared independent.")
  }
  m$ProfileID <- paste(m$Patient, m$Group, sep = "__")
  ids <- unique(m$ProfileID)
  p <- m[match(ids, m$ProfileID), c("Patient", "Group", "ProfileID"), drop = FALSE]
  p$Patient <- factor(p$Patient)
  p$Group <- factor(p$Group, levels = c(reference, test))
  p$N_biopsies <- vapply(ids, function(id) sum(m$ProfileID == id), integer(1))
  p$Biopsy_IDs <- vapply(ids, function(id) paste(m$SampleID[m$ProfileID == id], collapse = ";"), character(1))
  p <- p[order(p$Group, p$Patient), , drop = FALSE]
  rownames(p) <- p$ProfileID
  if (paired && any(table(p$Patient, p$Group) != 1L)) stop("Incomplete patient pairing.")
  if (!paired && anyDuplicated(as.character(p$Patient))) stop("Repeated patient in independent comparison.")
  list(biopsies = m, profiles = p)
}

aggregate_counts <- function(x, m, p) {
  if (!all(m$SampleID %in% colnames(x))) stop("Biopsy absent from count matrix.")
  result <- t(rowsum(t(x[, m$SampleID, drop = FALSE]), group = m$ProfileID, reorder = FALSE))
  result <- result[, rownames(p), drop = FALSE]
  stopifnot(identical(colnames(result), rownames(p)), sum(result) == sum(x[, m$SampleID, drop = FALSE]))
  result
}

mean_profiles <- function(x, m, p) {
  result <- rowsum(x[m$SampleID, , drop = FALSE], m$ProfileID, reorder = FALSE)
  result <- result[rownames(p), , drop = FALSE]
  sweep(result, 1, p$N_biopsies, "/")
}

aggregate_features <- function(x, labels) {
  keep <- !is.na(labels)
  if (!any(keep)) stop("No retained taxonomic labels.")
  rowsum(x[keep, , drop = FALSE], group = labels[keep], reorder = TRUE)
}

classified_genus_counts <- function(x, tt) {
  labels <- clean_plain_genus(tt[, "Genus"])
  organelle <- grepl("CHLOROPLAST|MITOCHONDR", toupper(apply(tt, 1, paste, collapse = " ")))
  labels[organelle] <- NA_character_
  z <- aggregate_features(x, labels)
  z[rowSums(z) > 0, , drop = FALSE]
}

species_counts <- function(x, tt, out) {
  if (!all(c("Genus", "Species") %in% colnames(tt))) stop("Species taxonomy is missing.")
  g <- clean_taxonomy_value(tt[, "Genus"])
  s <- clean_taxonomy_value(tt[, "Species"])
  a <- t(vapply(seq_along(g), function(i) classify_species_label(g[i], s[i]), character(2)))
  csv(data.frame(ASV = rownames(tt), Genus = g, Species = s, a, check.names = FALSE),
      file.path(out, "species_assignment_audit.csv"))
  aggregate_features(x, a[, "Species_Key"])
}

plot_theme <- function() {
  ggplot2::theme_classic(base_size = 10, base_family = "sans") +
    ggplot2::theme(axis.text = ggplot2::element_text(colour = "black"),
                   legend.title = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"))
}

colours <- c(H_ADJ_SSL_HP = "#E9C61D", SSL_HP = "#FF0000", H_COL = "#1313EC")

save_plot <- function(p, stem, width = 8, height = 5) {
  device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(paste0(stem, ".pdf"), p, width = width, height = height, device = device, bg = "white")
  ggplot2::ggsave(paste0(stem, ".png"), p, width = width, height = height, dpi = 600, bg = "white")
}

species_plot_font <- function() {
  if (.Platform$OS.type == "windows") {
    try(grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial")), silent = TRUE)
    return("Arial")
  }
  families <- if (requireNamespace("systemfonts", quietly = TRUE)) {
    tryCatch(unique(systemfonts::system_fonts()$family), error = function(e) character(0))
  } else character(0)
  if (any(tolower(families) == "arial")) "Arial" else "sans"
}

make_species_volcano <- function(tab, group_1, group_2) {
  font_family <- species_plot_font()
  alpha <- 0.05
  publication_have_ggrepel <- TRUE
  keep <- is.finite(tab$log2FoldChange) & is.finite(tab$padj)
  plot_data <- data.frame(
    Species = tab$Taxon[keep],
    Effect = tab$log2FoldChange[keep],
    FDR = tab$padj[keep],
    P_value = tab$pvalue[keep],
    Significant = tab$BH_significant[keep],
    stringsAsFactors = FALSE
  )
  plot_data$FDR_plot <- pmax(plot_data$FDR, .Machine$double.xmin)
  plot_data$Negative_log10_FDR <- -log10(plot_data$FDR_plot)
  plot_data$Volcano_status <- factor(
    ifelse(plot_data$Significant %in% TRUE & plot_data$Effect > 0, "Up",
           ifelse(plot_data$Significant %in% TRUE & plot_data$Effect < 0, "Down", "NS")),
    levels = c("Up", "Down", "NS")
  )
  label_data <- plot_data[
    plot_data$Significant %in% TRUE & plot_data$Volcano_status %in% c("Up", "Down") &
      !is.na(plot_data$Species) & nzchar(trimws(plot_data$Species)), , drop = FALSE
  ]
  label_data <- label_data[order(label_data$FDR, label_data$P_value, -abs(label_data$Effect)), , drop = FALSE]
  volcano_colours <- c(
    "Up" = "#A52A2A",
    "Down" = "green4",
    "NS" = "grey75"
  )

  plot_object <- ggplot(
    plot_data,
    aes(x = Effect, y = Negative_log10_FDR)
  ) +
    geom_point(
      aes(colour = Volcano_status),
      size = 1.45,
      alpha = 0.85,
      na.rm = TRUE
    ) +
    scale_colour_manual(
      values = volcano_colours,
      breaks = c("Up", "Down", "NS"),
      limits = c("Up", "Down", "NS"),
      drop = FALSE
    ) +
    scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.10, 0.14))
    ) +
    scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.14))
    ) +
    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dashed",
      linewidth = 0.35,
      colour = "grey50"
    ) +
    geom_hline(
      yintercept = -log10(alpha),
      linetype = "dashed",
      linewidth = 0.35,
      colour = "grey50"
    ) +
    labs(
      title = paste0(group_1, " vs ", group_2),
      x = "log2 fold change",
      y = "-log10(FDR)",
      colour = NULL
    ) +
    coord_cartesian(clip = "off") +
    theme_classic(
      base_family = font_family,
      base_size = 7
    ) +
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.text = element_text(
        family = font_family,
        size = 6,
        colour = "black"
      ),
      legend.key = element_blank(),
      legend.background = element_blank(),
      legend.key.size = grid::unit(3.2, "mm"),
      legend.spacing.x = grid::unit(2, "mm"),
      plot.title = element_text(
        family = font_family,
        face = "bold",
        size = 7,
        colour = "black",
        hjust = 0.5
      ),
      axis.title.x = element_text(
        family = font_family,
        size = 7,
        colour = "black",
        margin = margin(t = 4)
      ),
      axis.title.y = element_text(
        family = font_family,
        size = 7,
        colour = "black",
        margin = margin(r = 4)
      ),
      axis.text.x = element_text(
        family = font_family,
        size = 7,
        colour = "black"
      ),
      axis.text.y = element_text(
        family = font_family,
        size = 7,
        colour = "black"
      ),
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks.length = grid::unit(1.2, "mm"),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(5, 14, 5, 6, "mm")
    ) +
    guides(
      colour = guide_legend(
        override.aes = list(size = 2.2, alpha = 1)
      )
    )

  if (nrow(label_data) > 0L) {
    if (publication_have_ggrepel) {
      plot_object <- plot_object +
        ggrepel::geom_text_repel(
          data = label_data,
          aes(label = Species),
          family = font_family,
          size = 7 / 2.845276,
          colour = "black",
          max.overlaps = Inf,
          force = 5,
          force_pull = 0.08,
          box.padding = 0.30,
          point.padding = 0.20,
          segment.size = 0.20,
          min.segment.length = 0,
          seed = 1,
          show.legend = FALSE
        )
    } else {
      plot_object <- plot_object +
        geom_text(
          data = label_data,
          aes(label = Species),
          family = font_family,
          size = 7 / 2.845276,
          colour = "black",
          vjust = -0.55,
          check_overlap = TRUE,
          show.legend = FALSE
        )
    }
  }

  list(
    plot = plot_object,
    source_data = plot_data,
    labelled_species = label_data
  )
}

run_deseq <- function(counts, normalisation, p, paired, out, normalisation_label) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  counts <- counts[rowSums(counts) > 0, , drop = FALSE]
  normalisation <- normalisation[rowSums(normalisation) > 0, , drop = FALSE]
  stopifnot(identical(colnames(counts), rownames(p)), identical(colnames(normalisation), rownames(p)))
  design <- if (paired) ~ Patient + Group else ~ Group
  dm <- model.matrix(design, p)
  if (qr(dm)$rank != ncol(dm)) stop("DESeq2 design is not full rank.")
  csv(p, file.path(out, "profiles.csv"))
  matrix_csv(counts, file.path(out, "counts.csv"))
  matrix_csv(normalisation, file.path(out, "normalisation_counts.csv"))
  attempt <- capture_fit({
    size_dds <- DESeq2::DESeqDataSetFromMatrix(round(normalisation), p, ~ 1)
    size_dds <- DESeq2::estimateSizeFactors(size_dds, type = "poscounts")
    sf <- DESeq2::sizeFactors(size_dds)
    if (any(!is.finite(sf)) || any(sf <= 0)) stop("Invalid DESeq2 size factors.")
    dds <- DESeq2::DESeqDataSetFromMatrix(round(counts), p, design)
    DESeq2::sizeFactors(dds) <- sf
    fit <- fit_deseq_with_fallback(dds, basename(out))
    dds <- fit$dds
    res <- DESeq2::results(dds, contrast = c("Group", levels(p$Group)[2], levels(p$Group)[1]),
                          alpha = 0.05, independentFiltering = TRUE, cooksCutoff = TRUE)
    list(dds = dds, res = res, strategy = fit$strategy, attempts = fit$attempt_messages)
  })
  write_messages(attempt$warnings, file.path(out, "warnings.txt"))
  if (!is.null(attempt$error)) {
    writeLines(attempt$error, file.path(out, "ERROR.txt"))
    csv(data.frame(Status = "FAILED", Error = attempt$error), file.path(out, "status.csv"))
    return(list(status = "FAILED", results = NULL))
  }
  fit <- attempt$value
  tab <- data.frame(Taxon = rownames(fit$res), as.data.frame(fit$res), row.names = NULL)
  converged <- S4Vectors::mcols(fit$dds)$betaConv
  if (is.null(converged)) converged <- rep(NA, nrow(tab))
  tab$Coefficient_converged <- converged
  tab$Finite_estimate <- is.finite(tab$log2FoldChange) & is.finite(tab$lfcSE) & tab$lfcSE >= 0
  tab$BH_significant <- tab$Finite_estimate & converged %in% TRUE & !is.na(tab$padj) & tab$padj < 0.05
  tab$Holm_P <- p.adjust(tab$pvalue, method = "holm")
  status <- if (any(!converged %in% TRUE)) "REVIEW_REQUIRED" else "COMPLETED"
  csv(tab, file.path(out, "results.csv"))
  csv(tab[tab$BH_significant, , drop = FALSE], file.path(out, "significant.csv"))
  csv(data.frame(ProfileID = colnames(fit$dds), Size_factor = DESeq2::sizeFactors(fit$dds)), file.path(out, "size_factors.csv"))
  matrix_csv(DESeq2::counts(fit$dds, normalized = TRUE), file.path(out, "normalised_counts.csv"))
  saveRDS(fit$dds, file.path(out, "model.rds"))
  write_messages(c(fit$strategy, fit$attempts), file.path(out, "dispersion_fit.txt"))
  csv(data.frame(Status = status, Design = paste(deparse(design), collapse = ""),
                 Normalisation = normalisation_label, N_profiles = nrow(p), N_patients = nlevels(p$Patient),
                 N_taxa = nrow(tab), N_significant = sum(tab$BH_significant)), file.path(out, "status.csv"))
  plot_result <- make_species_volcano(tab, levels(p$Group)[2], levels(p$Group)[1])
  save_plot(plot_result$plot, file.path(out, "volcano"), 183 / 25.4, 110 / 25.4)
  csv(plot_result$source_data, file.path(out, "volcano_SOURCE_DATA.csv"))
  csv(plot_result$labelled_species, file.path(out, "volcano_LABELLED_SPECIES.csv"))
  list(status = status, results = tab)
}

run_ancom <- function(counts, p, paired, out) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  counts <- counts[rowSums(counts) > 0, , drop = FALSE]
  stopifnot(identical(colnames(counts), rownames(p)))
  formula <- if (paired) "Patient + Group" else "Group"
  if (qr(model.matrix(as.formula(paste("~", formula)), p))$rank != ncol(model.matrix(as.formula(paste("~", formula)), p))) stop("ANCOM-BC2 design is not full rank.")
  csv(p, file.path(out, "profiles.csv"))
  matrix_csv(counts, file.path(out, "counts.csv"))
  attempt <- capture_fit({
    if (any(colSums(counts) == 0)) stop("Empty taxonomic profile after filtering.")
    taxonomy <- matrix(rownames(counts), ncol = 1, dimnames = list(rownames(counts), "Feature"))
    phy <- phyloseq::phyloseq(phyloseq::otu_table(round(counts), taxa_are_rows = TRUE),
                              phyloseq::sample_data(p), phyloseq::tax_table(taxonomy))
    ANCOMBC::ancombc2(data = phy, tax_level = NULL, fix_formula = formula, rand_formula = NULL,
                     p_adj_method = "BH", pseudo_sens = TRUE, prv_cut = 0, lib_cut = 0,
                     s0_perc = 0.05, group = "Group", struc_zero = TRUE, neg_lb = FALSE,
                     alpha = 0.05, n_cl = 1L, verbose = TRUE, global = FALSE, pairwise = FALSE,
                     dunnet = FALSE, trend = FALSE,
                     iter_control = list(tol = 1e-2, max_iter = 20, verbose = TRUE),
                     em_control = list(tol = 1e-5, max_iter = 100), lme_control = lme4::lmerControl())
  })
  write_messages(attempt$warnings, file.path(out, "warnings.txt"))
  if (!is.null(attempt$error)) {
    writeLines(attempt$error, file.path(out, "ERROR.txt"))
    csv(data.frame(Status = "FAILED", Error = attempt$error), file.path(out, "status.csv"))
    return(list(status = "FAILED", results = NULL))
  }
  fit <- attempt$value
  saveRDS(fit, file.path(out, "full_output.rds"))
  raw <- as.data.frame(fit$res)
  if (!"taxon" %in% names(raw)) raw$taxon <- rownames(raw)
  csv(raw, file.path(out, "package_results.csv"))
  if (!is.null(fit$zero_ind)) {
    zero <- as.data.frame(fit$zero_ind)
    if (!"taxon" %in% names(zero)) zero$taxon <- rownames(zero)
    csv(zero, file.path(out, "structural_zeros.csv"))
  }
  term <- paste0("Group", levels(p$Group)[2])
  get_col <- function(prefix, logical_column = FALSE) {
    name <- paste0(prefix, "_", term)
    if (!name %in% names(raw) && prefix %in% names(raw)) name <- prefix
    if (!name %in% names(raw)) return(rep(if (logical_column) NA else NA_real_, nrow(raw)))
    if (logical_column) {
      z <- raw[[name]]
      answer <- toupper(as.character(z)) %in% c("TRUE", "1", "T")
      answer[is.na(z)] <- NA
      answer
    } else as.numeric(raw[[name]])
  }
  tab <- data.frame(Taxon = raw$taxon, LFC = get_col("lfc"), SE = get_col("se"),
                    W = get_col("W"), P = get_col("p"), BH_FDR = get_col("q"),
                    Pseudocount_stable = get_col("passed_ss", TRUE))
  if (all(is.na(tab$Pseudocount_stable))) tab$Pseudocount_stable <- get_col("diff_robust", TRUE)
  tab$Finite_estimate <- is.finite(tab$LFC) & is.finite(tab$SE) & tab$SE >= 0 &
    is.finite(tab$W) & is.finite(tab$P) & tab$P >= 0 & tab$P <= 1 &
    is.finite(tab$BH_FDR) & tab$BH_FDR >= 0 & tab$BH_FDR <= 1

  X <- model.matrix(as.formula(paste("~", formula)), p)
  ft <- as.matrix(fit$feature_table)[tab$Taxon, rownames(p), drop = FALSE]
  e <- rep(0, ncol(X)); e[match(term, colnames(X))] <- 1
  info <- t(vapply(seq_len(nrow(ft)), function(i) {
    Xi <- X[ft[i, ] > 0, , drop = FALSE]
    rank <- qr(Xi)$rank
    c(df = nrow(Xi) - rank, estimable = qr(rbind(Xi, e))$rank == rank)
  }, numeric(2)))
  tab$Residual_df <- info[, "df"]
  tab$Contrast_estimable <- info[, "estimable"] == 1
  tab$Test_available <- tab$Finite_estimate & tab$SE > 0 & tab$Residual_df > 0 & tab$Contrast_estimable
  tab$Reporting_P <- ifelse(tab$Test_available, tab$P, NA_real_)
  tab$Reporting_BH_FDR <- ifelse(tab$Test_available, tab$BH_FDR, NA_real_)

  tab$Log2_difference <- tab$LFC / log(2)
  tab$CI_low <- (tab$LFC - 1.96 * tab$SE) / log(2)
  tab$CI_high <- (tab$LFC + 1.96 * tab$SE) / log(2)
  tab$Holm_P <- p.adjust(tab$P, "holm")
  tab$Passes_numeric_criteria <- tab$Test_available & tab$BH_FDR < 0.05 & tab$Pseudocount_stable %in% TRUE
  informative <- grepl("group variable has < 3 categories", attempt$warnings, fixed = TRUE)
  other_warnings <- attempt$warnings[!informative]
  status <- if (length(other_warnings) || any(!tab$Finite_estimate) || anyNA(tab$Pseudocount_stable) || !nrow(tab)) "REVIEW_REQUIRED" else "COMPLETED"
  tab$Fit_status <- status
  tab$Reporting_eligible <- tab$Passes_numeric_criteria
  csv(tab, file.path(out, "results.csv"))
  csv(tab[tab$Reporting_eligible, , drop = FALSE], file.path(out, "reporting_eligible.csv"))
  csv(tab[!tab$Test_available, , drop = FALSE], file.path(out, "unavailable_estimates.csv"))
  csv(data.frame(Status = status, Design = formula, N_taxa_input = nrow(counts), N_model_rows = nrow(tab),
                 N_finite = sum(tab$Finite_estimate), N_test_available = sum(tab$Test_available),
                 N_unavailable = sum(!tab$Test_available),
                 N_numeric_candidates = sum(tab$Passes_numeric_criteria),
                 N_reporting_eligible = sum(tab$Reporting_eligible), N_warnings = length(other_warnings)),
      file.path(out, "status.csv"))
  finite <- tab[tab$Finite_estimate, , drop = FALSE]
  if (nrow(finite)) {
    finite$Taxon <- factor(finite$Taxon, levels = finite$Taxon[order(finite$Log2_difference)])
    g <- ggplot2::ggplot(finite, ggplot2::aes(Log2_difference, Taxon)) +
      ggplot2::geom_segment(ggplot2::aes(x = CI_low, xend = CI_high, yend = Taxon), linewidth = 0.3) +
      ggplot2::geom_point(ggplot2::aes(colour = Passes_numeric_criteria), size = 1.8) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2) +
      ggplot2::scale_colour_manual(values = c("FALSE" = "grey50", "TRUE" = "#A0522D")) +
      ggplot2::labs(title = paste(levels(p$Group)[2], "vs", levels(p$Group)[1]),
                    subtitle = if (status == "COMPLETED") "ANCOM-BC2 finite estimates" else "Diagnostic output: fit requires review",
                    x = expression(log[2]~abundance~difference), y = NULL) + plot_theme()
    save_plot(g, file.path(out, "effects"), 10, max(5, nrow(finite) * 0.17))
  }
  list(status = status, results = tab)
}

prepare_ancom <- function(out) {
  attempt <- capture_fit(patch_ancombc_bias_em_if_needed())
  if (!is.null(attempt$error)) stop("ANCOMBC preflight failed: ", attempt$error)
  csv(attempt$value, file.path(out, "ANCOMBC_patch_status.csv"))
  fun <- get(".bias_em", asNamespace("ANCOMBC"))
  writeLines(c(paste0("ANCOMBC version: ", packageVersion("ANCOMBC")), deparse(fun)),
             file.path(out, "ANCOMBC_bias_em_used.txt"))
}

finish_run <- function(out, started, notes) {
  writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
  writeLines(notes, file.path(out, "analysis_notes.txt"))
  csv(data.frame(Started = as.character(started), Finished = as.character(Sys.time()),
                 Elapsed_minutes = as.numeric(difftime(Sys.time(), started, units = "mins"))),
      file.path(out, "runtime.csv"))
  statuses <- list.files(out, pattern = "^status\\.csv$", full.names = TRUE, recursive = TRUE)
  if (length(statuses)) {
    inventory <- do.call(rbind, lapply(statuses, function(f) {
      z <- read.csv(f, check.names = FALSE)
      data.frame(Analysis = substring(dirname(f), nchar(out) + 2), Status = z$Status[1])
    }))
    csv(inventory, file.path(out, "RUN_STATUS.csv"))
    print(inventory, row.names = FALSE)
  }
  cat("\nSaved outputs to:\n", normalizePath(out), "\n", sep = "")
  invisible(out)
}
publication_wrap_text <- function(x, width) {
  if (is.null(x) || length(x) == 0L) return(x)
  x <- as.character(x)[1]
  if (is.na(x) || !nzchar(trimws(x))) return(x)
  paste(base::strwrap(x, width = width), collapse = "\n")
}

publication_squish <- function(x, range, only.finite = TRUE) {
  if (only.finite) {
    finite <- is.finite(x)
    x[finite] <- pmax(range[1], pmin(range[2], x[finite]))
  } else {
    x <- pmax(range[1], pmin(range[2], x))
  }
  x
}

make_species_heatmap <- function(
    heatmap_data,
    title_text,
    subtitle_text,
    fill_legend_title,
    fill_limit,
    column_labels,
    caption_text,
    axis_text_size = 8
) {
  font_family <- species_plot_font()
  if (is.null(heatmap_data) || nrow(heatmap_data) == 0L) return(NULL)

  title_text <- publication_wrap_text(title_text, width = 62L)
  subtitle_text <- publication_wrap_text(subtitle_text, width = 104L)
  caption_text <- publication_wrap_text(caption_text, width = 132L)

  heatmap_data$Significance_text_colour <- ifelse(
    !is.finite(heatmap_data$Effect), "black",
    ifelse(abs(heatmap_data$Effect) >= 0.62 * fill_limit, "white", "black")
  )

  ggplot(
    heatmap_data,
    aes(x = Column_key, y = Species, fill = Effect)
  ) +
    geom_tile(colour = "white", linewidth = 0.40) +
    geom_text(
      aes(
        label = Significance_stars,
        colour = Significance_text_colour
      ),
      family = font_family,
      fontface = "bold",
      size = 3.2,
      na.rm = TRUE
    ) +
    scale_colour_identity() +
    scale_x_discrete(labels = column_labels, drop = FALSE) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-fill_limit, fill_limit),
      oob = publication_squish,
      na.value = "grey93",
      name = fill_legend_title
    ) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = NULL,
      y = NULL,
      caption = caption_text
    ) +
    theme_classic(base_family = font_family, base_size = 10) +
    theme(
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.title = element_text(
        face = "bold",
        size = 13,
        hjust = 0,
        lineheight = 1.02,
        margin = margin(b = 4)
      ),
      plot.subtitle = element_text(
        size = 9,
        hjust = 0,
        lineheight = 1.08,
        margin = margin(b = 9)
      ),
      plot.caption = element_text(
        size = 7,
        hjust = 0,
        lineheight = 1.08,
        margin = margin(t = 10)
      ),
      axis.text.x = element_text(
        face = "bold",
        angle = 20,
        hjust = 1,
        vjust = 1,
        colour = "black",
        margin = margin(t = 7)
      ),
      axis.text.y = element_text(
        size = axis_text_size,
        colour = "black"
      ),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      plot.margin = margin(18, 20, 20, 20, unit = "pt")
    )
}

run_species <- function(work_dir) {
  message("Species analysis input folder: ", work_dir)
  started <- Sys.time()
  old_options <- options(stringsAsFactors = FALSE, contrasts = c("contr.treatment", "contr.poly"))
  on.exit(options(old_options), add = TRUE)
  set.seed(20260826)
  check_packages(c("phyloseq", "DESeq2", "ANCOMBC", "ggplot2", "ggrepel", "lme4"))
  suppressPackageStartupMessages({ library(phyloseq); library(DESeq2); library(ggplot2) })
  input <- file.path(work_dir, "phyloseq_final_March26_5pct_without7699.rds")
  out <- new_run(work_dir, "Species_rerun", input, "170_species_analysis.R")
  original_bias_em <- get(".bias_em", asNamespace("ANCOMBC"))
  on.exit(assignInNamespace(".bias_em", original_bias_em, ns = "ANCOMBC"), add = TRUE)
  prepare_ancom(out)
  d <- load_input(input, out)
  species <- species_counts(d$counts, d$taxonomy, out)
  species <- species[rowSums(species) > 0, , drop = FALSE]
  matrix_csv(species, file.path(out, "species_biopsy_counts.csv"))
  csv(data.frame(SampleID = colnames(d$counts), ASV_reads = colSums(d$counts),
                  Species_assigned_reads = colSums(species),
                  Fraction_assigned = colSums(species) / colSums(d$counts)),
      file.path(out, "species_resolution.csv"))
  matched <- exact_pairs(d$metadata, out)
  comparisons <- list(
    SSL_HP_vs_H_COL = list(reference = "H_COL", test = "SSL_HP", paired = FALSE),
    H_ADJ_SSL_HP_vs_H_COL = list(reference = "H_COL", test = "H_ADJ_SSL_HP", paired = FALSE),
    SSL_HP_vs_H_ADJ_SSL_HP = list(reference = "H_ADJ_SSL_HP", test = "SSL_HP", paired = TRUE)
  )
  run_one <- function(meta, comparison, settings, parent, prevalence_fraction) {
    location <- file.path(parent, comparison)
    dir.create(location, recursive = TRUE, showWarnings = FALSE)
    bundle <- patient_metadata(meta, settings$reference, settings$test, settings$paired)
    m <- bundle$biopsies
    p <- bundle$profiles
    csv(meta[meta$Group %in% c(settings$reference, settings$test) & !meta$SampleID %in% m$SampleID, , drop = FALSE],
        file.path(location, "biopsies_not_used.csv"))
    if (settings$paired && nlevels(p$Patient) != 10L) stop("Expected 10 paired patients for the lesion-adjacent species comparison.")
    counts <- aggregate_counts(species, m, p)
    asv <- aggregate_counts(d$counts, m, p)
    n_patients <- nlevels(p$Patient)
    detected <- vapply(seq_len(nrow(counts)), function(i) {
      length(unique(as.character(p$Patient)[counts[i, ] > 0]))
    }, integer(1))
    threshold <- max(1L, ceiling(prevalence_fraction * n_patients))
    audit <- data.frame(Taxon = rownames(counts), Detected_patients = detected, Total_patients = n_patients,
      Required_patients = threshold, Retained = detected >= threshold & rowSums(counts) >= 1)
    csv(audit, file.path(location, "prevalence_audit.csv"))
    counts <- counts[audit$Retained, , drop = FALSE]
    if (nrow(counts) < 2L) stop("Fewer than two eligible species for ", comparison)
    de <- run_deseq(counts, asv, p, settings$paired, file.path(location, "DESeq2"),
                    "poscounts on all retained ASVs, summed within patient and tissue group")
    an <- run_ancom(counts, p, settings$paired, file.path(location, "ANCOMBC2"))
    list(deseq2 = de, ancombc2 = an, metadata = p)
  }
  primary <- list()
  for (comparison in names(comparisons)) {
    message("Species primary: ", comparison)
    primary[[comparison]] <- run_one(d$metadata, comparison, comparisons[[comparison]],
                                     file.path(out, "01_primary"), 0)
  }
  for (comparison in names(comparisons)) {
    message("Species 10% patient-prevalence sensitivity: ", comparison)
    run_one(d$metadata, comparison, comparisons[[comparison]], file.path(out, "02_prevalence_10pct"), 0.10)
  }
  message("Species exact-site sensitivity: 10 paired patient profiles from all 12 site pairs")
  run_one(matched, "SSL_HP_vs_H_ADJ_SSL_HP", comparisons$SSL_HP_vs_H_ADJ_SSL_HP,
           file.path(out, "03_exact_site"), 0)
  design_table <- do.call(rbind, lapply(names(primary), function(n) {
    p <- primary[[n]]$metadata
    data.frame(Comparison = n, Paired = comparisons[[n]]$paired,
      Model = if (comparisons[[n]]$paired) "~ Patient + Group" else "~ Group",
      N_profiles = nrow(p), N_patients = nlevels(p$Patient),
      N_reference_patients = sum(p$Group == levels(p$Group)[1]),
      N_test_patients = sum(p$Group == levels(p$Group)[2]))
  }))
  csv(design_table, file.path(out, "comparison_designs.csv"))
  heatmap_dir <- file.path(out, "04_primary_heatmaps")
  dir.create(heatmap_dir)
  for (method in c("deseq2", "ancombc2")) {
    tables <- lapply(names(primary), function(n) {
      a <- primary[[n]][[method]]
      z <- a$results
      if (is.null(z)) return(NULL)
      if (method == "deseq2") {
        data.frame(Taxon = z$Taxon, Comparison = n, Effect = z$log2FoldChange,
                    Q = z$padj, Significant = z$BH_significant, Status = a$status)
      } else {
        data.frame(Taxon = z$Taxon, Comparison = n, Effect = ifelse(z$Test_available, z$Log2_difference, NA_real_),
                    Q = z$BH_FDR, Significant = z$Reporting_eligible, Status = a$status)
      }
    })
    h <- do.call(rbind, tables)
    if (is.null(h) || !nrow(h)) next
    csv(h, file.path(heatmap_dir, paste0(method, "_source.csv")))
    h$Stars <- ifelse(h$Significant & !is.na(h$Q), ifelse(h$Q < 0.001, "***", ifelse(h$Q < 0.01, "**", "*")), "")
    h$Effect[!is.finite(h$Effect)] <- NA_real_
    finite <- abs(h$Effect[is.finite(h$Effect)])
    if (!length(finite)) next
    limit <- unname(quantile(finite, 0.98))
    if (limit == 0) limit <- 1
    h$Display_effect <- pmin(pmax(h$Effect, -limit), limit)
    ranked <- unique(as.character(h$Taxon))
    rank_rows <- lapply(ranked, function(taxon) {
      z <- h[h$Taxon == taxon, , drop = FALSE]
      q <- z$Q[is.finite(z$Q) & is.finite(z$Effect)]
      effects <- abs(z$Effect[is.finite(z$Effect)])
      data.frame(Taxon = taxon, Significant = any(z$Significant %in% TRUE),
                 Best_FDR = if (length(q)) min(q) else Inf,
                 Largest_effect = if (length(effects)) max(effects) else -Inf)
    })
    rank_table <- do.call(rbind, rank_rows)
    rank_table <- rank_table[order(!rank_table$Significant, rank_table$Best_FDR,
                                  -rank_table$Largest_effect, rank_table$Taxon), , drop = FALSE]
    h$Taxon <- factor(h$Taxon, levels = rev(rank_table$Taxon))
    h$Comparison <- factor(h$Comparison, levels = names(comparisons))
    graph <- data.frame(Species = h$Taxon, Column_key = h$Comparison,
                        Effect = h$Effect, Significance_stars = h$Stars)
    column_labels <- setNames(gsub("_vs_", " vs ", names(comparisons), fixed = TRUE),
                              names(comparisons))
    method_label <- if (method == "deseq2") "DESeq2" else "ANCOM-BC2"
    caption <- paste0(
      "Fill shows the ", method_label, " log2 ",
      if (method == "deseq2") "fold change. " else "abundance difference. ",
      "Positive values indicate higher abundance in the first group named in each comparison. ",
      if (method == "ancombc2") "Asterisks require both BH FDR and pseudocount stability: " else "Asterisks indicate ",
      "* BH FDR < 0.05; ** < 0.01; *** < 0.001. ",
      "Grey cells indicate returned model rows without an available effect/test; blank cells have no returned model row. ",
      "The colour scale is truncated symmetrically at the 98th percentile of absolute available effects; exact values are retained in the source-data CSV."
    )
    axis_size <- if (nlevels(h$Taxon) <= 45L) 7.5 else if (nlevels(h$Taxon) <= 65L) 6.7 else if (nlevels(h$Taxon) <= 90L) 6.4 else 5.7
    g <- make_species_heatmap(
      graph,
      title_text = paste0(method_label, ": all model-output species"),
      subtitle_text = "Every species represented in at least one primary model output; rows are ordered by significance and absolute effect size",
      fill_legend_title = if (method == "deseq2") expression("DESeq2 log"[2] * "FC") else expression("ANCOM-BC2 log"[2] * "FC"),
      fill_limit = limit, column_labels = column_labels,
      caption_text = caption, axis_text_size = axis_size
    )
    if (method == "ancombc2") {
      width_in <- if (nlevels(h$Taxon) <= 65L) 11.69 else 16.54
      height_in <- if (nlevels(h$Taxon) <= 65L) 16.54 else 23.39
      device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
      ggplot2::ggsave(file.path(heatmap_dir, paste0(method, ".pdf")), g,
                      width = width_in, height = height_in, device = device, bg = "white")
      ggplot2::ggsave(file.path(heatmap_dir, paste0(method, ".png")), g,
                      width = width_in, height = height_in, dpi = 300, bg = "white")
    } else {
      save_plot(g, file.path(heatmap_dir, method), 10.4,
                max(6.2, 4.0 + 0.235 * nlevels(h$Taxon)))
    }

  }
  finish_run(out, started, c(
    "Input: final selectively subsampled, 5% ASV-prevalence-filtered phyloseq object; patient 7699 excluded.",
    "Existing genus/species labels are audited using the original conservative rules. No taxonomy reassignment is performed.",
    "Accepted species counts are summed within patient and tissue group; species labels remain provisional 16S assignments.",
    "SSL_HP vs H_COL and H_ADJ_SSL_HP vs H_COL use separate patients and model ~ Group.",
    "SSL_HP vs H_ADJ_SSL_HP uses only the 10 patients represented in both groups and model ~ Patient + Group.",
    "Unpaired lesion patients remain included in SSL_HP vs H_COL; they are excluded only from the paired lesion-adjacent comparison.",
    "The primary paired comparison sums all available biopsies within each shared patient and tissue group; sites need not match in that comparison.",
    "The exact-site sensitivity retains all 12 site pairs from 10 patients; counts are summed within patient and group using only the matched sites.",
    "DESeq2 size factors are estimated from all retained ASVs after the same patient aggregation and then applied to species counts.",
    "DESeq2: Wald tests; independent filtering; Cook's filtering; unshrunk log2 fold changes; BH within each comparison.",
    "The primary analysis requires detection in at least one patient. The existing 10% patient-prevalence sensitivity is retained.",
    "No new read-depth or fold-change cutoff is applied. Holm-adjusted P values are exported as a sensitivity analysis.",
    "Heatmap colours are capped at the 98th percentile of absolute effects for display only; source tables retain full values.",
    "ANCOM-BC2 retains the original parameter settings. The internal function and any modification are archived.",
    "ANCOM-BC2 low-taxon-count warnings and numerical warnings cause REVIEW_REQUIRED. They are not suppressed from the saved record.",
    "ANCOM-BC2 test availability uses finite coefficients, positive standard errors, positive residual degrees of freedom and an estimable tissue contrast. Warning diagnostics are retained separately.",
    "Species assignment does not establish strain identity or biological function. No species-specific contamination claim is inferred from these models.",
    "The rerun is complete only as recorded in RUN_STATUS.csv; COMPLETED describes execution, not guaranteed validity of every modelling assumption."
  ))
}

config <- new.env(parent = baseenv())
if (file.exists("config/config.R")) sys.source("config/config.R", envir = config)
work_dir <- if (exists("PATH_WORK", config, inherits = FALSE)) config$PATH_WORK else "outputs"
if (!dir.exists(work_dir)) stop("Input directory does not exist: ", work_dir)
run_species(work_dir)
})
