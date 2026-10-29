#!/usr/bin/env Rscript
# Validate the fixed gastric signature in GSE76987 and GSE198692.
# Normal runs use config/config.R; --self-test checks identifier and sample-selection logic without private inputs.
# Cohort selection, probe mapping and study-specific differential-expression models are recorded in the outputs.


self_test_only <- "--self-test" %in% commandArgs(trailingOnly = TRUE) ||
  identical(Sys.getenv("SSL24_SELF_TEST_ONLY"), "1")
if (self_test_only) {
  input_dir <- tempdir()
} else {
  source("config/config.R")
  input_dir <- PATH_EXTERNAL_VALIDATION
  if (nzchar(Sys.getenv("SSL24_INPUT_DIR"))) input_dir <- Sys.getenv("SSL24_INPUT_DIR")
}

ssl24_external_validation <- function(input_dir) {
  old_options <- options(stringsAsFactors = FALSE,
                         timeout = max(1800, getOption("timeout", 60)))
  on.exit(options(old_options), add = TRUE)
  signature <- c("AQP5", "ANXA10", "MUC5AC", "VSIG1", "TCN1", "MUC6",
                 "PSAPL1", "CLDN18", "SFTA2", "MUCL3", "GATA4", "KLK11",
                 "DUOXA1", "SULT1C2", "IRX3", "AADAC", "FOXQ1", "IGFALS",
                 "NR0B2", "C13ORF46", "FOLR1", "LYPD6B", "TFR2", "SPTSSB")
  aliases <- data.frame(old = c("DPCR1", "LOC100507747"),
                        new = c("MUCL3", "C13ORF46"),
                        entrez = c("135656", "100507747"),
                        probe = c("A_23_P122508", "A_33_P3283400"))
  fail <- function(...) stop(..., call. = FALSE)
  check <- function(ok, ...) if (!isTRUE(ok)) fail(...)
  need_cols <- function(x, columns, label) {
    check(all(columns %in% names(x)), label, " lacks: ",
          paste(setdiff(columns, names(x)), collapse = ", "))
  }
  clean_symbol <- function(x) vapply(as.character(x), function(s) {
    if (is.na(s)) return(NA_character_)
    p <- unique(toupper(trimws(strsplit(s, "\\s*///\\s*|\\s*//\\s*|[;|]",
                                        perl = TRUE)[[1]])))
    p <- p[nzchar(p) & !p %in% c("NA", "---", "UNMAPPED", "NULL")]
    if (length(p) == 1L) p else NA_character_
  }, character(1), USE.NAMES = FALSE)
  harmonise_probes <- function(a) {
    need_cols(a, c("ProbeID", "GeneSymbol", "GPL_Entrez"), "Probe annotation")
    a$GeneSymbol_before_correction <- a$GeneSymbol
    a$identifier_correction <- "none"
    ids <- sub("\\.0+$", "", trimws(as.character(a$GPL_Entrez)))
    for (i in seq_len(nrow(aliases))) {
      anchor <- !is.na(a$ProbeID) & a$ProbeID == aliases$probe[i]
      named <- a$GeneSymbol %in% c(aliases$old[i], aliases$new[i])
      check(any(anchor), "Expected identifier anchor not found: ", aliases$probe[i])
      check(all(!is.na(ids[anchor]) & ids[anchor] == aliases$entrez[i] & named[anchor]),
            "Probe/symbol/NCBI ID disagreement: ", aliases$probe[i])
      check(all(!is.na(ids[named]) & ids[named] == aliases$entrez[i]),
            "Conflicting NCBI ID for ", aliases$old[i], "/", aliases$new[i])
      changed <- named & a$GeneSymbol %in% aliases$old[i]
      a$GeneSymbol[named] <- aliases$new[i]
      a$identifier_correction[changed] <- paste(aliases$old[i], "->", aliases$new[i])
    }
    a
  }

  make_ranks <- function(de) {
    need_cols(de, c("gene_id", "gene_symbol", "stat", "pvalue", "baseMean"), "DE table")
    d <- de[is.finite(de$stat) & is.finite(de$pvalue) &
              !is.na(de$gene_symbol) & nzchar(de$gene_symbol), , drop = FALSE]
    d <- d[order(d$gene_symbol, -abs(d$stat), -d$baseMean, d$gene_id,
                 method = "radix", na.last = TRUE), , drop = FALSE]
    d <- d[!duplicated(d$gene_symbol), , drop = FALSE]
    d <- d[order(-d$stat, d$gene_symbol, method = "radix"), , drop = FALSE]
    r <- stats::setNames(as.numeric(d$stat), d$gene_symbol)
    check(length(r) > 0L && !anyDuplicated(names(r)), "Invalid ranked universe.")
    if (anyDuplicated(r)) r <- r - seq_along(r) *
      (.Machine$double.eps * max(1, max(abs(r))))
    sort(r, decreasing = TRUE)
  }
  enrichment_curve <- function(ranks, genes = signature) {
    hit <- names(ranks) %in% genes
    n <- length(ranks); k <- sum(hit)
    if (!k || k >= n) return(data.frame())
    w <- abs(ranks); total <- sum(w[hit])
    if (total == 0) { w[] <- 1; total <- k }
    es <- cumsum(ifelse(hit, w / total, -1 / (n - k)))
    data.frame(position = 100 * seq_len(n) / n, ES = es,
               hit = hit, gene_symbol = names(ranks))
  }

  geo_urls <- function(accession, filename, kind = "suppl") {
    series <- sub("[0-9]{3}$", "nnn", accession)
    a <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/", series, "/",
                accession, "/", kind, "/", filename)
    b <- if (kind == "soft") paste0(
      "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", accession,
      "&targ=all&form=text&view=full") else paste0(
      "https://www.ncbi.nlm.nih.gov/geo/download/?acc=", accession,
      "&format=file&file=", utils::URLencode(filename, reserved = TRUE))

    unique(c(a, b))
  }
  read_soft_lines <- function(path) {
    con <- file(path, "rb"); sig <- readBin(con, "raw", n = 2L); close(con)
    gz <- identical(as.integer(sig), c(31L, 139L))
    con <- if (gz) gzfile(path, "rt") else file(path, "rt")
    on.exit(close(con), add = TRUE)
    readLines(con, warn = FALSE)
  }
  valid_input <- function(path, type) {
    if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size < 100L)
      return(FALSE)
    isTRUE(tryCatch({
      if (type == "xlsx") {
        nm <- utils::unzip(path, list = TRUE)$Name
        all(c("[Content_Types].xml", "xl/workbook.xml") %in% nm)
      } else if (type == "soft") {
        x <- read_soft_lines(path)
        any(x == "^SERIES = GSE76987") && sum(grepl("^\\^SAMPLE = GSM[0-9]+", x)) == 86L
      } else if (type == "tar") {
        nm <- utils::untar(path, list = TRUE)
        all(vapply(paste0("GSM", 5955144:5955163), function(g)
          sum(grepl(paste0("(^|/)", g, "([^0-9]|$)"), nm)) == 1L, logical(1)))
      } else FALSE
    }, error = function(e) FALSE, warning = function(w) FALSE))
  }
  obtain_input <- function(accession, filename, type, kind = "suppl") {
    candidates <- file.path(c(input_dir, download_dir), filename)
    for (p in candidates[file.exists(candidates)]) {
      check(valid_input(p, type), "Existing file is not a valid ", type, " input: ", p,
            "\nIt has NOT been overwritten. Remove an incomplete download yourself before retrying.")
      return(normalizePath(p, winslash = "/", mustWork = TRUE))
    }
    target <- file.path(download_dir, filename)
    urls <- geo_urls(accession, filename, kind)
    errors <- character()
    for (u in urls) for (attempt in 1:2) {
      temp <- tempfile(pattern = "SSL24_download_", tmpdir = download_dir,
                       fileext = paste0(".", type))
      message("Downloading ", filename, " (attempt ", attempt, ") ...")
      ok <- tryCatch({
        method <- if (capabilities("libcurl")) "libcurl" else "auto"
        status <- utils::download.file(u, temp, mode = "wb", method = method, quiet = FALSE)
        check(identical(as.integer(status), 0L), "Download did not complete.")
        check(valid_input(temp, type), "Response was incomplete or was not the required data file.")
        TRUE
      }, error = function(e) {
        errors <<- c(errors, paste(u, conditionMessage(e), sep = " : "))
        FALSE
      })
      if (ok) {
        check(!file.exists(target) && file.rename(temp, target), "Cannot save ", target)
        return(normalizePath(target, winslash = "/", mustWork = TRUE))
      }
      if (file.exists(temp)) unlink(temp)
      if (attempt == 1L) Sys.sleep(2)
    }
    fail("Could not download ", filename, ". No analysis result was invented.\n",
         paste(unique(errors), collapse = "\n"), "\nThe same script can be rerun after the network recovers.")
  }
  normal_text <- function(x) tolower(gsub("[[:space:]]+", " ", trimws(x)))
  parse_geo_samples <- function(path) {
    x <- read_soft_lines(path)
    starts <- grep("^\\^SAMPLE = GSM[0-9]+", x)
    check(length(starts) == 86L, "GSE76987 SOFT does not contain the expected 86 sample records.")
    records <- lapply(seq_along(starts), function(i) {
      stop_at <- if (i == length(starts)) length(x) else starts[i + 1L] - 1L
      z <- x[starts[i]:stop_at]

      next_entity <- which(grepl("^\\^", z))
      if (length(next_entity) > 1L) z <- z[seq_len(next_entity[2L] - 1L)]
      get_field <- function(key) {
        p <- paste0("!Sample_", key, " = ")
        substring(z[startsWith(z, p)], nchar(p) + 1L)
      }
      chars <- get_field("characteristics_ch1")
      get_characteristic <- function(key, required = TRUE) {
        keys <- normal_text(sub(":.*$", "", chars))
        v <- unique(trimws(sub("^[^:]*:[[:space:]]*", "", chars[keys == key])))
        if (!length(v) && !required) return(NA_character_)
        check(length(v) == 1L, "Missing/conflicting ", key, " in ", z[1])
        v
      }
      desc <- get_field("description")
      title <- unique(get_field("title"))
      check(length(title) == 1L, "Missing/conflicting GEO title in ", z[1])
      data.frame(GSM = sub("^\\^SAMPLE = ", "", z[1]), title = title,
        cohort = get_characteristic("patient cohort", required = FALSE),
        tissue = get_characteristic("tissue", required = FALSE),
        patient_id = get_characteristic("patient id", required = FALSE),
        description = paste(desc, collapse = " || "), stringsAsFactors = FALSE)
    })
    m <- do.call(rbind, records)
    check(!anyDuplicated(m$GSM), "Duplicate GEO sample records.")
    co <- normal_text(m$cohort); ti <- normal_text(m$tissue)
    m$group <- NA_character_
    m$group[which(co == "sporadic" & ti == "sessile serrated adenoma/polyp")] <- "SSL"
    m$group[which(co == "control" & ti == "control right colon")] <- "Control_Right"
    m$included <- !is.na(m$group)
    m$reason <- ifelse(m$included, "Included by tissue/cohort rule before testing",
      ifelse(is.na(co) | is.na(ti), "Missing cohort/tissue metadata: not selected",
        ifelse(co == "serrated polyposis", "SPS specimens excluded: documented repeated lesions",
          "Outside the fixed sporadic-SSL versus control-right-colon contrast")))
    check(sum(m$group == "SSL", na.rm = TRUE) == 9L &&
            sum(m$group == "Control_Right", na.rm = TRUE) == 10L,
          "Metadata did not identify exactly 9 sporadic SSLs and 10 control-right samples.\n",
          "Observed SSL = ", sum(m$group == "SSL", na.rm = TRUE),
          "; control right = ", sum(m$group == "Control_Right", na.rm = TRUE),
          ". The script will not guess or substitute samples.")
    m
  }

  expected_rnaseq_samples <- function() data.frame(
    GSM = c(paste0("GSM", 2042070:2042078), paste0("GSM", 2042124:2042133)),
    author_label = c("110SSP_A", "100SSP_S", "101SSP_A", "118SSP_C1",
      "203SSP_C", "12SSP_T", "118SSP_C2", "127SSP_A", "133SSP_R",
      "71NA", "72NA", "74NA", "81NA", "84NA", "85NA", "86NA", "87NA", "88NA", "89NA"),
    expected_group = c(rep("SSL", 9L), rep("Control_Right", 10L)),
    stringsAsFactors = FALSE)
  repeated_pair <- c("GSM2042073", "GSM2042076")
  retained_lesion <- "GSM2042073"
  excluded_lesion <- "GSM2042076"
  linkage_status <- "INFERRED_FROM_AUTHOR_LABELS_AND_SUPPLEMENT_NOT_EXPLICIT_GSM_CROSSWALK"
  restriction_note <- c(
    "Source fact: Kanth Suppl Table 1a reports 9 sporadic SSL lesions from 8 patients.",
    "Source fact: supplement specimens 8 and 9 come from one patient; both are cecal lesions.",
    "GEO description: GSM2042073 / SSA/P-16 is 118SSP_C1.",
    "GEO description: GSM2042076 / SSA/P-19 is 118SSP_C2.",
    "INFERENCE: these two uniquely repeated 118SSP labels denote that repeated patient.",
    "The supplement does not explicitly link its row numbers to GEO accessions.",
    "Rule: retain C1 and exclude C2; fixed before calculating this eight-specimen result.",
    "No choice is made from expression, library size, fold changes, NES or P values.",
    "Do not describe the inferred label-to-patient cross-reference as author-verified.",
    "If the inference is wrong, independence of the eight retained specimens is not established.",
    "No patient IDs, sizes, sex, ages or lesion locations are assigned by guessed row order.",
    "Separate lesions are not summed as if they were technical sequencing replicates.")
  restrict_rnaseq_samples <- function(m) {
    need_cols(m, c("GSM", "title", "description", "group", "included", "reason"),
              "GSE76987 selection metadata")
    check(!anyNA(m$GSM) && !anyDuplicated(m$GSM), "Invalid or duplicate GEO accessions.")
    check(!anyNA(m$included), "Missing source inclusion flags.")
    e <- expected_rnaseq_samples()
    check(setequal(m$GSM[m$included], e$GSM),
          "The source's 9-SSL/10-control selection differs from the audited dataset.")
    j <- match(e$GSM, m$GSM)
    check(!anyNA(j) && identical(as.character(m$group[j]), e$expected_group),
          "GSE76987 GSM-to-group assignment changed.")
    ok <- vapply(seq_len(nrow(e)), function(i) {
      e$author_label[i] %in% sample_tokens(m$GSM[j[i]], m$title[j[i]], m$description[j[i]])
    }, logical(1))
    check(all(ok), "GEO GSM-to-author-label mapping changed: ",
          paste(e$GSM[!ok], collapse = ", "), ". No positional deletion was attempted.")

    if ("patient_id" %in% names(m)) {
      p <- trimws(as.character(m$patient_id[match(repeated_pair, m$GSM)]))
      known <- !is.na(p) & nzchar(p) & !toupper(p) %in% c("NA", "N/A", "UNKNOWN")
      check(!all(known) || length(unique(p)) == 1L,
            "Explicit patient identifiers contradict the inferred 118SSP linkage.")
    }
    m$eligible_before_single_lesion_rule <- m$included
    m$verified_GEO_author_label <- NA_character_
    m$verified_GEO_author_label[j] <- e$author_label
    m$patient_linkage_status <- "No new patient identity assigned"
    m$patient_linkage_status[m$GSM %in% repeated_pair] <- linkage_status
    m$included[m$GSM == excluded_lesion] <- FALSE
    m$reason[m$GSM == excluded_lesion] <- paste0(
      "Excluded by fixed C1-before-C2 rule for INFERRED repeated 118SSP patient; ",
      "not an expression/QC exclusion and not a confirmed GSM-to-patient crosswalk")
    m$reason[m$GSM == retained_lesion] <-
      "Retained C1 by fixed author-label rule for INFERRED repeated 118SSP patient"
    check(sum(m$included) == 18L &&
            sum(m$included & m$group == "SSL", na.rm = TRUE) == 8L &&
            sum(m$included & m$group == "Control_Right", na.rm = TRUE) == 10L &&
            retained_lesion %in% m$GSM[m$included] &&
            !excluded_lesion %in% m$GSM[m$included],
          "Eight-SSL/ten-control restriction failed.")
    m
  }

  sample_tokens <- function(gsm, title, description) {
    parts <- unlist(strsplit(description, "[[:space:]|,;:]+", perl = TRUE))
    codes <- parts[grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", parts) &
                     grepl("[A-Za-z]", parts) & grepl("[0-9]", parts) &
                     !grepl("\\.(xlsx?|gz|txt|csv)$", parts, ignore.case = TRUE)]
    unique(c(gsm, title, codes))
  }
  sample_key <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))
  regex_escape <- function(x) {
    ch <- strsplit(x, "", fixed = TRUE)[[1L]]
    special <- ch %in% c("\\", ".", "^", "$", "|", "(", ")", "[", "]", "{", "}", "*", "+", "?")
    paste0(ifelse(special, paste0("\\", ch), ch), collapse = "")
  }
  matches_sample <- function(header, tokens) {
    body <- trimws(sub("^counts?[[:space:]_.:-]+", "", header, ignore.case = TRUE))
    if (sample_key(body) %in% sample_key(tokens)) return(TRUE)
    any(vapply(tokens, function(t) grepl(paste0("(^|[^A-Za-z0-9])",
      regex_escape(t), "([^A-Za-z0-9]|$)"), body, ignore.case = TRUE, perl = TRUE), logical(1)))
  }
  read_geo_workbook <- function(path, metadata, audit_dir) {

    sheets <- readxl::excel_sheets(path)
    allowed <- normal_text(sheets) %in% c("analyzed genes", "analysed genes", "all genes", "all analyzed genes")
    check(sum(allowed) == 1L, "Cannot uniquely identify the unselected USeq gene table.\n",
          "Workbook sheets: ", paste(sheets, collapse = ", "),
          "\nRefusing a significance-filtered or guessed sheet.")
    sheet <- sheets[allowed]
    preview <- as.data.frame(readxl::read_excel(path, sheet = sheet, col_names = FALSE,
                    col_types = "text", n_max = 30L, .name_repair = "minimal"))
    count_header <- function(v) grepl("^counts?[[:space:]_.:-]+[^[:space:]]", trimws(v), ignore.case = TRUE)
    nh <- apply(as.matrix(preview), 1L, function(v) sum(!is.na(v) & count_header(v)))
    hi <- which(nh >= 19L)
    check(length(hi) == 1L, "Cannot identify a unique header with explicit raw Counts columns.\n",
          "No FPKM, RLog or inferred integer-valued columns will be used.")
    tab <- as.data.frame(readxl::read_excel(path, sheet = sheet, skip = hi - 1L,
                       col_types = "text", .name_repair = "minimal"), check.names = FALSE)
    headers <- trimws(names(tab)); names(tab) <- headers
    cc <- which(!is.na(headers) & count_header(headers))

    tokens <- lapply(seq_len(nrow(metadata)), function(i) sample_tokens(
      metadata$GSM[i], metadata$title[i], metadata$description[i]))
    owners <- lapply(headers[cc], function(h) which(vapply(tokens,
      function(t) matches_sample(h, t), logical(1))))
    audit <- data.frame(excel_column = cc, header = headers[cc],
      n_metadata_matches = lengths(owners),
      matched_GSM = vapply(owners, function(j) paste(metadata$GSM[j], collapse = ";"), character(1)))
    write_csv(audit, file.path(audit_dir, "RAW_COUNT_COLUMN_MAPPING.csv"))
    sel <- metadata[metadata$included, , drop = FALSE]

    sel$group <- factor(sel$group, levels = c("Control_Right", "SSL"))
    sel <- sel[order(sel$group, sel$GSM), , drop = FALSE]
    cols <- vapply(sel$GSM, function(g) {
      j <- which(audit$n_metadata_matches == 1L & audit$matched_GSM == g)
      check(length(j) == 1L, "No unique raw count column for ", g,
            "; inspect RAW_COUNT_COLUMN_MAPPING.csv. No positional matching was performed.")
      cc[j]
    }, integer(1))
    check(!anyDuplicated(cols), "Two samples mapped to the same raw count column.")
    sel$count_column <- headers[cols]
    sel$sample_id <- sel$GSM
    check(all(grepl("RightColonProcessed\\.xlsx", sel$description, ignore.case = TRUE)),
          "Selected GEO sample metadata do not all name the right-colon analysis workbook.")

    is_ensg <- function(x) !is.na(x) & grepl("^ENSG[0-9]+(\\.[0-9]+)?$", trimws(x))
    noncount <- setdiff(seq_len(ncol(tab)), cc)
    scores <- vapply(noncount, function(j) {
      v <- tab[[j]]; ok <- !is.na(v) & nzchar(trimws(v))
      if (!any(ok)) 0 else mean(is_ensg(v[ok]))
    }, numeric(1))
    gene_cols <- noncount[scores >= 0.90]
    check(length(gene_cols) >= 1L, "No credible Ensembl gene-ID column was found in the full gene table.")
    if (length(gene_cols) > 1L) check(all(vapply(gene_cols[-1L], function(j)
      identical(trimws(tab[[j]]), trimws(tab[[gene_cols[1L]]])), logical(1))),
      "Conflicting Ensembl identifier columns; refusing to guess.")
    gc <- gene_cols[1L]; ids <- trimws(tab[[gc]])
    blank <- is.na(ids) | !nzchar(ids)
    empty_counts <- apply(tab[, cols, drop = FALSE], 1L, function(v) all(is.na(v) | !nzchar(trimws(v))))
    discard <- blank & empty_counts
    if (any(discard)) { tab <- tab[!discard, , drop = FALSE]; ids <- ids[!discard] }
    check(all(is_ensg(ids)), "Malformed/missing gene ID in a non-empty count row.")
    ids <- sub("\\.[0-9]+$", "", ids)
    check(!anyDuplicated(ids), "Duplicated Ensembl genes in the source sheet. Counts are not blindly summed.")
    counts <- as.matrix(tab[, cols, drop = FALSE])
    suppressWarnings(storage.mode(counts) <- "numeric")
    check(all(is.finite(counts)) && all(counts >= 0) &&
            all(abs(counts - round(counts)) <= 1e-8) && max(counts) <= .Machine$integer.max,
          "Selected raw count columns contain missing, negative or non-integer values. No rounding was applied.")
    storage.mode(counts) <- "integer"
    rownames(counts) <- ids; colnames(counts) <- sel$sample_id
    check(nrow(counts) >= 1000L && all(colSums(counts) > 0), "Incomplete gene universe or zero-size library.")
    rownames(sel) <- sel$sample_id
    write_csv(sel, file.path(audit_dir, "SELECTED_SAMPLE_DESIGN.csv"))
    write_csv(data.frame(workbook = basename(path), sheet = sheet, header_row = hi,
      gene_column = headers[gc], n_source_genes = nrow(counts), n_selected_samples = ncol(counts),
      row_selection = "All deposited rows of Analyzed Genes; no DE/P/CPM selection"),
      file.path(audit_dir, "WORKBOOK_IMPORT_AUDIT.csv"))
    list(counts = counts, metadata = sel)
  }

  logical_tests <- function() {
    check(length(signature) == 24L && !anyDuplicated(signature), "Signature altered.")
    m <- data.frame(ProbeID = aliases$probe, GeneSymbol = aliases$old,
                    GPL_Entrez = paste0(aliases$entrez, ".0"))
    check(identical(harmonise_probes(m)$GeneSymbol, aliases$new), "Alias test failed.")
    bad <- m; bad$GPL_Entrez[1] <- "999"
    check(inherits(try(harmonise_probes(bad), silent = TRUE), "try-error"),
          "Conflicting annotation was not rejected.")
    d <- data.frame(gene_id = letters[1:5], gene_symbol = c("A", "B", "C", "A", "D"),
                    stat = c(2, 0, 9, 1, -1), pvalue = c(0.1, 1, NA, 0.4, 0.3),
                    baseMean = 1)
    check(identical(names(make_ranks(d)), c("A", "B", "D")), "Ranking test failed.")
    check(isTRUE(all.equal(stats::p.adjust(c(0.001, 1e-8), "BH", n = 2),
                           c(0.001, 2e-8))), "BH test failed.")
    cc <- enrichment_curve(stats::setNames(c(3, 2, -1, -2), LETTERS[1:4]), "A")
    check(abs(tail(cc$ES, 1)) < 1e-12 && abs(max(cc$ES) - 1) < 1e-12,
          "Enrichment-curve test failed.")
    tk <- sample_tokens("GSM2042070", "SSA/P - 13", "RightColonProcessed.xlsx || 110SSP_A")
    check(matches_sample("Counts 110SSP_A", tk) &&
            matches_sample("Counts Sporadic_110SSP_A_0", tk) &&
            !matches_sample("Counts 1110SSP_A", tk) &&
            !matches_sample("Counts 110SSP_AA", tk), "Exact sample-label matching test failed.")
    e <- expected_rnaseq_samples()
    sm <- data.frame(GSM = e$GSM, title = e$GSM,
      description = paste("RightColonProcessed.xlsx ||", e$author_label),
      group = e$expected_group, included = TRUE, reason = "eligible", patient_id = NA_character_)
    chosen <- restrict_rnaseq_samples(sm)
    check(sum(chosen$included) == 18L &&
            identical(chosen$GSM[!chosen$included], excluded_lesion),
          "Fixed single-lesion selection test failed.")

    reversed <- restrict_rnaseq_samples(sm[nrow(sm):1L, , drop = FALSE])
    check(setequal(chosen$GSM[chosen$included], reversed$GSM[reversed$included]),
          "Single-lesion restriction depends on metadata row order.")
    bad <- sm; bad$description[bad$GSM == excluded_lesion] <- "RightColonProcessed.xlsx || 999SSP_C2"
    check(inherits(try(restrict_rnaseq_samples(bad), silent = TRUE), "try-error"),
          "Mismatched repeated-specimen author label was not rejected.")
    bad <- sm; bad$patient_id[match(repeated_pair, bad$GSM)] <- c("PATIENT_A", "PATIENT_B")
    check(inherits(try(restrict_rnaseq_samples(bad), silent = TRUE), "try-error"),
          "Contradictory explicit patient metadata was not rejected.")
    check(inherits(try(restrict_rnaseq_samples(sm[-1L, , drop = FALSE]), silent = TRUE), "try-error"),
          "Incomplete source sample set was not rejected.")
    invisible(TRUE)
  }
  logical_tests()
  if ("--self-test" %in% commandArgs(trailingOnly = TRUE) ||
      identical(Sys.getenv("SSL24_SELF_TEST_ONLY"), "1")) {
    message("Logical self-tests passed; no real data or fitted results were tested.")
    return(invisible(NULL))
  }

  check(dir.exists(input_dir), "Input directory does not exist. Edit input_dir: ", input_dir)
  input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
  cran <- c("data.table", "ggplot2", "ggrepel", "pheatmap", "cowplot", "R.utils", "statmod", "readxl")
  bioc <- c("DESeq2", "AnnotationDbi", "org.Hs.eg.db", "SummarizedExperiment",
            "limma", "GEOquery", "fgsea", "singscore")
  missing <- function(p) p[!vapply(p, requireNamespace, logical(1), quietly = TRUE)]
  check(!length(missing(c(cran, bioc))), "Required packages missing: ",
        paste(missing(c(cran, bioc)), collapse = ", "))
  suppressPackageStartupMessages(library(ggplot2))
  out <- file.path(input_dir, paste0("SSL24_76987_198692_",
                                    format(Sys.time(), "%Y%m%d_%H%M%S")))
  check(!dir.exists(out) && dir.create(out, recursive = TRUE), "Cannot create: ", out)
  mk <- function(...) { p <- file.path(...); dir.create(p, recursive = TRUE, showWarnings = FALSE); p }
  rdir <- mk(out, "01_GSE76987"); sdir <- mk(out, "02_GSE198692")
  fdir <- mk(out, "03_FINAL_FIGURE"); adir <- mk(out, "00_METADATA")
  completed <- FALSE
  on.exit(if (!completed) writeLines("Execution did not finish; do not use partial outputs as final.",
                                    file.path(out, "RUN_INCOMPLETE.txt")), add = TRUE)
  log_path <- file.path(adir, "WARNINGS.txt")
  writeLines("Warnings encountered in this execution:", log_path)
  write_csv <- function(x, path) {
    x <- as.data.frame(x)
    for (nm in names(x)) if (is.list(x[[nm]]))
      x[[nm]] <- vapply(x[[nm]], paste, collapse = ";", FUN.VALUE = character(1))
    data.table::fwrite(x, path, na = "NA")
  }
  read_csv <- function(path) as.data.frame(data.table::fread(path, check.names = FALSE))
  matrix_csv <- function(x, path, id = "gene_symbol") {
    d <- data.frame(rownames(x), x, check.names = FALSE); names(d)[1] <- id
    write_csv(d, path)
  }
  download_dir <- mk(input_dir, "SSL24_DOWNLOADS")
  manifest <- data.frame(file = character(), size_bytes = numeric(), md5 = character())
  record_input <- function(path) {
    row <- data.frame(file = normalizePath(path, winslash = "/", mustWork = TRUE),
                       size_bytes = file.info(path)$size, md5 = unname(tools::md5sum(path)))
    manifest <<- unique(rbind(manifest, row))
    write_csv(manifest, file.path(adir, "RAW_INPUT_MANIFEST.csv"))
  }
  writeLines(c("Combined pipeline revision: single-lesion restriction, 2026-09-15",
    "RNAseq: GSM2042076 excluded under inferred shared-118SSP linkage; never called author-verified.",
    "The ten-pair microarray model, preprocessing and GSEA settings are unchanged.",
    "Previous nine-specimen RNAseq outputs are not reused as new eight-specimen results."),
    file.path(adir, "PIPELINE_REVISION.txt"))
  writeLines(signature, file.path(adir, "FIXED_24_GENE_SIGNATURE.txt"))
  write_csv(aliases, file.path(adir, "VERIFIED_IDENTIFIER_CORRECTIONS.csv"))
  write_csv(data.frame(package = c(cran, bioc), version = vapply(c(cran, bioc),
    function(p) as.character(utils::packageVersion(p)), character(1))),
    file.path(adir, "PACKAGE_VERSIONS.csv"))
  set.seed(20260827)

  fmt <- function(x) if (!length(x) || !is.finite(x)) "NA" else sprintf("%.3g", x)
  save_plot <- function(p, stem, w = 8, h = 6) {
    ggplot2::ggsave(paste0(stem, ".pdf"), p, width = w, height = h,
                    device = grDevices::pdf, useDingbats = FALSE, bg = "white")
    ggplot2::ggsave(paste0(stem, ".png"), p, width = w, height = h,
                    dpi = 600, bg = "white")
  }
  draw_base <- function(stem, draw, w = 8, h = 6) {
    for (ext in c("pdf", "png")) {
      if (ext == "pdf") grDevices::pdf(paste0(stem, ".pdf"), width = w, height = h)
      else grDevices::png(paste0(stem, ".png"), width = w, height = h, units = "in", res = 400)
      tryCatch(draw(), finally = grDevices::dev.off())
    }
  }
  heat <- function(x, stem, title, annotation = NULL, cluster_cols = FALSE,
                   row_z = FALSE, w = 10, h = 8, breaks = NULL) {
    if (row_z) x <- t(scale(t(x)))
    x <- x[apply(x, 1, function(v) all(is.finite(v))), , drop = FALSE]
    if (nrow(x) < 2L || ncol(x) < 2L) return(invisible(NULL))
    args <- list(mat = x, cluster_rows = TRUE, cluster_cols = cluster_cols,
                 annotation_col = annotation, main = title, border_color = NA,
                 fontsize_row = 8, fontsize_col = 8, angle_col = "45", silent = TRUE)
    if (!is.null(breaks)) {
      args$breaks <- breaks
      args$color <- grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(length(breaks)-1L)
    }
    ph <- do.call(pheatmap::pheatmap, args)
    draw_base(stem, function() { grid::grid.newpage(); grid::grid.draw(ph$gtable) }, w, h)
  }
  signature_table <- function(de, ranks) {
    d <- de[order(de$gene_symbol, -abs(de$stat), -de$baseMean, de$gene_id,
                   method = "radix", na.last = TRUE), , drop = FALSE]
    d <- d[!duplicated(d$gene_symbol), , drop = FALSE]
    a <- data.frame(gene_symbol = signature)
    other <- setdiff(names(d), "gene_symbol")
    a <- cbind(a, d[match(signature, d$gene_symbol), other, drop = FALSE])
    rownames(a) <- NULL
    a$in_ranked_universe <- signature %in% names(ranks)
    a$FDR_below_0_05 <- is.finite(a$padj) & a$padj < 0.05
    a
  }
  run_gsea <- function(ranks, label, seed, directory, stem) {
    k <- intersect(signature, names(ranks))
    result <- data.frame(comparison = label, ES = NA_real_, NES = NA_real_,
                          pval = NA_real_, padj_within_single_pathway = NA_real_,
                          log2err = NA_real_, size = length(k),
                          n_signature_genes_in_ranked_universe = length(k),
                          n_fixed_signature_genes = length(signature),
                          n_ranked_genes = length(ranks), leading_edge_genes = "",
                          missing_signature_genes = paste(setdiff(signature, k), collapse = ";"),
                          analysis_status = "NOT_RUN_FEWER_THAN_5_SIGNATURE_GENES")
    if (length(k) >= 5L) {
      check(length(k) < length(ranks) && sum(abs(ranks[names(ranks) %in% k])) > 0,
            "Uninformative ranked signature in ", label)
      set.seed(seed)
      fg <- fgsea::fgseaMultilevel(
        pathways = list(Serrated_Lesion_Gastric_Metaplasia_24_gene_signature = signature),
        stats = ranks, minSize = 5L, maxSize = 500L, eps = 0,
        scoreType = "std", nproc = 1L)
      check(nrow(fg) == 1L && fg$size[1] == length(k), "Unexpected fgsea coverage: ", label)
      for (nm in c("ES", "NES", "pval", "log2err")) result[[nm]] <- fg[[nm]][1]
      result$padj_within_single_pathway <- fg$padj[1]
      result$leading_edge_genes <- paste(fg$leadingEdge[[1]], collapse = ";")
      result$analysis_status <- if (all(is.finite(c(result$ES, result$NES, result$pval))))
        "RUN" else "FAILED_NONFINITE_GSEA_RESULT"
    }
    write_csv(result, file.path(directory, paste0(stem, "__GSEA.csv")))
    check(result$analysis_status != "FAILED_NONFINITE_GSEA_RESULT",
          "fgsea returned non-finite values for ", label, "; inspect the saved GSEA CSV and warnings.")
    write_csv(data.frame(rank = seq_along(ranks), gene_symbol = names(ranks),
                         ranking_statistic = unname(ranks), in_signature = names(ranks) %in% signature),
              file.path(directory, paste0(stem, "__RANKED_GENES.csv")))
    write_csv(data.frame(gene_symbol = signature, represented = signature %in% names(ranks)),
              file.path(directory, paste0(stem, "__COVERAGE.csv")))
    list(result = result, ranks = ranks)
  }
  gsea_plot <- function(g, title) {
    if (g$result$analysis_status != "RUN") return(ggplot() + theme_void() +
      annotate("text", x = 0, y = 0, label = paste(title, g$result$analysis_status, sep = "\n")))
    fgsea::plotEnrichment(signature, g$ranks) + theme_classic(base_size = 11) +
      labs(title = title, subtitle = paste0("NES = ", fmt(g$result$NES),
        "; nominal GSEA P = ", fmt(g$result$pval), "; coverage = ", g$result$size, "/24"),
        x = "Genes ordered by test statistic", y = "Running enrichment score")
  }
  volcano <- function(de, stem, title, adjusted = FALSE) {
    d <- de; d$in_signature <- d$gene_symbol %in% signature
    pp <- if (adjusted) d$padj else d$pvalue
    d$y <- -log10(pmax(pp, .Machine$double.xmin))
    label <- d[d$in_signature & is.finite(d$y) & is.finite(d$log2FoldChange), , drop = FALSE]
    p <- ggplot(d, aes(log2FoldChange, y)) +
      geom_point(aes(shape = in_signature), size = 1.2, alpha = 0.55, na.rm = TRUE) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
      ggrepel::geom_text_repel(data = label, aes(label = gene_symbol),
                              size = 2.8, max.overlaps = Inf, seed = 1, show.legend = FALSE) +
      theme_classic(base_size = 11) + labs(title = title, shape = "SSL24 gene",
        x = "log2 fold change", y = if (adjusted) "-log10(BH-adjusted P)" else "-log10(nominal P)")
    save_plot(p, stem, 8.5, 6.5)
  }
  save_pca <- function(coords, annotation, stem, title, variance) {
    d <- data.frame(sample_id = rownames(coords), PC1 = coords[, 1], PC2 = coords[, 2],
                     group = annotation$group[match(rownames(coords), annotation$sample_id)])
    write_csv(d, paste0(stem, "_COORDINATES.csv"))
    p <- ggplot(d, aes(PC1, PC2, shape = group)) + geom_point(size = 3) +
      ggrepel::geom_text_repel(aes(label = sample_id), size = 3, seed = 1, max.overlaps = Inf) +
      theme_classic(base_size = 11) + labs(title = title, shape = NULL,
        x = sprintf("PC1 (%.1f%%)", variance[1]), y = sprintf("PC2 (%.1f%%)", variance[2]))
    save_plot(p, stem)
  }

  reuse_microarray <- function() {
    roots <- list.dirs(input_dir, recursive = FALSE, full.names = TRUE)
    roots <- roots[grepl("^(SSL24_EXTERNAL_VALIDATION_FINAL_|SSL24_76987_198692_)", basename(roots)) & roots != out]
    if (!length(roots)) return(NULL)
    roots <- roots[order(file.info(roots)$mtime, decreasing = TRUE)]
    failures <- character()
    for (root in roots) {
      src <- file.path(root, "02_GSE198692")
      rel <- c("01_QC_AND_METADATA/SAMPLE_AND_RAW_FILE_MAP.csv",
        "01_QC_AND_METADATA/PROBE_ANNOTATION_AND_FILTERING_AUDIT.csv",
        "01_QC_AND_METADATA/PAIRED_LIMMA_DESIGN_MATRIX.csv",
        "02_PAIRED_LIMMA_DE/SSL_VS_PAIRED_NORMAL_ALL_GENES.csv",
        "03_SSL24_GSEA/SSL_VS_PAIRED_NORMAL__RANKED_GENES.csv",
        "03_SSL24_GSEA/SSL_VS_PAIRED_NORMAL__GSEA.csv",
        "04_SSL24_SAMPLE_SCORES/PAIRED_SCORE_TESTS.csv",
        "06_LEAVE_ONE_PAIR_OUT/LEAVE_ONE_PAIR_OUT_SUMMARY.csv")
      if (!file.exists(file.path(root, "RUN_COMPLETED.txt")) ||
          !all(file.exists(file.path(src, rel)))) next
      result <- tryCatch({
        check(identical(readLines(file.path(root, "00_METADATA/FIXED_24_GENE_SIGNATURE.txt")), signature),
              "Fixed signature differs from requested signature.")
        m <- read_csv(file.path(src, rel[1])); a <- read_csv(file.path(src, rel[2]))
        expected_pairs <- c("01", "03", "04", "08", "10", "12", "13", "14", "17", "21")
        expected_samples <- paste0(rep(expected_pairs, each = 2), rep(c("N", "T"), 10))
        check(identical(as.character(m$GSM), paste0("GSM", 5955144:5955163)) &&
                identical(as.character(m$Sample), expected_samples) &&
                identical(as.character(m$Group), rep(c("Adjacent_Normal", "SSL"), 10)), "Cache sample map mismatch.")
        m$PairID <- factor(rep(expected_pairs, each = 2), levels = expected_pairs)
        m$Group <- factor(m$Group, levels = c("Adjacent_Normal", "SSL"))
        dd <- read_csv(file.path(src, rel[3])); design <- stats::model.matrix(~ PairID + Group, m)
        check(identical(names(dd)[-1L], colnames(design)) &&
          isTRUE(all.equal(unname(as.matrix(dd[, -1L, drop = FALSE])), unname(design), check.attributes = FALSE)),
          "Saved paired design does not match the expected model.")

        a <- harmonise_probes(a)
        check(all(aliases$new %in% a$GeneSymbol[a$retained]), "Corrected gene identifiers not retained.")
        d <- read_csv(file.path(src, rel[4])); ranks <- make_ranks(d)
        rt <- read_csv(file.path(src, rel[5])); g <- read_csv(file.path(src, rel[6]))
        check(identical(as.character(rt$gene_symbol), names(ranks)) &&
          isTRUE(all.equal(as.numeric(rt$ranking_statistic), unname(ranks), tolerance = 1e-12)),
          "Saved GSEA ranks disagree with the DE table.")
        overlap <- intersect(signature, names(ranks))
        check(nrow(g) == 1L && g$analysis_status == "RUN" &&
          all(is.finite(c(g$ES, g$NES, g$pval))) && g$pval > 0 && g$pval <= 1 &&
          g$size == length(overlap) && length(overlap) == 23L &&
          identical(setdiff(signature, overlap), "SPTSSB") &&
          !any(aliases$old %in% names(ranks)), "Cache signature/GSEA mismatch.")
        cv <- enrichment_curve(ranks)
        es <- if (max(cv$ES) > -min(cv$ES)) max(cv$ES) else min(cv$ES)
        check(abs(es - g$ES) < 1e-8, "Saved enrichment score disagrees with ranks.")
        list(gsea = list(result = g, ranks = ranks), metadata = m,
             score_tests = read_csv(file.path(src, rel[7])), lopo = read_csv(file.path(src, rel[8])))
      }, error = function(e) {
        failures <<- c(failures, paste(root, conditionMessage(e), sep = " : "))
        NULL
      })
      if (!is.null(result)) {
        files <- list.files(src, recursive = TRUE, full.names = FALSE, all.files = FALSE)
        files <- files[!startsWith(files, "00_EXTRACTED_RAW_FILES/")]
        for (f in files) {
          dest <- file.path(sdir, f); mk(dirname(dest))
          check(file.copy(file.path(src, f), dest, overwrite = FALSE), "Could not copy cache output: ", f)
        }
        write_csv(data.frame(source = file.path(src, rel),
          md5 = unname(tools::md5sum(file.path(src, rel)))), file.path(adir, "MICROARRAY_CACHE_PROVENANCE.csv"))
        session_candidates <- file.path(root, c("00_METADATA/CACHED_MICROARRAY_SESSION_INFO.txt",
                                                "00_METADATA/SESSION_INFO.txt"))
        prior_session <- session_candidates[file.exists(session_candidates)][1L]
        if (!is.na(prior_session)) file.copy(prior_session, file.path(adir, "CACHED_MICROARRAY_SESSION_INFO.txt"))
        writeLines(c("Microarray outputs reused from:", root,
          "Checks: fixed signature, sample map, paired design, alias ownership, DE/rank agreement, coverage and ES.",
          "This does not constitute refitting the model. The new final BH adjustment is recalculated."),
          file.path(sdir, "CACHE_REUSE.txt"))
        message("Reusing completed corrected GSE198692 output: ", root)
        return(result)
      }
    }
    if (length(failures)) writeLines(failures, file.path(adir, "MICROARRAY_CACHE_NOT_USED.txt"))
    NULL
  }

  withCallingHandlers({

    message("1/3: GSE76987 -- 8 sporadic SSL specimens versus 10 control-right specimens")
    rq <- mk(rdir, "00_QC_AND_METADATA"); rd <- mk(rdir, "01_DESEQ2_TABLES")
    rg <- mk(rdir, "02_SSL24_GSEA"); rf <- mk(rdir, "03_FIGURES")
    re <- mk(rdir, "04_EXPRESSION")
    soft_path <- obtain_input("GSE76987", "GSE76987_family.soft.gz", "soft", "soft")
    record_input(soft_path)
    all_samples <- parse_geo_samples(soft_path)
    write_csv(all_samples, file.path(rq, "SOURCE_COHORT_BEFORE_SINGLE_LESION_RULE.csv"))
    all_samples <- restrict_rnaseq_samples(all_samples)
    write_csv(all_samples, file.path(rq, "ALL_GEO_SAMPLE_SELECTION_AUDIT.csv"))
    write_csv(all_samples[all_samples$eligible_before_single_lesion_rule, , drop = FALSE],
              file.path(rq, "SINGLE_LESION_SELECTION_AUDIT.csv"))
    writeLines(restriction_note, file.path(rq, "PATIENT_LINKAGE_EVIDENCE_AND_ASSUMPTION.txt"))
    message("Panel A restriction: retaining 118SSP_C1; excluding 118SSP_C2.")
    message("Patient linkage is inferred from author labels plus Supplementary Table 1a; see the evidence note.")
    workbook <- obtain_input("GSE76987", "GSE76987_RightColonProcessed.xlsx", "xlsx")
    record_input(workbook)
    imported <- read_geo_workbook(workbook, all_samples, rq)
    counts <- imported$counts; rn_meta <- imported$metadata

    check(ncol(counts) == 18L && sum(rn_meta$group == "SSL") == 8L &&
            sum(rn_meta$group == "Control_Right") == 10L,
          "The imported analysis matrix must contain 8 SSL and 10 control specimens.")
    if (any(!is.na(rn_meta$patient_id))) {
      known <- !is.na(rn_meta$patient_id) & nzchar(rn_meta$patient_id)
      check(!anyDuplicated(rn_meta$patient_id[known]),
            "GEO explicitly identifies repeated patients in the selected subset.\n",
            "A group-only independent-sample model is not fitted without resolving this.")
    }
    check(identical(colnames(counts), rownames(rn_meta)), "RNA-seq sample alignment failed.")
    matrix_csv(counts, file.path(re, "RAW_COUNTS_INCLUDED_SPECIMENS.csv"), id = "gene_id")
    gene_ids <- rownames(counts)
    db <- getExportedValue("org.Hs.eg.db", "org.Hs.eg.db")
    available <- AnnotationDbi::keys(db, keytype = "ENSEMBL")
    query <- intersect(gene_ids, available)
    check(length(query) > 1000L, "Too few valid Ensembl identifiers in the downloaded gene table.")
    map <- AnnotationDbi::select(db, keys = query, keytype = "ENSEMBL", columns = c("SYMBOL", "ENTREZID"))
    map <- unique(as.data.frame(map))
    map$gene_symbol <- clean_symbol(map$SYMBOL)
    for (i in seq_len(nrow(aliases))) {
      named <- map$gene_symbol %in% c(aliases$old[i], aliases$new[i])
      check(all(!is.na(map$ENTREZID[named]) & map$ENTREZID[named] == aliases$entrez[i]),
            "Unexpected gene-ID mapping for ", aliases$new[i])
      map$gene_symbol[named] <- aliases$new[i]
    }
    write_csv(map, file.path(rq, "ENSEMBL_SYMBOL_MAPPING_ALL_RECORDS.csv"))
    groups <- split(map$gene_symbol, map$ENSEMBL)
    symbol_map <- vapply(groups, function(v) {
      v <- unique(v[!is.na(v) & nzchar(v)])
      if (length(v) == 1L) v else NA_character_
    }, character(1))
    gm <- data.frame(gene_id = gene_ids, gene_symbol = unname(symbol_map[gene_ids]))

    gm$mapping_status <- ifelse(is.na(gm$gene_symbol), "Unmapped_or_ambiguous", "Unique_symbol")
    write_csv(gm, file.path(rq, "GENE_ID_TO_SYMBOL.csv"))
    multi_sig <- table(gm$gene_symbol[gm$gene_symbol %in% signature])
    check(!any(multi_sig > 1L), "A signature gene maps to multiple counted Ensembl rows.\n",
          "See the mapping audit; the script does not sum potentially overlapping gene models.")
    allzero <- rowSums(counts) == 0
    write_csv(data.frame(gm, total_selected_count = rowSums(counts),
      n_samples_detected = rowSums(counts > 0), retained_for_DESeq2 = !allzero),
      file.path(rq, "COUNT_FILTER_AUDIT.csv"))
    write_csv(data.frame(sample_id = colnames(counts), group = rn_meta$group,
      total_deposited_gene_counts = colSums(counts), detected_genes = colSums(counts > 0)),
      file.path(rq, "SAMPLE_LIBRARY_QC.csv"))

    dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts[!allzero, , drop = FALSE],
      colData = rn_meta[, "group", drop = FALSE], design = ~ group)
    dds <- DESeq2::DESeq(dds, test = "Wald", fitType = "parametric", quiet = FALSE)
    res <- DESeq2::results(dds, contrast = c("group", "SSL", "Control_Right"),
      alpha = 0.05, independentFiltering = TRUE, cooksCutoff = TRUE)
    rn_de <- data.frame(gene_id = rownames(res), as.data.frame(res), row.names = NULL,
                         check.names = FALSE)
    rn_de$gene_symbol <- gm$gene_symbol[match(rn_de$gene_id, gm$gene_id)]
    rn_de$comparison <- "Sporadic SSL versus control right-colon mucosa"
    write_csv(rn_de, file.path(rd, "SSL_VS_CONTROL_RIGHT_ALL_GENES.csv"))
    saveRDS(dds, file.path(re, "DESEQ2_MODEL.rds"))
    write_csv(data.frame(sample_id = rn_meta$sample_id,
      stats::model.matrix(~ group, rn_meta), check.names = FALSE),
      file.path(rq, "DESEQ2_DESIGN_MATRIX.csv"))
    nr <- DESeq2::counts(dds, normalized = TRUE)
    matrix_csv(nr, file.path(re, "NORMALIZED_COUNTS.csv"), id = "gene_id")
    primary <- run_gsea(make_ranks(rn_de), "Sporadic SSL versus control right-colon mucosa",
                         76987L, rg, "SSL_VS_CONTROL_RIGHT")
    rn_sig <- signature_table(rn_de, primary$ranks)
    rn_sig$n_SSL_samples_with_counts <- vapply(signature, function(s) {
      j <- which(gm$gene_symbol %in% s)
      if (!length(j)) return(NA_integer_)
      as.integer(sum(colSums(counts[j, rn_meta$group == "SSL", drop = FALSE]) > 0))
    }, integer(1))
    rn_sig$n_control_samples_with_counts <- vapply(signature, function(s) {
      j <- which(gm$gene_symbol %in% s)
      if (!length(j)) return(NA_integer_)
      as.integer(sum(colSums(counts[j, rn_meta$group == "Control_Right", drop = FALSE]) > 0))
    }, integer(1))
    rn_sig$leading_edge <- rn_sig$gene_symbol %in% strsplit(primary$result$leading_edge_genes, ";", fixed = TRUE)[[1L]]
    write_csv(rn_sig, file.path(rd, "SSL24_GENE_RESULTS_AND_COUNT_AUDIT.csv"))
    save_plot(gsea_plot(primary, "GSE76987: sporadic SSL versus control right colon"), file.path(rf, "SSL24_GSEA"))
    volcano(rn_de, file.path(rf, "DESEQ2_VOLCANO"), "GSE76987: sporadic SSL versus control right colon")

    vsd <- DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
    vx <- SummarizedExperiment::assay(vsd)
    vj <- head(order(apply(vx, 1, stats::var), decreasing = TRUE), min(2000L, nrow(vx)))
    rpca <- stats::prcomp(t(vx[vj, , drop = FALSE]), scale. = FALSE)
    save_pca(rpca$x[, 1:2, drop = FALSE], rn_meta, file.path(rf, "PCA"), "GSE76987 selected samples",
              100 * rpca$sdev^2 / sum(rpca$sdev^2))
    rann <- data.frame(Group = rn_meta$group, row.names = rn_meta$sample_id)
    heat(as.matrix(stats::dist(t(vx[vj, , drop = FALSE]))), file.path(rf, "SAMPLE_DISTANCE"),
         "GSE76987 sample distances", rann, cluster_cols = TRUE, w = 9, h = 8)
    sj <- match(rn_sig$gene_id, rownames(vx)); sk <- !is.na(sj)
    if (sum(sk) >= 2L) {
      sx <- vx[sj[sk], , drop = FALSE]; rownames(sx) <- rn_sig$gene_symbol[sk]
      matrix_csv(sx, file.path(re, "SSL24_VST_EXPRESSION.csv"))
      heat(sx, file.path(rf, "SSL24_HEATMAP"), "GSE76987: fixed SSL24 signature", rann,
           row_z = TRUE, w = 11, h = 9)
    }
    writeLines(c("GSE76987: eight-specimen analysis with documented patient-linkage inference",
      restriction_note, "",
      "Analysis: 8 included sporadic SSL specimens versus 10 control right-colon specimens.",
      "The specimen-to-patient cross-reference is inferred, not directly supplied in the supplement.",
      "SSA/P is referred to as SSL without claiming a new histology review.",
      "SPS, uninvolved mucosa, HP, adenomas, cancers and control left-colon specimens are excluded.",
      "Unpaired DESeq2 with tissue group only; no patient-specific covariate is invented.",
      "Lesions were drawn from both right and left colon; controls are right-colon specimens.",
      "Main-paper and supplemental side counts disagree; no exact side counts are assigned here.",
      "RightColonProcessed is an ANALYSIS workbook name, not proof of each lesion's location.",
      "The contrast is not fully site-matched or site-adjusted.",
      "All available Analyzed Genes rows are imported; upstream source exclusions remain.",
      "Only all-zero rows across the 18 included specimens are removed before DESeq2.",
      "Ranks require finite Wald statistics and finite nominal P values, not a significance cutoff.",
      "No log2CPM >= 2 filter is applied to the signature or the ranked background.",
      "Reference: Kanth et al. (2016), doi:10.1158/1940-6207.CAPR-15-0363; Suppl Table 1a; GEO GSE76987."),
      file.path(rdir, "DESIGN_AND_LIMITATIONS.txt"))

    run_microarray <- function() {

    tar_path <- obtain_input("GSE198692", "GSE198692_RAW.tar", "tar")
    record_input(tar_path)
    message("2/3: GSE198692 -- identifier-corrected paired analysis from raw arrays")
    sr <- mk(sdir, "00_EXTRACTED_RAW_FILES"); sq <- mk(sdir, "01_QC_AND_METADATA")
    sde <- mk(sdir, "02_PAIRED_LIMMA_DE"); sg <- mk(sdir, "03_SSL24_GSEA")
    ss <- mk(sdir, "04_SSL24_SAMPLE_SCORES"); sf <- mk(sdir, "05_FIGURES")
    sl <- mk(sdir, "06_LEAVE_ONE_PAIR_OUT")
    pairs <- c("01", "03", "04", "08", "10", "12", "13", "14", "17", "21")
    meta <- data.frame(GSM = paste0("GSM", 5955144:5955163),
      Sample = paste0(rep(pairs, each = 2), rep(c("N", "T"), 10)),
      PairID = factor(rep(pairs, each = 2), levels = pairs),
      Group = factor(rep(c("Adjacent_Normal", "SSL"), 10), levels = c("Adjacent_Normal", "SSL")))
    check(all(table(meta$PairID, meta$Group) == 1L), "SSL pairs are incomplete or duplicated.")
    members <- utils::untar(tar_path, list = TRUE)
    check(length(members) > 0 && !any(grepl("(^[/\\\\])|(^[A-Za-z]:)|(^|[/\\\\])\\.\\.([/\\\\]|$)", members)),
          "Unsafe or empty paths in the Agilent archive.")
    utils::untar(tar_path, exdir = sr)
    gz <- list.files(sr, pattern = "\\.txt\\.gz$", recursive = TRUE, full.names = TRUE)
    for (f in gz) R.utils::gunzip(f, destname = sub("\\.gz$", "", f), remove = FALSE, overwrite = FALSE)
    txt <- list.files(sr, pattern = "\\.txt$", recursive = TRUE, full.names = TRUE)
    fgsm <- vapply(basename(txt), function(x) {
      m <- regexpr("GSM[0-9]+", x); if (m[1] > 0) regmatches(x, m) else ""
    }, character(1))
    use <- fgsm %in% meta$GSM
    check(sum(use) == 20L && !anyDuplicated(fgsm[use]), "Expected one raw array for each of 20 GSMs.")
    meta$File <- txt[match(meta$GSM, fgsm)]
    write_csv(meta, file.path(sq, "SAMPLE_AND_RAW_FILE_MAP.csv"))
    lines <- readLines(meta$File[1], n = 100L, warn = FALSE)
    features <- lines[startsWith(lines, "FEATURES\t")]
    check(length(features) == 1L, "Agilent FEATURES header missing.")
    fc <- strsplit(features, "\t", fixed = TRUE)[[1]][-1]
    flags <- c("gIsWellAboveBG", "gIsSaturated", "gIsFeatNonUnifOL", "gIsFeatPopnOL")
    check(all(c("ProbeName", "ControlType", "gMedianSignal", "gBGMedianSignal", flags) %in% fc),
          "Raw Agilent columns needed by the original workflow are missing.")
    raw <- limma::read.maimages(meta$File, source = "agilent.median", green.only = TRUE,
      annotation = intersect(c("ProbeName", "ControlType", "GeneName", "SystematicName"), fc),
      other.columns = flags)
    colnames(raw$E) <- meta$Sample
    bc <- limma::backgroundCorrect(raw, method = "normexp", offset = 50)
    norm <- limma::normalizeBetweenArrays(bc, method = "quantile")
    colnames(norm$E) <- meta$Sample
    boxplot_intensity <- function(m, title, stem) draw_base(stem, function() {
      graphics::par(mar = c(8, 4, 3, 1)); graphics::boxplot(m, names = meta$Sample,
        las = 2, outline = FALSE, main = title, ylab = "Log2 intensity")
    }, 12, 6)
    boxplot_intensity(log2(pmax(raw$E, 1)), "Raw Agilent foreground", file.path(sf, "RAW_INTENSITIES"))
    boxplot_intensity(norm$E, "Normalized Agilent expression", file.path(sf, "NORMALIZED_INTENSITIES"))
    flag <- function(nm) {
      x <- raw$other[[nm]]; check(!is.null(x), "Missing flag: ", nm)
      y <- toupper(trimws(as.character(x)))
      check(all(!is.na(y) & y %in% c("0", "1", "TRUE", "FALSE", "T", "F", "YES", "NO", "Y", "N")),
            "Unrecognised flag values in ", nm)
      matrix(y %in% c("1", "TRUE", "T", "YES", "Y"), nrow(x), ncol(x))
    }
    good <- flag("gIsWellAboveBG") & !flag("gIsSaturated") &
      !flag("gIsFeatNonUnifOL") & !flag("gIsFeatPopnOL")
    nN <- rowSums(good[, meta$Group == "Adjacent_Normal", drop = FALSE])
    nS <- rowSums(good[, meta$Group == "SSL", drop = FALSE])
    control <- toupper(trimws(as.character(raw$genes$ControlType)))
    noncontrol <- is.na(control) | control %in% c("0", "FALSE", "F", "NO", "N")

    local_candidates <- file.path(c(input_dir, download_dir), "GPL21185.rds")
    local_gpl <- local_candidates[file.exists(local_candidates)][1L]
    if (is.na(local_gpl)) local_gpl <- file.path(download_dir, "GPL21185.rds")
    gpl <- if (file.exists(local_gpl)) as.data.frame(readRDS(local_gpl)) else {
      message("Obtaining GPL21185 annotation from GEO ...")
      as.data.frame(GEOquery::Table(GEOquery::getGEO("GPL21185", destdir = sq, AnnotGPL = FALSE)))
    }
    need_cols(gpl, c("ID", "GENE_SYMBOL", "LOCUSLINK_ID"), "GPL21185")
    saveRDS(gpl, file.path(sq, "GPL21185.rds"))
    if (!file.exists(file.path(download_dir, "GPL21185.rds")))
      saveRDS(gpl, file.path(download_dir, "GPL21185.rds"))
    record_input(file.path(sq, "GPL21185.rds"))

    ga <- unique(gpl[, c("ID", "GENE_SYMBOL", "LOCUSLINK_ID"), drop = FALSE])
    check(!anyDuplicated(ga$ID), "GPL21185 has conflicting annotation rows for one probe.")
    gi <- match(raw$genes$ProbeName, ga$ID)
    from_gpl <- clean_symbol(ga$GENE_SYMBOL[gi])
    from_raw <- if ("GeneName" %in% names(raw$genes)) clean_symbol(raw$genes$GeneName) else rep(NA_character_, nrow(raw$E))
    a <- data.frame(ProbeID = as.character(raw$genes$ProbeName),
      GPL_GeneSymbol = as.character(ga$GENE_SYMBOL[gi]), GPL_Entrez = as.character(ga$LOCUSLINK_ID[gi]),
      GeneSymbol = ifelse(!is.na(from_gpl), from_gpl, from_raw),
      passes_noncontrol = noncontrol, n_good_normal = nN, n_good_SSL = nS,
      passes_detection = nN >= 3L | nS >= 3L)
    a <- harmonise_probes(a)
    a$retained <- a$passes_noncontrol & a$passes_detection &
      !is.na(a$GeneSymbol) & nzchar(a$GeneSymbol)
    write_csv(a, file.path(sq, "PROBE_ANNOTATION_AND_FILTERING_AUDIT.csv"))
    check(sum(a$retained) >= 1000L, "Too few retained probes; inspect the audit.")
    expr <- as.matrix(limma::avereps(norm$E[a$retained, , drop = FALSE], ID = a$GeneSymbol[a$retained]))
    check(all(is.finite(expr)) && !anyDuplicated(rownames(expr)), "Invalid gene-level expression.")
    check(identical(colnames(expr), meta$Sample), "Expression/sample alignment failure.")
    probe_n <- table(a$GeneSymbol[a$retained])
    matrix_csv(expr, file.path(sq, "GENE_LEVEL_NORMALIZED_EXPRESSION.csv"))
    write_csv(data.frame(gene_symbol = names(probe_n), n_retained_probes = as.integer(probe_n)),
              file.path(sq, "RETAINED_PROBES_PER_GENE.csv"))

    paired_fit <- function(e, md) {
      md <- droplevels(md); design <- stats::model.matrix(~ PairID + Group, data = md)
      check(qr(design)$rank == ncol(design) && "GroupSSL" %in% colnames(design), "Invalid paired design.")
      fit <- limma::eBayes(limma::lmFit(e, design), trend = TRUE, robust = TRUE)
      t <- limma::topTable(fit, coef = "GroupSSL", number = Inf, sort.by = "P", adjust.method = "BH")
      d <- data.frame(gene_id = rownames(t), gene_symbol = rownames(t),
        log2FoldChange = t$logFC, baseMean = t$AveExpr, stat = t$t, pvalue = t$P.Value,
        padj = t$adj.P.Val, log_odds_DE = t$B, row.names = NULL)
      list(de = d, design = design, fit = fit)
    }
    fitted <- paired_fit(expr, meta); de <- fitted$de
    de$n_retained_probes <- as.integer(probe_n[de$gene_symbol])
    saveRDS(fitted$fit, file.path(sde, "PAIRED_LIMMA_MODEL.rds"))
    write_csv(data.frame(Sample = meta$Sample, fitted$design, check.names = FALSE),
              file.path(sq, "PAIRED_LIMMA_DESIGN_MATRIX.csv"))
    write_csv(de, file.path(sde, "SSL_VS_PAIRED_NORMAL_ALL_GENES.csv"))
    sig <- is.finite(de$padj) & de$padj < 0.05
    write_csv(data.frame(genes_tested = nrow(de), FDR_below_0_05 = sum(sig),
      higher_in_SSL = sum(sig & de$log2FoldChange > 0), lower_in_SSL = sum(sig & de$log2FoldChange < 0)),
      file.path(sde, "PAIRED_LIMMA_SUMMARY.csv"))
    sslg <- run_gsea(make_ranks(de), "SSL versus paired adjacent normal mucosa", 198692L, sg, "SSL_VS_PAIRED_NORMAL")
    write_csv(signature_table(de, sslg$ranks), file.path(sde, "SSL24_GENE_RESULTS.csv"))
    save_plot(gsea_plot(sslg, "SSL versus paired adjacent normal mucosa"), file.path(sf, "SSL24_GSEA"))
    volcano(de, file.path(sf, "PAIRED_LIMMA_VOLCANO"), "SSL versus paired adjacent normal mucosa")

    variable <- names(sort(apply(expr, 1, stats::var), decreasing = TRUE))
    variable <- head(variable, min(2000L, length(variable)))
    pca <- stats::prcomp(t(expr[variable, , drop = FALSE]), scale. = FALSE)
    smeta <- data.frame(sample_id = meta$Sample, group = meta$Group)
    save_pca(pca$x[, 1:2], smeta, file.path(sf, "PCA"), "GSE198692 paired samples",
              100 * pca$sdev^2 / sum(pca$sdev^2))
    sann <- data.frame(Group = meta$Group, PairID = meta$PairID, row.names = meta$Sample)
    heat(as.matrix(stats::dist(t(expr[variable, , drop = FALSE]))), file.path(sf, "SAMPLE_DISTANCE"),
         "GSE198692 sample distances", sann, cluster_cols = TRUE, w = 9, h = 8)
    measured <- intersect(signature, rownames(expr))
    check(length(measured) >= 5L, "Insufficient measured genes for the original sample-score analysis.")
    scores <- as.data.frame(singscore::simpleScore(singscore::rankGenes(expr),
      upSet = measured, centerScore = TRUE, knownDirection = TRUE))
    scores$Sample <- rownames(scores)
    score_order <- match(meta$Sample, scores$Sample)
    check(!anyNA(score_order), "Singscore sample names do not match metadata.")
    scores <- cbind(meta[, c("Sample", "PairID", "Group")],
                    scores[score_order, setdiff(names(scores), "Sample"), drop = FALSE])
    rownames(scores) <- NULL
    normal_id <- paste0(pairs, "N"); ssl_id <- paste0(pairs, "T")
    sw <- data.frame(PairID = pairs,
      Adjacent_Normal = scores$TotalScore[match(normal_id, scores$Sample)],
      SSL = scores$TotalScore[match(ssl_id, scores$Sample)])
    sw$difference <- sw$SSL - sw$Adjacent_Normal
    wt <- function(w, side) stats::wilcox.test(w$SSL, w$Adjacent_Normal,
      paired = TRUE, alternative = side, exact = FALSE, correct = TRUE)
    score_stats <- data.frame(n_pairs = nrow(sw), higher_in_SSL = sum(sw$difference > 0),
      median_paired_difference = median(sw$difference), mean_paired_difference = mean(sw$difference),
      two_sided_P = wt(sw, "two.sided")$p.value, one_sided_greater_P = wt(sw, "greater")$p.value)
    write_csv(scores, file.path(ss, "SINGSCORES_ALL_SAMPLES.csv"))
    write_csv(sw, file.path(ss, "PAIRED_SCORE_DIFFERENCES.csv"))
    write_csv(score_stats, file.path(ss, "PAIRED_SCORE_TESTS.csv"))
    save_plot(ggplot(scores, aes(Group, TotalScore, group = PairID)) +
      geom_line(alpha = 0.6, linewidth = 0.5) + geom_point(size = 2.5) + theme_classic() +
      scale_x_discrete(labels = c(Adjacent_Normal = "Paired normal", SSL = "SSL")) +
      labs(title = "Serrated Lesion Gastric Metaplasia Gene Signature", x = NULL, y = "Singscore",
           subtitle = paste0("Two-sided paired Wilcoxon P = ", fmt(score_stats$two_sided_P))),
      file.path(sf, "SSL24_PAIRED_SINGSCORE"), 7, 6)
    se <- expr[measured, meta$Sample, drop = FALSE]
    heat(se, file.path(sf, "SSL24_ZSCORE_HEATMAP"), "SSL24 sample expression", sann,
         row_z = TRUE, w = 12, h = 9, breaks = seq(-2.5, 2.5, length.out = 101))
    diff <- se[, ssl_id, drop = FALSE] - se[, normal_id, drop = FALSE]
    colnames(diff) <- pairs
    matrix_csv(diff, file.path(ss, "SSL24_PAIRED_GENE_LOG2_DIFFERENCES.csv"))
    lim <- max(abs(diff)); if (lim == 0) lim <- 1
    heat(diff, file.path(sf, "SSL24_PAIRED_LOG2FC_HEATMAP"), "SSL24 paired SSL-minus-normal log2 differences",
         w = 9, h = 9, breaks = seq(-lim, lim, length.out = 101))
    long <- data.frame(gene_symbol = rep(rownames(se), ncol(se)),
      Sample = rep(colnames(se), each = nrow(se)), Expression = as.vector(se))
    long$PairID <- meta$PairID[match(long$Sample, meta$Sample)]
    long$Group <- meta$Group[match(long$Sample, meta$Sample)]
    save_plot(ggplot(long, aes(Group, Expression, group = PairID)) +
      geom_line(linewidth = 0.3, alpha = 0.5) + geom_point(size = 1.2) +
      facet_wrap(~ gene_symbol, scales = "free_y", ncol = 4) + theme_classic(base_size = 9) +
      scale_x_discrete(labels = c(Adjacent_Normal = "Normal", SSL = "SSL")) +
      labs(x = NULL, y = "Normalized log2 expression", title = "Fixed SSL24 genes"),
      file.path(sf, "SSL24_GENES_PAIRED_EXPRESSION"), 12, 16)

    lopo <- list()
    for (pair in pairs) {
      kp <- as.character(meta$PairID) != pair
      df_i <- paired_fit(expr[, kp, drop = FALSE], meta[kp, , drop = FALSE])$de
      gi <- run_gsea(make_ranks(df_i), paste("SSL versus paired normal; omitted pair", pair),
                      198692L + as.integer(pair), sl, paste0("OMIT_PAIR_", pair))
      wi <- sw[sw$PairID != pair, , drop = FALSE]
      lopo[[pair]] <- data.frame(excluded_pair = pair, gi$result,
        score_pairs_higher_in_SSL = sum(wi$difference > 0),
        score_median_paired_difference = median(wi$difference),
        score_one_sided_P = wt(wi, "greater")$p.value)
    }
    lopo <- do.call(rbind, lopo); rownames(lopo) <- NULL
    lopo$BH_across_10_LOPO_runs <- stats::p.adjust(lopo$pval, "BH", n = 10)
    write_csv(lopo, file.path(sl, "LEAVE_ONE_PAIR_OUT_SUMMARY.csv"))
    save_plot(ggplot(lopo, aes(excluded_pair, NES)) + geom_point(size = 3) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) + theme_classic() +
      labs(x = "Omitted matched pair", y = "Normalized enrichment score", title = "SSL24 leave-one-pair-out analysis"),
      file.path(sf, "LEAVE_ONE_PAIR_OUT_NES"), 8, 5.5)

      list(gsea = sslg, metadata = meta, score_tests = score_stats, lopo = lopo)
    }
    message("2/3: checking for an already completed corrected microarray run")
    microarray <- reuse_microarray()
    if (is.null(microarray)) microarray <- run_microarray()
    sslg <- microarray$gsea; meta <- microarray$metadata
    score_stats <- microarray$score_tests; lopo <- microarray$lopo

    message("3/3: final Figure 3.8.2 and reporting tables")
    final <- rbind(data.frame(cohort = "GSE76987", primary$result),
                     data.frame(cohort = "GSE198692", sslg$result))
    final$BH_FDR_across_external_cohorts <- stats::p.adjust(final$pval, method = "BH", n = 2L)
    final$ranking_statistic <- c("DESeq2 Wald statistic", "Paired limma moderated t statistic")
    final$model <- c("~ group; 8 sporadic SSL specimens versus 10 control-right-colon specimens; C2 omitted under inferred 118SSP linkage, unpaired",
                     "~ PairID + Group; 10 matched pairs")
    final$n_lesion_specimens <- c(sum(rn_meta$group == "SSL"), 10L)
    final$n_control_specimens <- c(sum(rn_meta$group == "Control_Right"), 10L)
    final$patient_linkage_status <- c(linkage_status, "Original ten matched pairs retained")
    final$excluded_repeated_lesion_GSM <- c(excluded_lesion, NA_character_)
    final$excluded_repeated_lesion_author_label <- c("118SSP_C2", NA_character_)
    final$FDR_family <- "Two displayed SSL-cohort GSEA tests; prior UACL exploration archived separately"
    write_csv(final, file.path(fdir, "FINAL_SSL24_EXTERNAL_VALIDATION_GSEA_SUMMARY.csv"))

    curves <- list(enrichment_curve(primary$ranks), enrichment_curve(sslg$ranks))
    yy <- unlist(lapply(curves, function(x) x$ES))
    yy <- yy[is.finite(yy)]; if (!length(yy)) yy <- c(0, 1)
    yy[abs(yy) < 1e-12] <- 0
    lo <- floor(min(c(0, yy)) * 4) / 4; hi <- ceiling(max(c(0, yy)) * 4) / 4
    if (lo == hi) hi <- lo + 0.25
    span <- hi - lo; ylim <- c(lo - 0.06 * span, hi + 0.10 * span)
    titles <- c("A  GSE76987", "B  GSE198692")
    subs <- c("Bulk RNA-seq: 8 sporadic SSL specimens\nversus 10 control-right-colon specimens",
              "Independent SSL cohort\nSSL versus paired adjacent normal mucosa")
    cores <- lapply(1:2, function(i) {
      cv <- curves[[i]]; rr <- final[i, ]
      note <- if (rr$analysis_status == "RUN") paste0("NES = ", sprintf("%.2f", rr$NES),
        "\nBH FDR = ", fmt(rr$BH_FDR_across_external_cohorts), "\nGenes = ", rr$size, "/24") else
        paste0("GSEA not evaluable\nGenes = ", rr$size, "/24")
      if (!nrow(cv)) cv <- data.frame(position = c(0, 100), ES = NA_real_, hit = FALSE)
      ggplot(cv, aes(position, ES)) +
        geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
        geom_line(linewidth = 0.7, na.rm = TRUE) +
        geom_rug(data = cv[cv$hit, , drop = FALSE], aes(x = position),
                 inherit.aes = FALSE, sides = "b", linewidth = 0.18) +
        annotate("text", x = 98, y = ylim[2] - 0.025 * diff(ylim), label = note,
                 hjust = 1, vjust = 1, size = 2.4, lineheight = 1.05) +
        scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25), expand = c(0.01, 0)) +
        scale_y_continuous(limits = ylim, breaks = seq(lo, hi, 0.25), expand = c(0, 0)) +
        theme_classic(base_size = 8) + theme(plot.margin = margin(2, 4, 4, 4, "mm")) +
        labs(x = "Ranked gene list position (%)", y = "Running enrichment score")
    })
    cores <- cowplot::align_plots(plotlist = cores, align = "hv", axis = "tblr")
    panels <- lapply(1:2, function(i) {
      header <- cowplot::ggdraw() +
        cowplot::draw_label(titles[i], x = 0.035, y = 0.96, hjust = 0, vjust = 1, fontface = "bold", size = 9) +
        cowplot::draw_label(subs[i], x = 0.035, y = 0.49, hjust = 0, vjust = 0.5, size = 7.5, lineheight = 1)
      cowplot::plot_grid(header, cores[[i]], ncol = 1, rel_heights = c(0.21, 0.79))
    })
    figure <- cowplot::plot_grid(panels[[1]], NULL, panels[[2]], nrow = 1, rel_widths = c(1, 0.035, 1))
    save_plot(figure, file.path(fdir, "Figure_3_8_2"), 183 / 25.4, 85 / 25.4)
    ggplot2::ggsave(file.path(fdir, "Figure_3_8_2.tif"), figure, width = 183, height = 85,
                    units = "mm", dpi = 600, compression = "lzw", bg = "white")
    for (i in 1:2) save_plot(panels[[i]], file.path(fdir, paste0("Panel_", LETTERS[i])), 90/25.4, 85/25.4)
    for (i in 1:2) write_csv(curves[[i]], file.path(fdir, paste0(final$cohort[i], "_PLOTTED_CURVE.csv")))
    notes <- c("SSL24 EXTERNAL VALIDATION -- COMPLETED EXECUTION", paste("Output:", out),
      "", capture.output(print(final, row.names = FALSE)), "",
      "Panel A: GSE76987 bulk RNA sequencing, 8 sporadic SSL specimens versus 10 control-right specimens.",
      restriction_note,
      "Patient cohort and tissue metadata determine inclusion before any enrichment result is calculated.",
      "SPS specimens with documented repeated lesions are excluded. No random sample selection.",
      "Lesions include both colonic sides whereas controls are right-sided; this is not fully site-matched.",
      "Main text and supplement disagree on exact side counts, so no exact count is asserted here.",
      "RightColonProcessed.xlsx is an analysis-workbook name, not a lesion-location label.",
      "Use all available raw Counts rows from Analyzed Genes, not FPKM/RLog or DE-filtered sheets.",
      "The original deposited workbook already reflects upstream processing; missing source genes stay missing.",
      "Additional RNA-seq prefilter: only genes with zero counts across the 18 included specimens are removed.",
      "DESeq2 rankings require finite Wald statistics AND finite nominal P values, not finite padj.",
      "No log2CPM cutoff and no individual-gene significance cutoff select the GSEA universe.",
      "Panel B: GSE198692, 10 matched SSL-normal pairs, corrected identifiers before probe averaging.",
      "Agilent normexp offset 50 + quantile normalization; good detection in >=3 arrays in either group.",
      "Paired limma design ~ PairID + Group, eBayes trend=TRUE and robust=TRUE.",
      "See CACHE_REUSE.txt when previously completed microarray results were reused rather than refitted.",
      "The existing microarray singscore, paired tests and ten LOPO analyses are retained.",
      "LOPO keeps full-cohort preprocessing fixed and refits the paired model; its BH family is separate.",
      "Displayed BH family: two current SSL-cohort GSEA tests. Do not call the historical exploration prespecified.",
      "The archived UACL output is not modified or used to choose a significant subset of the new dataset.",
      "A nonsignificant or negative result is exported unchanged; no automatic result-driven sample switching.",
      "Coverage is actual ranked-list overlap, not proof that every signature gene is differentially expressed.",
      "Use the NEW numerical summary and Figure_3_8_2; old BH values used a different second cohort.",
      "References: Kanth et al. (2016), doi:10.1158/1940-6207.CAPR-15-0363; Ohki et al. (2022), doi:10.1002/cam4.4810.")
    writeLines(notes, file.path(out, "READ_RESULTS_FIRST.txt"))
    writeLines(capture.output(sessionInfo()), file.path(adir, "SESSION_INFO.txt"))
    writeLines(c(
      "Figure 3.8.2 | External evaluation of the Serrated Lesion Gastric Metaplasia Gene Signature.",
      "(A) GSE76987 bulk RNA sequencing: eight sporadic SSL specimens versus ten control right-colon specimens.",
      "(B) GSE198692 Agilent microarray: ten SSL specimens versus patient-matched surrounding normal mucosa.",
      "Genes were ranked by the DESeq2 Wald statistic (A) and paired limma moderated t statistic (B).",
      "Positive NES indicates enrichment in SSL. Vertical marks denote signature genes in each ranked list.",
      "Coverage is the overlap of the fixed 24-gene set with each complete valid ranked universe.",
      "Displayed P values are Benjamini-Hochberg adjusted across the two displayed cohort GSEA tests.",
      "Panel A excludes 118SSP_C2 under a fixed C1-before-C2 rule; same-patient linkage is inferred from author labels plus the supplement.",
      "The supplementary table confirms a repeated patient but does not explicitly map its rows to GEO accessions.",
      "Panel A is not patient-matched or fully site-matched; lesions came from both colonic sides.",
      "SSL, sessile serrated lesion; NES, normalised enrichment score; FDR, false discovery rate."),
      file.path(fdir, "FIGURE_CAPTION_AND_DESIGN_NOTE.txt"))
    saveRDS(list(signature = signature, rnaseq_design = rn_meta,
      rnaseq_selection_audit = all_samples, patient_linkage_note = restriction_note,
      patient_linkage_status = linkage_status, rnaseq_result = primary,
      rnaseq_signature_audit = rn_sig, ssl_design = meta, ssl_result = sslg,
      paired_score_tests = score_stats, lopo = lopo, final_summary = final),
      file.path(out, "ANALYSIS_RESULTS.rds"))
    required <- file.path(fdir, c("Figure_3_8_2.pdf", "Figure_3_8_2.png",
                                  "FINAL_SSL24_EXTERNAL_VALIDATION_GSEA_SUMMARY.csv"))
    check(all(file.exists(required)), "Required final outputs are missing.")
    writeLines(paste("Execution finished:", Sys.time()), file.path(out, "RUN_COMPLETED.txt"))
    completed <- TRUE
    print(final, row.names = FALSE)
    message("Finished. All outputs are in: ", out)
    message("Main figure: ", file.path(fdir, "Figure_3_8_2.png"))
    message("Use the new numerical summary and read DESIGN_AND_LIMITATIONS.txt for panel A.")
    invisible(final)
  }, warning = function(w) {
    cat(conditionMessage(w), "\n", file = log_path, append = TRUE)

  })
}

ssl24_external_validation(input_dir)
