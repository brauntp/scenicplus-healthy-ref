#!/usr/bin/env Rscript
# Inspect an ArchRProject for SCENIC+ readiness.
#
# READ-ONLY: loads the project with showLogo=FALSE and never calls saveArchRProject.
# Reports the peak set geometry, genome build, per-cell-type cell counts,
# fragment-file availability, and whether marker peaks / a peak matrix exist.
#
# Usage:
#   Rscript inspect_archr.R --proj /path/to/ArchRProject --out report_archr
#
# Emits <out>.json and <out>.md. Summary statistics only; no count data.

suppressPackageStartupMessages({
  library(ArchR)
})

# ---------------------------------------------------------------- args
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) && length(args) > i[1]) args[i[1] + 1] else default
}
# Accept both spellings: --proj here, --archr-project in
# 01_cistopic/export_from_archr.R. Same thing; taking either avoids a
# needless failure when copying a path between the two stages.
proj_path <- getarg("--proj")
if (is.null(proj_path)) proj_path <- getarg("--archr-project")
out_stem  <- getarg("--out", "report_archr")
if (is.null(proj_path))
  stop("path to the ArchRProject is required: --proj (or --archr-project) /path/to/ArchRProject")

# json writer: prefer jsonlite, fall back to a minimal serializer so this
# script never fails on a cluster module that lacks it.
have_jsonlite <- requireNamespace("jsonlite", quietly = TRUE)
rep <- list()
md  <- c("# ArchRProject inspection (SCENIC+ readiness)", "")

add_md <- function(...) md <<- c(md, ...)

# ---------------------------------------------------------------- load
message("Loading ArchRProject (read-only) ...")
proj <- loadArchRProject(path = proj_path, showLogo = FALSE)

rep$project_path <- proj_path
rep$archr_version <- as.character(utils::packageVersion("ArchR"))
rep$n_cells <- nrow(proj@cellColData)
rep$genome <- tryCatch(as.character(getGenome(proj)), error = function(e) NA_character_)
rep$sample_names <- as.character(unique(proj$Sample))
rep$n_samples <- length(rep$sample_names)

add_md(
  sprintf("- project: `%s`", proj_path),
  sprintf("- ArchR version: `%s`", rep$archr_version),
  sprintf("- genome: **%s**", rep$genome),
  sprintf("- cells: **%s**", format(rep$n_cells, big.mark = ",")),
  sprintf("- samples: **%d** (%s)", rep$n_samples,
          paste(utils::head(rep$sample_names, 8), collapse = ", ")),
  ""
)

# ---------------------------------------------------------------- cellColData
ccd <- as.data.frame(proj@cellColData)
rep$cellColData_columns <- colnames(ccd)

# Candidate grouping columns: low cardinality, complete, character/factor.
cand <- list()
for (cn in colnames(ccd)) {
  v <- ccd[[cn]]
  if (!(is.character(v) || is.factor(v) || is.logical(v))) next
  n <- length(unique(v[!is.na(v)]))
  if (n < 2 || n > 60) next
  if (any(is.na(v))) next
  tb <- sort(table(v), decreasing = TRUE)
  cand[[cn]] <- list(
    n_levels = n,
    min_cells = as.integer(min(tb)),
    levels = as.list(setNames(as.integer(tb), names(tb)))
  )
}
rep$candidate_group_keys <- cand

add_md("## Candidate grouping keys in cellColData", "")
if (!length(cand)) {
  add_md("**None found** (no complete, low-cardinality categorical column).", "")
} else {
  for (cn in names(cand)) {
    d <- cand[[cn]]
    add_md(sprintf("- **`%s`** -- %d levels, smallest has %d cells",
                   cn, d$n_levels, d$min_cells))
    small <- names(d$levels)[unlist(d$levels) < 50]
    if (length(small)) {
      add_md(sprintf("    - under 50 cells: %s",
                     paste(utils::head(small, 20), collapse = ", ")))
    }
  }
  add_md("")
}

# ---------------------------------------------------------------- peak set
add_md("## Peak set", "")
ps <- tryCatch(getPeakSet(proj), error = function(e) NULL)
if (is.null(ps) || !length(ps)) {
  rep$peakset <- list(present = FALSE)
  add_md("**No peak set in this project.** addReproduciblePeakSet() has not been run,",
         "or the peaks live in a different project object.", "")
} else {
  w <- BiocGenerics::width(ps)
  sn <- as.character(GenomeInfoDb::seqnames(ps))
  uniq_chr <- sort(unique(sn))
  std <- paste0("chr", c(1:22, "X", "Y"))
  rep$peakset <- list(
    present = TRUE,
    n_peaks = length(ps),
    width_median = as.numeric(stats::median(w)),
    width_min = as.integer(min(w)),
    width_max = as.integer(max(w)),
    width_n_distinct = length(unique(w)),
    fixed_width = length(unique(w)) == 1L,
    n_chroms = length(uniq_chr),
    has_chr_prefix = all(grepl("^chr", uniq_chr)),
    nonstandard_contigs = setdiff(uniq_chr, std),
    max_end_chr1 = if ("chr1" %in% sn) as.integer(max(BiocGenerics::end(ps)[sn == "chr1"])) else NA_integer_
  )
  # peakType / nearestGene are added by ArchR's peak annotation
  if (!is.null(S4Vectors::mcols(ps)) && "peakType" %in% colnames(S4Vectors::mcols(ps))) {
    tb <- table(S4Vectors::mcols(ps)$peakType)
    rep$peakset$peakType <- as.list(setNames(as.integer(tb), names(tb)))
  }
  rep$peakset$mcols_columns <- colnames(S4Vectors::mcols(ps))

  add_md(
    sprintf("- **%s peaks**", format(length(ps), big.mark = ",")),
    sprintf("- width: median **%.0f bp**, min %d, max %d, %d distinct -> fixed width: **%s**",
            rep$peakset$width_median, rep$peakset$width_min, rep$peakset$width_max,
            rep$peakset$width_n_distinct, rep$peakset$fixed_width),
    sprintf("- chroms: %d, `chr` prefix: **%s**",
            rep$peakset$n_chroms, rep$peakset$has_chr_prefix),
    sprintf("- non-standard contigs: %s",
            if (length(rep$peakset$nonstandard_contigs))
              paste(utils::head(rep$peakset$nonstandard_contigs, 20), collapse = ", ") else "none"),
    sprintf("- max end on chr1: %s (hg38 chr1 = 248,956,422; hg19 = 249,250,621)",
            rep$peakset$max_end_chr1),
    ""
  )
  if (!is.null(rep$peakset$peakType)) {
    add_md(sprintf("- peakType: %s", paste(
      sprintf("%s=%s", names(rep$peakset$peakType),
              format(unlist(rep$peakset$peakType), big.mark = ",")),
      collapse = ", ")), "")
  }
  # SCENIC+ needs a BED of the peak set; write one next to the report.
  bed <- data.frame(
    chrom = sn,
    start = format(BiocGenerics::start(ps) - 1L, scientific = FALSE, trim = TRUE),
    end   = format(BiocGenerics::end(ps),        scientific = FALSE, trim = TRUE),
    name  = sprintf("%s:%d-%d", sn, BiocGenerics::start(ps) - 1L, BiocGenerics::end(ps)),
    stringsAsFactors = FALSE
  )
  bed_path <- paste0(out_stem, "_consensus_peaks.bed")
  utils::write.table(bed, bed_path, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = FALSE)
  rep$peakset$bed_written <- bed_path
  add_md(sprintf("- consensus peak BED written to `%s` (0-based start, for cisTarget DB build)",
                 bed_path), "")
}

# ---------------------------------------------------------------- matrices
add_md("## Available matrices", "")
avail <- tryCatch(getAvailableMatrices(proj), error = function(e) character(0))
rep$available_matrices <- as.character(avail)
add_md(sprintf("- %s", if (length(avail)) paste(avail, collapse = ", ") else "none"), "")
rep$has_peak_matrix <- "PeakMatrix" %in% avail
rep$has_gene_score_matrix <- "GeneScoreMatrix" %in% avail
rep$has_gene_integration_matrix <- any(grepl("GeneIntegrationMatrix", avail))
if (!rep$has_peak_matrix) {
  add_md("> WARNING: no `PeakMatrix`. addPeakMatrix() must run before a",
         "> cell x peak count matrix can be exported for pycisTopic.", "")
}

# ---------------------------------------------------------------- fragments
add_md("## Fragment files", "")
af <- tryCatch(getArrowFiles(proj), error = function(e) character(0))
rep$arrow_files <- as.character(af)
rep$n_arrow_files <- length(af)
rep$arrow_files_exist <- as.logical(file.exists(af))
add_md(sprintf("- %d Arrow files, %d present on disk",
               length(af), sum(rep$arrow_files_exist)))
if (length(af)) {
  sz <- file.size(af[rep$arrow_files_exist])
  rep$arrow_total_gb <- round(sum(sz, na.rm = TRUE) / 1024^3, 2)
  add_md(sprintf("- total Arrow size: **%.2f GB**", rep$arrow_total_gb))
}
add_md("",
       "pycisTopic's canonical path wants per-pseudobulk fragments. Arrow files hold",
       "the fragments, so `export_fragments.R` in this repo can emit per-group",
       "fragment TSVs without needing the original CellRanger output.", "")

# ---------------------------------------------------------------- embeddings
add_md("## Embeddings / reduced dims", "")
emb <- tryCatch(names(proj@embeddings), error = function(e) character(0))
rd  <- tryCatch(names(proj@reducedDims), error = function(e) character(0))
rep$embeddings <- as.character(emb)
rep$reducedDims <- as.character(rd)
add_md(sprintf("- embeddings: %s", if (length(emb)) paste(emb, collapse = ", ") else "none"),
       sprintf("- reducedDims: %s", if (length(rd)) paste(rd, collapse = ", ") else "none"),
       "")
glue_like <- grep("glue|GLUE", c(emb, rd), value = TRUE)
rep$glue_like_embeddings <- as.character(glue_like)
if (length(glue_like)) {
  add_md(sprintf("- GLUE-like entries: %s", paste(glue_like, collapse = ", ")), "")
} else {
  add_md("> The GLUE latent space is not stored in this ArchRProject.",
         "> Confirm where the integrated embedding lives (an .h5ad, or a separate file);",
         "> metacell pairing reads it from there.", "")
}

# ---------------------------------------------------------------- write
json_path <- paste0(out_stem, ".json")
if (have_jsonlite) {
  writeLines(jsonlite::toJSON(rep, auto_unbox = TRUE, null = "null",
                              na = "null", pretty = TRUE), json_path)
} else {
  # Minimal fallback so the run still produces a machine-readable file.
  esc <- function(s) gsub('"', '\\\\"', as.character(s))
  ser <- function(x) {
    if (is.null(x) || (length(x) == 1 && is.na(x))) return("null")
    if (is.list(x)) {
      if (!is.null(names(x)))
        return(paste0("{", paste(sprintf('"%s":%s', esc(names(x)),
                                         vapply(x, ser, "")), collapse = ","), "}"))
      return(paste0("[", paste(vapply(x, ser, ""), collapse = ","), "]"))
    }
    if (is.logical(x) && length(x) == 1) return(if (x) "true" else "false")
    if (is.numeric(x) && length(x) == 1) return(as.character(x))
    if (length(x) == 1) return(sprintf('"%s"', esc(x)))
    paste0("[", paste(vapply(x, ser, ""), collapse = ","), "]")
  }
  writeLines(ser(rep), json_path)
}
writeLines(md, paste0(out_stem, ".md"))
message(sprintf("wrote %s and %s.md", json_path, out_stem))
