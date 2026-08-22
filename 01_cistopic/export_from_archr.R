#!/usr/bin/env Rscript
# =============================================================================
# 01_cistopic/export_from_archr.R
#
# Export the inputs pycisTopic needs out of an ArchRProject, without touching
# the project on disk.
#
# WHAT THIS WRITES (into --out-dir)
#   peak_matrix.mtx        Matrix Market, integer, REGIONS x CELLS.
#                          Orientation is deliberate: pycisTopic's
#                          create_cistopic_object() documents
#                          "sparse.csr_matrix containing cells as columns and
#                          regions as rows". Do not transpose it.
#   barcodes.tsv           One cell name per line, in .mtx COLUMN order.
#                          Names are ArchR cell names, i.e. "Sample#Barcode".
#   peaks.bed              4-column BED (chrom, start, end, name) in .mtx ROW
#                          order. `name` is the exact region string
#                          "chr:start-end" that must be handed to
#                          create_cistopic_object(region_names=...).
#   regions.txt            Just the region-name column of peaks.bed, one per
#                          line, in .mtx row order. Convenience for Python.
#   cell_metadata.tsv      cellColData as TSV, index column "cell_id",
#                          guaranteed to contain --group-col and "Sample".
#   export_manifest.tsv    Dimensions, nnz, checksums of intent, coordinate
#                          convention, ArchR/genome provenance. Read this
#                          before trusting anything else here.
#   fragments/             Only with --export-fragments. Per-group, per-sample
#                          unsorted BED-like fragment parts, plus
#                          finalize_fragments.sh which sorts + bgzips + tabixes
#                          them. See TRAP 4.
#
# TRAPS THIS SCRIPT EXISTS TO NAVIGATE
#
# TRAP 1 - COORDINATE BASE. ArchR keeps its peak set as a GRanges: 1-based,
#   both ends inclusive. A 501 bp fixed-width ArchR peak therefore has
#   end - start == 500. BED, and every tool downstream of pycisTopic
#   (pyranges blacklist overlap inside create_cistopic_object, the region-set
#   BEDs that SCENIC+ / pycistarget consume, bedtools, the cisTarget database
#   builder), is 0-based half-open, where the same peak has
#   end - start == 501. If you export the ArchR numbers verbatim you get
#   region names that are silently shifted 1 bp relative to every BED file you
#   later intersect them with, and 501 bp peaks that report as 500 bp.
#   Default here is --coord-base bed0, which emits start = start(gr) - 1 and
#   end = end(gr). The manifest records which base was used. --coord-base
#   archr1 is available only for reproducing an older export; it is not what
#   you want for a fresh SCENIC+ run.
#
# TRAP 2 - REGION ROW ORDER MUST BE FIXED ACROSS CHUNKS. We read the
#   PeakMatrix one Arrow file at a time (see TRAP 3), and ArchR's
#   .getFeatureDF() re-orders features by a split on seqnames. getMatrixFromArrow
#   returns rows in that featureDF order, NOT sorted, and NOT necessarily the
#   order of getPeakSet(). getMatrixFromProject() re-sorts afterwards, but we
#   are not calling it. So this script pins the region order from the FIRST
#   chunk it reads and then asserts, chunk by chunk, that the region names are
#   identical and in the same order. A mismatch is a hard error, never a
#   reorder-and-hope: a silent row permutation between samples would scramble
#   the count matrix in a way no downstream QC would catch.
#
# TRAP 3 - MEMORY. A hematopoietic reference is routinely 100k+ cells x 300k
#   peaks. getMatrixFromProject() builds every per-sample SummarizedExperiment
#   and then cbind()s the assays, so peak RSS is roughly the whole matrix twice
#   over; on ARC that is an OOM kill with no useful message. This script never
#   calls getMatrixFromProject(). It walks Arrow files, and within each Arrow
#   file walks cells in blocks of --chunk-size, converts each block's dgCMatrix
#   to (i, j, x) triplets, appends them to a plain text body file, and frees the
#   block. Only one chunk is resident at a time.
#
# TRAP 4 - MATRIX MARKET HEADER COMES LAST. The .mtx format wants
#   "nrow ncol nnz" on line 2, but nnz is unknown until the final chunk is
#   written. We stream triplets to peak_matrix.mtx.body, count nnz as we go,
#   then write the 2-line header to peak_matrix.mtx and cat the body onto it.
#   If the script dies mid-export you are left with a .body and no .mtx, which
#   is the correct failure mode - a truncated .mtx with a plausible header
#   would load fine and be wrong.
#
# TRAP 5 - FRAGMENTS ARE NOT SORTED AND MUST NOT BE bgzipped HERE.
#   getFragmentsFromArrow() returns a GRanges per Arrow file whose cell name
#   lives in mcols(gr)$RG as "Sample#Barcode". Concatenating those across
#   samples gives you coordinate-unsorted output. tabix requires
#   sort -k1,1 -k2,2n before bgzip. Sorting 10^9 fragment rows inside R would
#   blow memory, so this script writes unsorted parts and emits
#   finalize_fragments.sh to do the sort/bgzip/tabix with coreutils, which
#   spills to disk. Run that script (set TMPDIR to node-local scratch) before
#   pointing anything at the fragment files. Fragment export is also the
#   slowest thing here by an order of magnitude, which is why it is opt-in.
#
# TRAP 6 - THIS SCRIPT NEVER WRITES TO THE ArchRProject. saveArchRProject() is
#   never called, no addCellColData(), no addPeakMatrix(). The project is
#   opened read-only in spirit; if a collaborator's job is reading the same
#   project directory concurrently, nothing here will corrupt it. The only
#   thing ArchR itself may write is its own log file, which we redirect into
#   --out-dir.
#
# TRAP 7 - --group-col IS VALIDATED UP FRONT, BEFORE THE EXPENSIVE PART.
#   Getting 40 minutes into a matrix export and then discovering the cell-type
#   column is called "CellType" and not "Cell_Type" is the single most common
#   way this step is wasted. The column is checked, and its NA count and level
#   table are logged, in the first few seconds.
#
# TYPICAL INVOCATION (inside a SLURM job on ARC)
#   Rscript 01_cistopic/export_from_archr.R \
#     --archr-project /path/to/ArchR_hema \
#     --out-dir       /path/to/scenicplus/01_cistopic/archr_export \
#     --group-col     CellType \
#     --chunk-size    5000 \
#     --threads       8
#
#   add --export-fragments (and expect it to dominate the runtime) if you also
#   want per-cell-type fragment files for pseudobulk / bigwig work.
#
# EXIT CODES: 0 ok. 1 any validation or IO failure, always with a message on
# stderr prefixed "FATAL:". There are no silent fallbacks and no partial
# successes that return 0.
# =============================================================================

suppressWarnings(suppressMessages({
  library(ArchR)
  library(Matrix)
  library(SummarizedExperiment)
  library(GenomicRanges)
  library(S4Vectors)
}))

# ---------------------------------------------------------------------------
# argument parsing: tiny getarg helper over commandArgs(), same style as the
# other scripts in this repo
# ---------------------------------------------------------------------------
.ARGV <- commandArgs(trailingOnly = TRUE)

die <- function(...) {
  cat("FATAL: ", paste0(..., collapse = ""), "\n", sep = "", file = stderr())
  quit(save = "no", status = 1L)
}

logmsg <- function(...) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), paste0(..., collapse = ""),
      "\n", sep = "", file = stderr())
  flush(stderr())
}

#' Fetch "--flag value" from the command line.
#' required=TRUE and no value -> hard error (never a silent default).
getarg <- function(flag, default = NULL, required = FALSE) {
  hit <- which(.ARGV == flag)
  if (length(hit) > 1L) die(flag, " given ", length(hit), " times; give it once")
  if (length(hit) == 0L) {
    if (required) die("missing required argument ", flag)
    return(default)
  }
  if (hit + 1L > length(.ARGV)) die(flag, " given with no value")
  val <- .ARGV[hit + 1L]
  if (startsWith(val, "--")) die(flag, " given with no value (next token is ", val, ")")
  val
}

#' Presence-only flag.
hasflag <- function(flag) any(.ARGV == flag)

getarg_int <- function(flag, default) {
  v <- getarg(flag, NULL)
  if (is.null(v)) return(as.integer(default))
  n <- suppressWarnings(as.integer(v))
  if (is.na(n)) die(flag, " must be an integer, got '", v, "'")
  n
}

if (hasflag("--help") || hasflag("-h")) {
  cat(
    "Usage: Rscript export_from_archr.R --archr-project DIR --out-dir DIR --group-col COL [options]\n\n",
    "Required:\n",
    "  --archr-project DIR   ArchRProject directory (containing Save-ArchR-Project.rds)\n",
    "                        or a direct path to a .rds ArchRProject.\n",
    "  --out-dir DIR         output directory (created if absent)\n",
    "  --group-col COL       cellColData column holding the cell-type label\n\n",
    "Options:\n",
    "  --use-matrix NAME     matrix to export (default: PeakMatrix)\n",
    "  --chunk-size N        cells per export chunk (default: 5000)\n",
    "  --threads N           ArchR threads (default: 1; only used for fragment export)\n",
    "  --coord-base MODE     bed0 (default, 0-based half-open) | archr1 (1-based inclusive)\n",
    "  --extra-cols A,B,C    additional cellColData columns to hard-require\n",
    "  --export-fragments    also export per-group fragment parts (slow)\n",
    "  --fragment-groups A,B only export fragments for these groups (default: all)\n",
    "  --keep-body           keep peak_matrix.mtx.body after assembling the .mtx\n",
    "  --help                this message\n",
    sep = ""
  )
  quit(save = "no", status = 0L)
}

archr_path  <- getarg("--archr-project", required = TRUE)
out_dir     <- getarg("--out-dir",       required = TRUE)
group_col   <- getarg("--group-col",     required = TRUE)
use_matrix  <- getarg("--use-matrix", "PeakMatrix")
chunk_size  <- getarg_int("--chunk-size", 5000L)
threads     <- getarg_int("--threads", 1L)
coord_base  <- getarg("--coord-base", "bed0")
extra_cols  <- getarg("--extra-cols", "")
do_frags    <- hasflag("--export-fragments")
frag_groups <- getarg("--fragment-groups", "")
keep_body   <- hasflag("--keep-body")

if (!coord_base %in% c("bed0", "archr1")) {
  die("--coord-base must be 'bed0' or 'archr1', got '", coord_base, "'")
}
if (chunk_size < 1L) die("--chunk-size must be >= 1")
if (threads < 1L)    die("--threads must be >= 1")

extra_cols <- if (nzchar(extra_cols)) trimws(strsplit(extra_cols, ",")[[1]]) else character(0)
frag_groups <- if (nzchar(frag_groups)) trimws(strsplit(frag_groups, ",")[[1]]) else character(0)

# ---------------------------------------------------------------------------
# 0. environment / output dir
# ---------------------------------------------------------------------------
logmsg("export_from_archr.R starting")
logmsg("  R           : ", R.version.string)
logmsg("  ArchR       : ", as.character(utils::packageVersion("ArchR")))
logmsg("  Matrix      : ", as.character(utils::packageVersion("Matrix")))
logmsg("  coord base  : ", coord_base,
       if (coord_base == "bed0") "  (start-1, BED half-open -- recommended)"
       else "  (verbatim ArchR 1-based inclusive -- NOT BED compatible)")

if (!dir.exists(out_dir)) {
  if (!dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)) {
    die("could not create --out-dir ", out_dir)
  }
}
if (file.access(out_dir, mode = 2L) != 0L) die("--out-dir not writable: ", out_dir)

# Keep ArchR's own logs next to the export instead of the submit directory.
addArchRLogging(useLogs = TRUE)
addArchRThreads(threads = threads)
archr_log <- file.path(out_dir, "ArchR_export.log")

f_mtx      <- file.path(out_dir, "peak_matrix.mtx")
f_body     <- file.path(out_dir, "peak_matrix.mtx.body")
f_barcodes <- file.path(out_dir, "barcodes.tsv")
f_peaks    <- file.path(out_dir, "peaks.bed")
f_regions  <- file.path(out_dir, "regions.txt")
f_meta     <- file.path(out_dir, "cell_metadata.tsv")
f_manifest <- file.path(out_dir, "export_manifest.tsv")

# Refuse to half-overwrite a previous export: a stale .mtx next to a fresh
# barcodes.tsv is undetectable downstream.
for (f in c(f_mtx, f_body)) if (file.exists(f)) {
  logmsg("removing stale ", basename(f))
  if (!file.remove(f)) die("could not remove stale ", f)
}

# ---------------------------------------------------------------------------
# 1. load the project (read-only; saveArchRProject is never called)
# ---------------------------------------------------------------------------
logmsg("loading ArchRProject from ", archr_path)
if (!file.exists(archr_path)) die("--archr-project does not exist: ", archr_path)

proj <- tryCatch({
  if (dir.exists(archr_path)) {
    loadArchRProject(path = archr_path, showLogo = FALSE)
  } else if (grepl("\\.rds$", archr_path, ignore.case = TRUE)) {
    readRDS(archr_path)
  } else {
    die("--archr-project must be a project directory or a .rds file: ", archr_path)
  }
}, error = function(e) die("failed to load ArchRProject: ", conditionMessage(e)))

if (!inherits(proj, "ArchRProject")) {
  die("loaded object is class '", paste(class(proj), collapse = "/"),
      "', not ArchRProject")
}

cell_names_all <- getCellNames(proj)
n_cells_total  <- length(cell_names_all)
arrow_files    <- getArrowFiles(proj)
logmsg("  cells       : ", n_cells_total)
logmsg("  samples     : ", length(arrow_files), " (", paste(names(arrow_files), collapse = ", "), ")")
logmsg("  genome      : ", tryCatch(getGenome(proj), error = function(e) "<unknown>"))

missing_arrows <- arrow_files[!file.exists(arrow_files)]
if (length(missing_arrows)) {
  die("Arrow files referenced by the project do not exist on this filesystem:\n  ",
      paste(missing_arrows, collapse = "\n  "),
      "\n  (an ArchRProject moved between machines keeps absolute Arrow paths; ",
      "re-point them before exporting)")
}

# ---------------------------------------------------------------------------
# 2. validate --group-col and metadata BEFORE the expensive export (TRAP 7)
# ---------------------------------------------------------------------------
ccd <- as.data.frame(getCellColData(proj), stringsAsFactors = FALSE)
if (!identical(rownames(ccd), cell_names_all)) {
  die("cellColData rownames do not match getCellNames(); refusing to guess the ",
      "cell order")
}

if (!group_col %in% colnames(ccd)) {
  cat("FATAL: --group-col '", group_col, "' is not a column of cellColData.\n",
      sep = "", file = stderr())
  cat("Available columns (", length(colnames(ccd)), "):\n", sep = "", file = stderr())
  cat(paste0("  - ", colnames(ccd), collapse = "\n"), "\n", sep = "", file = stderr())
  near <- colnames(ccd)[agrepl(group_col, colnames(ccd), ignore.case = TRUE, max.distance = 0.35)]
  if (length(near)) {
    cat("Did you mean: ", paste(near, collapse = ", "), " ?\n", sep = "", file = stderr())
  }
  quit(save = "no", status = 1L)
}

if (!"Sample" %in% colnames(ccd)) {
  die("cellColData has no 'Sample' column. Every ArchRProject built by ",
      "createArrowFiles() has one; this project looks hand-modified. ",
      "pycisTopic needs a per-sample key for split_pattern handling.")
}

missing_extra <- setdiff(extra_cols, colnames(ccd))
if (length(missing_extra)) {
  die("--extra-cols not present in cellColData: ", paste(missing_extra, collapse = ", "))
}

grp <- as.character(ccd[[group_col]])
n_na_grp <- sum(is.na(grp) | !nzchar(grp))
logmsg("group column '", group_col, "': ", length(unique(grp[!is.na(grp)])),
       " levels, ", n_na_grp, " cells with NA/empty label")
grp_tab <- sort(table(grp, useNA = "ifany"), decreasing = TRUE)
for (i in seq_along(grp_tab)) {
  logmsg("    ", format(names(grp_tab)[i], width = 30), " ", grp_tab[[i]])
}
if (n_na_grp > 0L) {
  logmsg("NOTE: ", n_na_grp, " cells carry no ", group_col, " label. They are ",
         "exported anyway (the matrix must stay complete); pycisTopic's ",
         "find_diff_features() drops NA groups itself via dropna().")
}
if (length(unique(grp[!is.na(grp)])) < 2L) {
  die("--group-col '", group_col, "' has fewer than 2 non-NA levels; ",
      "per-cell-type DARs downstream would be undefined")
}

available <- getAvailableMatrices(proj)
if (!use_matrix %in% available) {
  die("--use-matrix '", use_matrix, "' not in the Arrow files. Available: ",
      paste(available, collapse = ", "),
      if (use_matrix == "PeakMatrix")
        "  (run addPeakMatrix(proj) in your ArchR session first -- this script will not modify the project)"
      else "")
}

peak_set <- tryCatch(getPeakSet(proj), error = function(e) NULL)
if (is.null(peak_set) || length(peak_set) == 0L) {
  die("getPeakSet() is empty. A reproducible peak set is required before ",
      "PeakMatrix export.")
}
pk_widths <- unique(width(peak_set))
logmsg("peak set     : ", length(peak_set), " peaks, width(s): ",
       paste(head(pk_widths, 5), collapse = ","),
       if (length(pk_widths) > 5) ",..." else "")
if (length(pk_widths) == 1L) {
  logmsg("  fixed-width peak set confirmed (", pk_widths, " bp as ArchR GRanges width)")
} else {
  logmsg("WARNING: peak set is not fixed-width. SCENIC+ tolerates this, but an ",
         "ArchR reproducible peak set normally is fixed-width; check that this ",
         "is the peak set you meant.")
}

# ---------------------------------------------------------------------------
# 3. helpers
# ---------------------------------------------------------------------------

#' GRanges/featureDF rows -> pycisTopic region names + BED coordinates.
#' Applies the --coord-base decision in exactly one place (TRAP 1).
.make_region_frame <- function(chrom, start1, end1) {
  if (coord_base == "bed0") {
    bed_start <- as.integer(start1) - 1L
    bed_end   <- as.integer(end1)
  } else {
    bed_start <- as.integer(start1)
    bed_end   <- as.integer(end1)
  }
  if (any(bed_start < 0L)) die("negative start after coordinate conversion; ",
                               "a peak starts at position 0/1 -- inspect the peak set")
  data.frame(
    chrom = as.character(chrom),
    start = bed_start,
    end   = bed_end,
    name  = paste0(as.character(chrom), ":", bed_start, "-", bed_end),
    stringsAsFactors = FALSE
  )
}

#' Region names from the rowRanges/rowData of a chunk SummarizedExperiment.
.se_region_frame <- function(se) {
  rr <- SummarizedExperiment::rowRanges(se)
  if (inherits(rr, "GRanges") && length(rr) > 0L) {
    return(.make_region_frame(seqnames(rr), start(rr), end(rr)))
  }
  # Fallback: ArchR returns a plain featureDF when makeGRangesFromDataFrame
  # fails. Columns are seqnames/idx/start/end (see ArchR MatrixFeatures.R).
  rd <- as.data.frame(SummarizedExperiment::rowData(se), stringsAsFactors = FALSE)
  if (!all(c("seqnames", "start", "end") %in% colnames(rd))) {
    die("cannot derive region coordinates: chunk rowData has columns ",
        paste(colnames(rd), collapse = ","),
        " and rowRanges is not a GRanges")
  }
  .make_region_frame(rd$seqnames, rd$start, rd$end)
}

# ---------------------------------------------------------------------------
# 4. stream the matrix out, chunked over cells (TRAP 2, TRAP 3, TRAP 4)
# ---------------------------------------------------------------------------
logmsg("exporting ", use_matrix, " in chunks of ", chunk_size,
       " cells; streaming triplets to ", basename(f_body))

region_frame   <- NULL      # pinned on first chunk, then asserted
n_regions      <- NA_integer_
nnz_total      <- 0
col_offset     <- 0L
barcode_order  <- character(0)
body_con <- file(f_body, open = "wt")
on.exit({ if (isOpen(body_con)) close(body_con) }, add = TRUE)

t_start <- Sys.time()
chunk_i <- 0L

for (s in seq_along(arrow_files)) {
  sample_id <- names(arrow_files)[s]
  arrow     <- arrow_files[[s]]

  # Cells of this sample that are actually in the project. ArchR cell names are
  # "Sample#Barcode", so the sample prefix is the authority here, not
  # ccd$Sample, in case someone relabelled Sample post hoc.
  cells_s <- cell_names_all[startsWith(cell_names_all, paste0(sample_id, "#"))]
  if (length(cells_s) == 0L) {
    logmsg("  [", s, "/", length(arrow_files), "] ", sample_id,
           ": 0 cells retained in the project, skipping")
    next
  }
  logmsg("  [", s, "/", length(arrow_files), "] ", sample_id, ": ",
         length(cells_s), " cells")

  starts <- seq(1L, length(cells_s), by = chunk_size)
  for (st in starts) {
    en        <- min(st + chunk_size - 1L, length(cells_s))
    cells_blk <- cells_s[st:en]
    chunk_i   <- chunk_i + 1L

    se <- tryCatch(
      getMatrixFromArrow(
        ArrowFile = arrow,
        useMatrix = use_matrix,
        cellNames = cells_blk,
        ArchRProj = proj,
        verbose   = FALSE,
        binarize  = FALSE,
        logFile   = archr_log
      ),
      error = function(e)
        die("getMatrixFromArrow failed for ", sample_id, " cells ", st, "-", en,
            ": ", conditionMessage(e))
    )
    if (is.null(se)) {
      die("getMatrixFromArrow returned NULL for ", sample_id, " cells ", st, "-", en,
          " (none of the requested cells are in ", use_matrix,
          " for this Arrow file)")
    }

    assay_names <- names(SummarizedExperiment::assays(se))
    if (!use_matrix %in% assay_names) {
      die("chunk SummarizedExperiment has assays [",
          paste(assay_names, collapse = ","), "] but not '", use_matrix,
          "'. A Sparse.Assays.Matrix (e.g. MotifMatrix) cannot be exported by ",
          "this script; export PeakMatrix or TileMatrix instead.")
    }
    m <- SummarizedExperiment::assays(se)[[use_matrix]]
    if (!inherits(m, "dgCMatrix")) m <- methods::as(m, "dgCMatrix")

    # --- region order: pin once, assert forever after (TRAP 2) -------------
    rf <- .se_region_frame(se)
    if (is.null(region_frame)) {
      region_frame <- rf
      n_regions    <- nrow(rf)
      if (anyDuplicated(region_frame$name)) {
        dups <- region_frame$name[duplicated(region_frame$name)]
        die("duplicated region names in ", use_matrix, " (e.g. ",
            paste(head(dups, 3), collapse = ", "),
            "). pycisTopic indexes regions by name and will mis-subset.")
      }
      if (n_regions != nrow(m)) {
        die("region frame has ", n_regions, " rows but the matrix has ", nrow(m))
      }
      logmsg("    region axis pinned from first chunk: ", n_regions, " regions, ",
             "first = ", region_frame$name[1], ", last = ",
             region_frame$name[n_regions])
      if (length(peak_set) != n_regions) {
        logmsg("    NOTE: getPeakSet() has ", length(peak_set), " peaks but ",
               use_matrix, " has ", n_regions, " rows. This is expected when ",
               "the peak set contains chromosomes excluded at matrix build ",
               "time; the matrix is the authority for row order.")
      }
    } else {
      if (nrow(rf) != n_regions || !identical(rf$name, region_frame$name)) {
        die("region order/content differs between chunks. Chunk ", chunk_i,
            " (", sample_id, " cells ", st, "-", en, ") has ", nrow(rf),
            " regions vs ", n_regions, " pinned",
            if (nrow(rf) == n_regions)
              paste0("; first mismatch at row ",
                     which(rf$name != region_frame$name)[1], ": '",
                     rf$name[which(rf$name != region_frame$name)[1]], "' vs '",
                     region_frame$name[which(rf$name != region_frame$name)[1]], "'")
            else "",
            ". Refusing to reorder rows silently -- rebuild ", use_matrix,
            " so all Arrow files share one featureDF.")
      }
    }

    # --- column identity: trust the matrix, not the request order ----------
    blk_cells <- colnames(m)
    if (is.null(blk_cells)) die("chunk matrix has no column names")
    if (!setequal(blk_cells, cells_blk)) {
      die("chunk ", chunk_i, " returned ", length(blk_cells), " cells for ",
          length(cells_blk), " requested; cannot align barcodes")
    }
    barcode_order <- c(barcode_order, blk_cells)

    # --- triplets ---------------------------------------------------------
    trip <- Matrix::summary(m)          # data.frame(i, j, x) on nonzeros only
    if (nrow(trip) > 0L) {
      if (any(trip$x <= 0)) {
        die("chunk ", chunk_i, " contains non-positive stored values; ",
            use_matrix, " should be non-negative counts")
      }
      utils::write.table(
        data.frame(i = trip$i,
                   j = trip$j + col_offset,
                   x = as.integer(trip$x)),
        file = body_con, sep = " ", row.names = FALSE, col.names = FALSE,
        quote = FALSE
      )
      nnz_total <- nnz_total + nrow(trip)
    }
    col_offset <- col_offset + ncol(m)

    if (chunk_i %% 5L == 0L || en == length(cells_s)) {
      el <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
      logmsg("    chunk ", chunk_i, ": ", col_offset, "/", n_cells_total,
             " cells, ", format(nnz_total, big.mark = ","), " nonzeros, ",
             sprintf("%.1f", el), " min elapsed")
    }
    rm(se, m, trip); invisible(gc(verbose = FALSE, full = FALSE))
  }
}

close(body_con)

if (is.null(region_frame)) die("no chunk produced any data; nothing exported")
if (col_offset == 0L)      die("0 cells exported; nothing to write")
if (nnz_total == 0)        die("matrix is entirely zero; refusing to write it")
if (length(barcode_order) != col_offset) {
  die("barcode bookkeeping mismatch: ", length(barcode_order), " names vs ",
      col_offset, " columns")
}
if (anyDuplicated(barcode_order)) {
  die("duplicated cell names in the export (e.g. ",
      barcode_order[duplicated(barcode_order)][1], ")")
}
if (col_offset != n_cells_total) {
  logmsg("NOTE: exported ", col_offset, " of ", n_cells_total,
         " project cells. The difference is cells absent from ", use_matrix,
         " in their Arrow file. cell_metadata.tsv is subset to match.")
}

# ---------------------------------------------------------------------------
# 5. assemble the .mtx (header first, then body) - TRAP 4
# ---------------------------------------------------------------------------
logmsg("assembling ", basename(f_mtx), ": ", n_regions, " x ", col_offset,
       ", ", format(nnz_total, big.mark = ","), " nonzeros")
hdr <- file(f_mtx, open = "wt")
writeLines("%%MatrixMarket matrix coordinate integer general", hdr)
writeLines(paste("% regions(rows) x cells(cols); region names in",
                 "regions.txt/peaks.bed, cell names in barcodes.tsv;",
                 "coord base:", coord_base), hdr)
writeLines(paste(n_regions, col_offset, nnz_total), hdr)
close(hdr)

# Append the body with a shell append (>>) so the 3-line header survives.
# Note: system2(stdout=<path>) TRUNCATES, so it cannot be used here.
rc <- system2("/bin/sh", args = c("-c",
  shQuote(paste("cat", shQuote(f_body), ">>", shQuote(f_mtx)))))
if (rc != 0L) die("failed to append matrix body to ", f_mtx, " (cat exit ", rc, ")")

# Cheap structural check: line count must be 3 header lines + nnz.
n_lines <- suppressWarnings(as.numeric(
  strsplit(trimws(system2("wc", c("-l", shQuote(f_mtx)), stdout = TRUE)), "\\s+")[[1]][1]))
if (!is.na(n_lines) && n_lines != nnz_total + 3) {
  die("assembled ", basename(f_mtx), " has ", n_lines, " lines, expected ",
      nnz_total + 3, " (3 header + nnz). Do not use this file.")
}
if (!keep_body) invisible(file.remove(f_body))

# ---------------------------------------------------------------------------
# 6. barcodes, peaks.bed, regions.txt, metadata
# ---------------------------------------------------------------------------
logmsg("writing ", basename(f_barcodes))
writeLines(barcode_order, f_barcodes)

logmsg("writing ", basename(f_peaks), " and ", basename(f_regions))
utils::write.table(region_frame[, c("chrom", "start", "end", "name")],
                   file = f_peaks, sep = "\t", quote = FALSE,
                   row.names = FALSE, col.names = FALSE)
writeLines(region_frame$name, f_regions)

logmsg("writing ", basename(f_meta))
meta <- ccd[barcode_order, , drop = FALSE]
if (nrow(meta) != length(barcode_order)) {
  die("cellColData subset lost rows; barcode/metadata mismatch")
}
# Guarantee the columns pycisTopic will be asked for, first and by name.
front <- unique(c(group_col, "Sample", extra_cols))
meta  <- meta[, c(front, setdiff(colnames(meta), front)), drop = FALSE]
meta  <- cbind(cell_id = barcode_order, meta)
# List columns (ArchR sometimes stores them) break read_table on the Python side.
listy <- names(meta)[vapply(meta, function(x) is.list(x) || !is.null(dim(x)), logical(1))]
if (length(listy)) {
  logmsg("NOTE: dropping non-atomic cellColData columns (unrepresentable in TSV): ",
         paste(listy, collapse = ", "))
  meta <- meta[, setdiff(colnames(meta), listy), drop = FALSE]
}
utils::write.table(meta, file = f_meta, sep = "\t", quote = FALSE,
                   row.names = FALSE, col.names = TRUE, na = "NA")

# ---------------------------------------------------------------------------
# 7. optional per-group fragment export (TRAP 5)
# ---------------------------------------------------------------------------
frag_dir <- file.path(out_dir, "fragments")
if (do_frags) {
  logmsg("--export-fragments set: exporting per-group fragment parts (slow)")
  dir.create(frag_dir, recursive = TRUE, showWarnings = FALSE)

  grp_vec <- as.character(meta[[group_col]])
  names(grp_vec) <- meta$cell_id
  groups_present <- sort(unique(grp_vec[!is.na(grp_vec) & nzchar(grp_vec)]))
  groups_to_do <- if (length(frag_groups)) {
    bad <- setdiff(frag_groups, groups_present)
    if (length(bad)) die("--fragment-groups not found in ", group_col, ": ",
                         paste(bad, collapse = ", "),
                         ". Present: ", paste(groups_present, collapse = ", "))
    frag_groups
  } else groups_present
  logmsg("  groups: ", paste(groups_to_do, collapse = ", "))

  # Sanitize group names for filenames but keep a mapping file so the link back
  # to the label is never guesswork.
  safe <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
  utils::write.table(
    data.frame(group = groups_to_do, file_stem = safe(groups_to_do)),
    file = file.path(frag_dir, "group_file_map.tsv"), sep = "\t",
    quote = FALSE, row.names = FALSE
  )

  for (s in seq_along(arrow_files)) {
    sample_id <- names(arrow_files)[s]
    arrow     <- arrow_files[[s]]
    cells_s   <- names(grp_vec)[startsWith(names(grp_vec), paste0(sample_id, "#"))]
    if (!length(cells_s)) next
    logmsg("  fragments from ", sample_id, " (", length(cells_s), " cells)")

    gr <- tryCatch(
      getFragmentsFromArrow(ArrowFile = arrow, cellNames = cells_s,
                            verbose = FALSE, logFile = archr_log),
      error = function(e)
        die("getFragmentsFromArrow failed for ", sample_id, ": ", conditionMessage(e))
    )
    if (is.null(gr) || length(gr) == 0L) {
      logmsg("    no fragments returned, skipping")
      next
    }
    if (!"RG" %in% colnames(mcols(gr))) {
      die("fragment GRanges from ", sample_id, " has no RG metadata column; ",
          "cannot assign fragments to cells")
    }
    rg  <- as.character(mcols(gr)$RG)
    gsp <- grp_vec[rg]

    for (g in groups_to_do) {
      idx <- which(!is.na(gsp) & gsp == g)
      if (!length(idx)) next
      part <- file.path(frag_dir, paste0(safe(g), "__", safe(sample_id), ".part.tsv"))
      # BED-style: 0-based start regardless of --coord-base. Fragment files are
      # consumed by tabix/scatac_fragment_tools, which are unconditionally BED.
      utils::write.table(
        data.frame(chrom = as.character(seqnames(gr))[idx],
                   start = start(gr)[idx] - 1L,
                   end   = end(gr)[idx],
                   cell  = rg[idx],
                   count = 1L),
        file = part, sep = "\t", quote = FALSE, row.names = FALSE,
        col.names = FALSE
      )
      logmsg("    ", g, ": ", length(idx), " fragments -> ", basename(part))
    }
    rm(gr, rg, gsp); invisible(gc(verbose = FALSE, full = FALSE))
  }

  sh <- file.path(frag_dir, "finalize_fragments.sh")
  writeLines(c(
    "#!/usr/bin/env bash",
    "# Generated by export_from_archr.R -- sort, bgzip and tabix the fragment parts.",
    "#",
    "# Fragments were written UNSORTED, one part per (group, sample). tabix needs",
    "# coordinate-sorted input, and sorting inside R would not fit in memory for a",
    "# 100k-cell reference. GNU sort spills to disk, so it does. Point TMPDIR at",
    "# node-local scratch (e.g. $SLURM_TMPDIR) or sort will fill /tmp.",
    "#",
    "# Requires: sort (coreutils), bgzip + tabix (htslib/samtools module).",
    "set -euo pipefail",
    "cd \"$(dirname \"$0\")\"",
    ": \"${TMPDIR:=/tmp}\"",
    "echo \"using TMPDIR=$TMPDIR\" >&2",
    "command -v bgzip >/dev/null || { echo 'FATAL: bgzip not on PATH (module load htslib)' >&2; exit 1; }",
    "command -v tabix >/dev/null || { echo 'FATAL: tabix not on PATH (module load htslib)' >&2; exit 1; }",
    "shopt -s nullglob",
    "stems=$(for f in *.part.tsv; do echo \"${f%%__*}\"; done | sort -u)",
    "if [[ -z \"$stems\" ]]; then echo 'FATAL: no *.part.tsv files here' >&2; exit 1; fi",
    "for g in $stems; do",
    "  echo \"[fragments] $g\" >&2",
    "  parts=( \"$g\"__*.part.tsv )",
    "  LC_ALL=C sort -T \"$TMPDIR\" -S 4G -k1,1 -k2,2n \"${parts[@]}\" \\",
    "    | bgzip -c > \"${g}.fragments.tsv.gz\"",
    "  tabix -p bed -f \"${g}.fragments.tsv.gz\"",
    "  echo \"[fragments] $g done: $(zcat \"${g}.fragments.tsv.gz\" | wc -l) fragments\" >&2",
    "done",
    "echo 'all groups finalized; remove *.part.tsv when satisfied' >&2"
  ), sh)
  Sys.chmod(sh, mode = "0755")
  logmsg("wrote ", sh, " -- RUN IT before using the fragment files")
} else {
  logmsg("--export-fragments not set: skipping fragment export")
}

# ---------------------------------------------------------------------------
# 8. manifest
# ---------------------------------------------------------------------------
manifest <- data.frame(
  key = c("script", "generated", "archr_project", "archr_version", "genome",
          "use_matrix", "coord_base", "region_name_format",
          "n_regions", "n_cells_exported", "n_cells_in_project", "nnz",
          "matrix_orientation", "group_col", "n_groups", "n_cells_na_group",
          "peak_widths_archr", "chunk_size", "fragments_exported",
          "mtx_file", "barcodes_file", "peaks_file", "regions_file",
          "metadata_file"),
  value = c("01_cistopic/export_from_archr.R",
            format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
            normalizePath(archr_path, mustWork = FALSE),
            as.character(utils::packageVersion("ArchR")),
            tryCatch(getGenome(proj), error = function(e) "unknown"),
            use_matrix, coord_base,
            if (coord_base == "bed0") "chr:start0-end (BED half-open)"
            else "chr:start1-end (ArchR 1-based inclusive)",
            n_regions, col_offset, n_cells_total, nnz_total,
            "regions_rows_x_cells_cols", group_col,
            length(unique(grp[!is.na(grp)])), n_na_grp,
            paste(head(pk_widths, 5), collapse = ","), chunk_size,
            as.character(do_frags),
            basename(f_mtx), basename(f_barcodes), basename(f_peaks),
            basename(f_regions), basename(f_meta)),
  stringsAsFactors = FALSE
)
utils::write.table(manifest, file = f_manifest, sep = "\t", quote = FALSE,
                   row.names = FALSE, col.names = TRUE)

logmsg("wrote ", basename(f_manifest))
logmsg("DONE in ", sprintf("%.1f", as.numeric(difftime(Sys.time(), t_start, units = "mins"))),
       " min. ", n_regions, " regions x ", col_offset, " cells, ",
       format(nnz_total, big.mark = ","), " nonzeros.")
logmsg("Next: run_cistopic.py --matrix ", f_mtx, " --barcodes ", f_barcodes,
       " --regions ", f_regions, " --cell-metadata ", f_meta,
       " --group-col ", group_col)
logmsg("saveArchRProject() was NOT called; the project on disk is unchanged.")
quit(save = "no", status = 0L)
