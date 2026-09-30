#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(edgeR))

option_list <- list(
  make_option(c("--counts"), type = "character", help = "nf-core mature_counts.csv"),
  make_option(c("--logtpm"), type = "character", help = "nf-core mature_logtpm.csv"),
  make_option(c("--metadata"), type = "character", help = "Metadata CSV with sample_id and condition"),
  make_option(c("--controls"), type = "character", default = "", help = "Optional comma-separated control sample IDs"),
  make_option(c("--outdir"), type = "character", default = ".")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(opt$counts)) stop("Counts file not found: ", opt$counts)
if (!file.exists(opt$logtpm)) stop("logTPM file not found: ", opt$logtpm)
if (!file.exists(opt$metadata)) stop("Metadata file not found: ", opt$metadata)

# ---------------------------------------------------------
# Helpers
# ---------------------------------------------------------
normalize_id <- function(x) chartr(".", "-", trimws(as.character(x)))

read_metadata <- function(file, sample_ids, controls = "") {
  md <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)

  required <- c("sample_id", "condition")
  missing <- setdiff(required, colnames(md))
  if (length(missing) > 0) {
    stop("Metadata must contain columns: sample_id, condition. Missing: ", paste(missing, collapse = ", "))
  }

  md$sample_id <- trimws(as.character(md$sample_id))
  md$condition <- trimws(as.character(md$condition))

  if (anyDuplicated(md$sample_id)) {
    stop("Metadata contains duplicated sample_id values.")
  }
  if (any(!nzchar(md$condition))) {
    stop("Metadata contains empty condition values.")
  }

  idx <- match(sample_ids, md$sample_id)

  # Fallback for R/CSV name conversion such as sample.01 -> sample-01.
  if (anyNA(idx)) {
    md_norm <- normalize_id(md$sample_id)
    ids_norm <- normalize_id(sample_ids)
    if (anyDuplicated(md_norm)) {
      stop("Metadata sample IDs become duplicated after '.' -> '-' normalization.")
    }
    idx_norm <- match(ids_norm, md_norm)
    if (anyNA(idx_norm)) {
      missing_ids <- sample_ids[is.na(idx_norm)]
      stop("These count-matrix sample IDs are missing from metadata: ", paste(missing_ids, collapse = ", "))
    }
    idx <- idx_norm
    cat("Used fallback sample-ID matching with '.' converted to '-'.\n")
  }

  md <- md[idx, , drop = FALSE]
  rownames(md) <- sample_ids

  conditions <- unique(md$condition)
  if (length(conditions) != 2) {
    stop("This differential-expression pipeline requires exactly TWO conditions. Found: ",
         paste(conditions, collapse = ", "))
  }

  if (nzchar(trimws(controls))) {
    control_ids <- trimws(unlist(strsplit(controls, ",", fixed = TRUE)))
    control_ids <- control_ids[nzchar(control_ids)]

    control_idx <- match(control_ids, sample_ids)
    if (anyNA(control_idx)) {
      control_idx_norm <- match(normalize_id(control_ids), normalize_id(sample_ids))
      if (anyNA(control_idx_norm)) {
        missing_controls <- control_ids[is.na(control_idx_norm)]
        stop("Control sample IDs were not found in the count matrix: ", paste(missing_controls, collapse = ", "))
      }
      control_idx <- control_idx_norm
      cat("Used fallback control-ID matching with '.' converted to '-'.\n")
    }

    control_condition <- unique(md$condition[control_idx])
    if (length(control_condition) != 1) {
      stop("All samples listed in --controls must belong to the same metadata condition.")
    }

    other_condition <- setdiff(conditions, control_condition)
    conditions <- c(control_condition, other_condition)
  }

  md$condition <- factor(md$condition, levels = conditions)
  md
}

# ---------------------------------------------------------
# Load nf-core output
#
# The matrix orientation can vary:
#
#   1. mature miRNAs x samples
#      rows    = miRNAs
#      columns = samples
#
#   2. samples x mature miRNAs
#      rows    = samples
#      columns = miRNAs
#
# We determine the orientation using the user-provided
# metadata sample_id column rather than assuming it.
# ---------------------------------------------------------

counts_raw <- read.csv(
  opt$counts,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

logtpm_raw <- read.csv(
  opt$logtpm,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------
# Read metadata sample IDs for orientation detection
# ---------------------------------------------------------

metadata_raw <- read.csv(
  opt$metadata,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_metadata <- c("sample_id", "condition")

missing_metadata <- setdiff(
  required_metadata,
  colnames(metadata_raw)
)

if (length(missing_metadata) > 0) {
  stop(
    "Metadata must contain columns: sample_id, condition. Missing: ",
    paste(missing_metadata, collapse = ", ")
  )
}

metadata_sample_ids <- trimws(
  as.character(metadata_raw$sample_id)
)

metadata_sample_ids <- normalize_id(metadata_sample_ids)

# Normalize matrix names.
rownames(counts_raw) <- normalize_id(rownames(counts_raw))
colnames(counts_raw) <- normalize_id(colnames(counts_raw))

rownames(logtpm_raw) <- normalize_id(rownames(logtpm_raw))
colnames(logtpm_raw) <- normalize_id(colnames(logtpm_raw))

# ---------------------------------------------------------
# Helper: identify sample IDs in matrix names
# ---------------------------------------------------------
#
# First try exact matching.
# If that fails, allow the metadata sample ID to occur
# inside the matrix name.
#
# Example:
#   metadata: PVL02
#   matrix:   PVL02_S1
#
# This still identifies PVL02 as the corresponding sample.
# ---------------------------------------------------------

match_sample_names <- function(matrix_ids, metadata_ids) {

  result <- rep(NA_character_, length(matrix_ids))

  # Exact matches first
  exact <- match(matrix_ids, metadata_ids)
  result[!is.na(exact)] <- metadata_ids[exact[!is.na(exact)]]

  # Partial matching for remaining IDs
  for (i in which(is.na(result))) {

    hits <- metadata_ids[
      vapply(
        metadata_ids,
        function(id) {
          grepl(
            paste0("(^|[_\\.-])", id, "($|[_\\.-])"),
            matrix_ids[i],
            ignore.case = FALSE
          )
        },
        logical(1)
      )
    ]

    if (length(hits) == 1) {
      result[i] <- hits
    } else if (length(hits) > 1) {
      stop(
        "Ambiguous sample ID match for matrix ID '",
        matrix_ids[i],
        "'. Matches: ",
        paste(hits, collapse = ", ")
      )
    }
  }

  result
}

row_matches <- match_sample_names(
  rownames(counts_raw),
  metadata_sample_ids
)

column_matches <- match_sample_names(
  colnames(counts_raw),
  metadata_sample_ids
)

row_sample_matches <- sum(!is.na(row_matches))
column_sample_matches <- sum(!is.na(column_matches))

cat("\n--- Matrix orientation detection ---\n")
cat("Metadata sample IDs:", length(metadata_sample_ids), "\n")
cat("Sample IDs found in matrix rows:", row_sample_matches, "\n")
cat("Sample IDs found in matrix columns:", column_sample_matches, "\n")

# ---------------------------------------------------------
# Determine orientation
# ---------------------------------------------------------

if (
  column_sample_matches > row_sample_matches &&
  column_sample_matches >= 2
) {

  cat("Detected orientation: mature miRNAs x samples\n")
  cat("Transposing matrix: samples x mature miRNAs\n")

  # Keep only columns that correspond to metadata samples.
  selected_columns <- which(!is.na(column_matches))

  counts_selected <- counts_raw[
    ,
    selected_columns,
    drop = FALSE
  ]

  logtpm_selected <- logtpm_raw[
    ,
    selected_columns,
    drop = FALSE
  ]

  # Rename matrix sample columns to the metadata sample IDs.
  colnames(counts_selected) <- column_matches[selected_columns]
  colnames(logtpm_selected) <- column_matches[selected_columns]

  # Match counts/logTPM features.
  common_mirnas <- intersect(
    rownames(counts_selected),
    rownames(logtpm_selected)
  )

  if (length(common_mirnas) < 1) {
    stop("No common mature miRNAs between counts and logTPM matrices.")
  }

  counts_selected <- counts_selected[
    common_mirnas,
    ,
    drop = FALSE
  ]

  logtpm_selected <- logtpm_selected[
    common_mirnas,
    ,
    drop = FALSE
  ]

  # Transpose ONCE.
  counts <- as.data.frame(
    t(as.matrix(counts_selected)),
    check.names = FALSE
  )

  logtpm <- as.data.frame(
    t(as.matrix(logtpm_selected)),
    check.names = FALSE
  )

} else if (
  row_sample_matches > column_sample_matches &&
  row_sample_matches >= 2
) {

  cat("Detected orientation: samples x mature miRNAs\n")
  cat("No transposition required.\n")

  selected_rows <- which(!is.na(row_matches))

  counts_selected <- counts_raw[
    selected_rows,
    ,
    drop = FALSE
  ]

  logtpm_selected <- logtpm_raw[
    selected_rows,
    ,
    drop = FALSE
  ]

  # Rename matrix sample rows to metadata sample IDs.
  rownames(counts_selected) <- row_matches[selected_rows]
  rownames(logtpm_selected) <- row_matches[selected_rows]

  common_mirnas <- intersect(
    colnames(counts_selected),
    colnames(logtpm_selected)
  )

  if (length(common_mirnas) < 1) {
    stop("No common mature miRNAs between counts and logTPM matrices.")
  }

  counts <- counts_selected[
    ,
    common_mirnas,
    drop = FALSE
  ]

  logtpm <- logtpm_selected[
    ,
    common_mirnas,
    drop = FALSE
  ]

} else {

  cat("\nMetadata sample IDs:\n")
  print(metadata_sample_ids)

  cat("\nFirst matrix row IDs:\n")
  print(head(rownames(counts_raw), 10))

  cat("\nFirst matrix column IDs:\n")
  print(head(colnames(counts_raw), 10))

  stop(
    paste0(
      "\nCould not determine matrix orientation from metadata.\n",
      "Sample IDs found in rows: ", row_sample_matches, "\n",
      "Sample IDs found in columns: ", column_sample_matches, "\n"
    )
  )
}

# ---------------------------------------------------------
# Validate numeric content explicitly
# ---------------------------------------------------------

for (j in seq_len(ncol(counts))) {

  x <- suppressWarnings(
    as.numeric(as.character(counts[[j]]))
  )

  if (anyNA(x) && any(!is.na(counts[[j]]))) {
    stop(
      "Non-numeric values found in counts column: ",
      colnames(counts)[j]
    )
  }

  counts[[j]] <- x
}

for (j in seq_len(ncol(logtpm))) {

  x <- suppressWarnings(
    as.numeric(as.character(logtpm[[j]]))
  )

  if (anyNA(x) && any(!is.na(logtpm[[j]]))) {
    stop(
      "Non-numeric values found in logTPM column: ",
      colnames(logtpm)[j]
    )
  }

  logtpm[[j]] <- x
}

cat("\n--- Final matrix orientation ---\n")
cat("Samples:", nrow(counts), "\n")
cat("Mature miRNAs:", ncol(counts), "\n")

# ---------------------------------------------------------
# Metadata matching
# ---------------------------------------------------------
metadata_matched <- read_metadata(opt$metadata, rownames(counts), opt$controls)
group_factors <- metadata_matched$condition

# Ensure metadata order exactly follows the matrix row order.
metadata_matched <- metadata_matched[rownames(counts), , drop = FALSE]

# ---------------------------------------------------------
# edgeR group-wise filtering
# ---------------------------------------------------------
counts_matrix <- t(as.matrix(counts))
dge <- DGEList(counts = counts_matrix, group = group_factors)
design <- model.matrix(~ 0 + group_factors)

h <- hat(design)
MinSampleSize <- 1/max(h)

cat("Minimum sample size used by filterByExpr:", MinSampleSize, "\n")
cat("Median library size:", median(colSums(dge$counts)), "\n")

keep_miRNAs <- filterByExpr(dge, design = design)

filtered_counts <- counts[, keep_miRNAs, drop = FALSE]
filtered_logtpm <- logtpm[, keep_miRNAs, drop = FALSE]

cat("\n--- edgeR filterByExpr SUMMARY ---\n")
cat("Original number of mature miRNAs:", ncol(counts), "\n")
cat("Filtered number of mature miRNAs:", ncol(filtered_counts), "\n")
cat("Removed:", ncol(counts) - ncol(filtered_counts), "low-expressed mature miRNAs\n")
cat("\nDetected conditions:\n")
print(table(group_factors))

write.csv(filtered_logtpm,
          file.path(opt$outdir, "mature_logtpm_filtered.csv"),
          row.names = TRUE)
write.csv(filtered_counts,
          file.path(opt$outdir, "mature_counts_filtered.csv"),
          row.names = TRUE)
write.csv(metadata_matched,
          file.path(opt$outdir, "metadata.csv"),
          row.names = FALSE)
