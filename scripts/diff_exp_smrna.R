#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(edgeR))

action_safe_numeric <- function(x) {
  out <- suppressWarnings(as.numeric(as.character(x)))
  if (anyNA(out) && any(!is.na(x))) stop("Non-numeric values detected in counts matrix.")
  out
}

normalize_id <- function(x) chartr(".", "-", trimws(as.character(x)))

read_metadata <- function(file, sample_ids, controls = "") {
  md <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)

  required <- c("sample_id", "condition")
  missing <- setdiff(required, colnames(md))
  if (length(missing) > 0) stop("Metadata must contain sample_id and condition. Missing: ", paste(missing, collapse = ", "))

  md$sample_id <- trimws(as.character(md$sample_id))
  md$condition <- trimws(as.character(md$condition))
  if (anyDuplicated(md$sample_id)) stop("Metadata contains duplicated sample_id values.")

  idx <- match(sample_ids, md$sample_id)
  if (anyNA(idx)) {
    md_norm <- normalize_id(md$sample_id)
    ids_norm <- normalize_id(sample_ids)
    if (anyDuplicated(md_norm)) stop("Metadata sample IDs collide after '.' -> '-' normalization.")
    idx_norm <- match(ids_norm, md_norm)
    if (anyNA(idx_norm)) {
      stop("Count-matrix samples missing from metadata: ", paste(sample_ids[is.na(idx_norm)], collapse = ", "))
    }
    idx <- idx_norm
    cat("Used fallback sample-ID matching with '.' converted to '-'.\n")
  }

  md <- md[idx, , drop = FALSE]
  rownames(md) <- sample_ids

  conditions <- unique(md$condition)
  if (length(conditions) != 2) {
    stop("Differential expression requires exactly two conditions. Found: ", paste(conditions, collapse = ", "))
  }

  if (nzchar(trimws(controls))) {
    control_ids <- trimws(unlist(strsplit(controls, ",", fixed = TRUE)))
    control_ids <- control_ids[nzchar(control_ids)]
    ctrl_idx <- match(control_ids, sample_ids)
    if (anyNA(ctrl_idx)) {
      ctrl_idx_norm <- match(normalize_id(control_ids), normalize_id(sample_ids))
      if (anyNA(ctrl_idx_norm)) stop("Control IDs not found: ", paste(control_ids[is.na(ctrl_idx_norm)], collapse = ", "))
      ctrl_idx <- ctrl_idx_norm
      cat("Used fallback control-ID matching with '.' converted to '-'.\n")
    }
    control_condition <- unique(md$condition[ctrl_idx])
    if (length(control_condition) != 1) stop("All --controls samples must share one metadata condition.")
    conditions <- c(control_condition, setdiff(conditions, control_condition))
  }

  md$condition <- factor(md$condition, levels = conditions)
  md
}

option_list <- list(
  make_option(c("--counts"), type = "character", help = "Filtered mature counts CSV"),
  make_option(c("--metadata"), type = "character", help = "Metadata CSV"),
  make_option(c("--controls"), type = "character", default = "", help = "Optional comma-separated control sample IDs"),
  make_option(c("--significance"), type = "character", default = "pvalue",
            help = "Significance metric: pvalue or fdr"),
  make_option(c("--cutoff"), type = "double", default = 0.05,
            help = "Significance cutoff"),
  make_option(c("--outdir"), type = "character", default = ".")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

opt <- parse_args(OptionParser(option_list = option_list))

if (!opt$significance %in% c("pvalue", "fdr")) {
  stop("--significance must be either 'pvalue' or 'fdr'.")
}

if (opt$cutoff <= 0 || opt$cutoff >= 1) {
  stop("--cutoff must be greater than 0 and less than 1.")
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(opt$counts)) stop("Filtered counts file not found.")
if (!file.exists(opt$metadata)) stop("Metadata file not found.")

filtered_counts <- read.csv(opt$counts, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(filtered_counts) < 2 || ncol(filtered_counts) < 1) stop("Filtered counts matrix is empty or too small.")
if (anyDuplicated(rownames(filtered_counts))) stop("Filtered counts contains duplicated sample IDs.")

for (j in seq_len(ncol(filtered_counts))) filtered_counts[[j]] <- action_safe_numeric(filtered_counts[[j]])

metadata_matched <- read_metadata(opt$metadata, rownames(filtered_counts), opt$controls)
group_factors <- metadata_matched$condition

cat("Detected experimental groups:\n")
print(table(group_factors))
cat("Control/reference condition:", levels(group_factors)[1], "\n")
cat("Comparison condition:", levels(group_factors)[2], "\n")

# Samples are rows in our downstream tables; edgeR requires
# features (mature miRNAs) in rows and samples in columns.
counts_matrix <- t(as.matrix(filtered_counts))
dge <- DGEList(counts = counts_matrix, group = group_factors)

dge <- calcNormFactors(dge, method = "TMM")
dge <- estimateDisp(dge)

grp_levels <- levels(group_factors)
et <- exactTest(dge, pair = grp_levels[c(1, 2)])

top_tags <- topTags(et, n = Inf)$table

significance_values <- if (opt$significance == "pvalue") {
  top_tags$PValue
} else {
  top_tags$FDR
}

significant_miRNAs <- top_tags[
  !is.na(significance_values) &
  significance_values < opt$cutoff,
  ,
  drop = FALSE
]

cat("\n--- DIFFERENTIAL EXPRESSION SUMMARY ---\n")
cat("Significance metric:", opt$significance, "\n")
cat("Significance cutoff:", opt$cutoff, "\n")
cat("Number of significant mature miRNAs:", nrow(significant_miRNAs), "\n")

write.csv(top_tags,
          file.path(opt$outdir, "differential_expression_results.csv"),
          row.names = TRUE)

saveRDS(dge, file.path(opt$outdir, "dge.rds"))

png(file.path(opt$outdir, "BCV_plot.png"), width = 2400, height = 1800, res = 300)
plotBCV(dge, main = "Biological Coefficient of Variation")
dev.off()
