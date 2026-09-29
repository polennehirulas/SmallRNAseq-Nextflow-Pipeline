#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(pheatmap))

option_list <- list(
  make_option(
    c("--dge-rds"),
    type = "character",
    dest = "dge_rds",
    help = "edgeR DGEList RDS"
  ),
  make_option(
    c("--dge"),
    type = "character",
    dest = "dge",
    help = "Differential expression results CSV"
  ),
  make_option(
    c("--metadata"),
    type = "character",
    dest = "metadata",
    help = "Metadata CSV"
  ),
  make_option(
    c("--outdir"),
    type = "character",
    dest = "outdir",
    default = "."
  ),
  make_option(
  c("--significance"),
  type = "character",
  default = "pvalue",
  help = "Significance metric: pvalue or fdr"
  ),
  make_option(
  c("--cutoff"),
  type = "double",
  default = 0.05,
  help = "Significance cutoff"
 )
)
opt <- parse_args(OptionParser(option_list = option_list))
if (!opt$significance %in% c("pvalue", "fdr")) {
  stop("--significance must be either 'pvalue' or 'fdr'.")
}

if (opt$cutoff <= 0 || opt$cutoff >= 1) {
  stop("--cutoff must be greater than 0 and less than 1.")
}
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(opt$dge_rds)) stop("dge.rds not found.")
if (!file.exists(opt$dge)) stop("DE results file not found.")
if (!file.exists(opt$metadata)) stop("Metadata file not found.")

dge <- readRDS(opt$dge_rds)
de <- read.csv(opt$dge, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
md <- read.csv(opt$metadata, check.names = FALSE, stringsAsFactors = FALSE)

if (!all(c("sample_id", "condition") %in% colnames(md))) stop("Metadata must contain sample_id and condition.")
md$sample_id <- trimws(as.character(md$sample_id))
md$condition <- trimws(as.character(md$condition))

sample_ids <- colnames(dge$counts)
idx <- match(sample_ids, md$sample_id)
if (anyNA(idx)) idx <- match(chartr(".", "-", sample_ids), chartr(".", "-", md$sample_id))
if (anyNA(idx)) stop("Some DGE samples are missing from metadata.")
md <- md[idx, , drop = FALSE]
rownames(md) <- sample_ids

significance_values <- if (opt$significance == "pvalue") {
  de$PValue
} else {
  de$FDR
}

sig_mirnas <- rownames(
  de[
    !is.na(significance_values) &
    significance_values < opt$cutoff,
    ,
    drop = FALSE
  ]
)
logCPM_matrix <- cpm(dge, log = TRUE)
valid_mirnas <- sig_mirnas[sig_mirnas %in% rownames(logCPM_matrix)]
heatmap_matrix <- logCPM_matrix[valid_mirnas, , drop = FALSE]

out_file <- file.path(opt$outdir, "Significant_miRNAs_Heatmap.png")

if (nrow(heatmap_matrix) == 0) {
  png(out_file, width = 3000, height = 2400, res = 300)
  plot.new()
  title(main = "Expression Profiles of Significant Mature miRNAs - No significant miRNAs")
  dev.off()
  quit(save = "no", status = 0)
}

# Avoid undefined z-scores for rows with zero variance.
row_sd <- apply(heatmap_matrix, 1, sd)
heatmap_matrix <- heatmap_matrix[row_sd > 0 & !is.na(row_sd), , drop = FALSE]

if (nrow(heatmap_matrix) == 0) {
  png(out_file, width = 3000, height = 2400, res = 300)
  plot.new()
  title(main = "Expression Profiles of Significant Mature miRNAs - No variable miRNAs")
  dev.off()
  quit(save = "no", status = 0)
}

annotation_columns <- data.frame(Group = md$condition)
rownames(annotation_columns) <- colnames(heatmap_matrix)

ordered_samples <- order(annotation_columns$Group)
heatmap_matrix <- heatmap_matrix[, ordered_samples, drop = FALSE]
annotation_columns <- annotation_columns[ordered_samples, , drop = FALSE]

pheatmap(
  mat = heatmap_matrix,
  annotation_col = annotation_columns,
  cluster_rows = nrow(heatmap_matrix) > 1,
  cluster_cols = FALSE,
  scale = "row",
  show_colnames = TRUE,
  show_rownames = TRUE,
  color = colorRampPalette(c("#313695", "#f4f6f7", "#a50026"))(100),
  border_color = "white",
  main = "Expression Profiles of Significant Mature miRNAs",
  filename = out_file,
  width = 10,
  height = 8
)
