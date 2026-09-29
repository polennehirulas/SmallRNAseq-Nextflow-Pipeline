#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))

option_list <- list(
  make_option(
    c("--dge-rds"),
    type = "character",
    dest = "dge_rds",
    help = "edgeR DGEList RDS"
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
  )
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(opt$dge_rds)) stop("dge.rds not found.")
if (!file.exists(opt$metadata)) stop("Metadata file not found.")

dge <- readRDS(opt$dge_rds)
md <- read.csv(opt$metadata, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("sample_id", "condition") %in% colnames(md))) stop("Metadata must contain sample_id and condition.")

md$sample_id <- trimws(as.character(md$sample_id))
md$condition <- trimws(as.character(md$condition))

sample_ids <- colnames(dge$counts)
idx <- match(sample_ids, md$sample_id)
if (anyNA(idx)) {
  idx <- match(chartr(".", "-", sample_ids), chartr(".", "-", md$sample_id))
}
if (anyNA(idx)) stop("Some DGE samples are missing from metadata: ", paste(sample_ids[is.na(idx)], collapse = ", "))

md <- md[idx, , drop = FALSE]
rownames(md) <- sample_ids

group_factors <- factor(md$condition, levels = unique(md$condition))

logCPM_matrix <- cpm(dge, log = TRUE)
pca_input <- t(logCPM_matrix)
pca_results <- prcomp(pca_input, center = TRUE, scale. = TRUE)
var_explained <- pca_results$sdev^2 / sum(pca_results$sdev^2) * 100

pca_df <- data.frame(
  SampleID = rownames(pca_input),
  PC1 = pca_results$x[, 1],
  PC2 = pca_results$x[, 2],
  Group = group_factors
)

p <- ggplot(pca_df, aes(PC1, PC2, color = Group)) +
  geom_point(size = 4.5, alpha = 0.85) +
  geom_text_repel(
    aes(label = SampleID),
    size = 3.2,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.linewidth = 0.4,
    min.segment.length = 0,
    show.legend = FALSE
  ) +
  labs(
    title = "PCA of normalized miRNA expression",
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
    color = "Group"
  ) +
  theme_classic(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    legend.background = element_rect(fill = "white", colour = NA),
    plot.title = element_text(face = "plain", size = 16, margin = margin(b = 15)),
    legend.position = "right",
    legend.title = element_text(face = "plain", size = 13),
    legend.text = element_text(size = 12)
   )

ggsave(file.path(opt$outdir, "PCA_normalized_miRNA_expression.png"), p, width = 10, height = 8, dpi = 300)
