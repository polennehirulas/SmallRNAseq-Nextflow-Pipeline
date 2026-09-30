#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))

option_list <- list(
  make_option(c("--targets"), type = "character", help = "CSV containing target_entrez"),
  make_option(c("--outdir"), type = "character", default = "."),
  make_option(c("--direction"), type = "character", default = "All_Significant"),
  make_option(c("--ont"), type = "character", default = NULL, help = "GO ontology: BP, CC or MF"),
  make_option(c("--kegg"), action = "store_true", default = FALSE),
  make_option(c("--significance"), type = "character", default = "pvalue",
            help = "Enrichment significance metric: pvalue or qvalue"),
  make_option(c("--cutoff"), type = "double", default = 0.05,
            help = "Enrichment significance cutoff")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (!opt$significance %in% c("pvalue", "qvalue")) {
  stop("--significance must be either 'pvalue' or 'qvalue'.")
}

if (opt$cutoff <= 0 || opt$cutoff >= 1) {
  stop("--cutoff must be greater than 0 and less than 1.")
}
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(opt$targets)) stop("Target gene file not found: ", opt$targets)

target_genes_df <- read.csv(opt$targets, check.names = FALSE, stringsAsFactors = FALSE)
if (!"target_entrez" %in% colnames(target_genes_df)) stop("Input file must contain target_entrez column.")

entrez_ids <- unique(as.character(target_genes_df$target_entrez))
entrez_ids <- entrez_ids[!is.na(entrez_ids) & nzchar(trimws(entrez_ids))]

save_blank_plot <- function(path, title_text) {
  png(path, width = 2700, height = 2100, res = 300)
  plot.new()
  title(main = title_text)
  dev.off()
}

cat("Context profile:", opt$direction, "\n")
cat("Enrichment significance metric:", opt$significance, "\n")
cat("Enrichment significance cutoff:", opt$cutoff, "\n")
cat("Loaded", length(entrez_ids), "valid Entrez IDs.\n")

if (!is.null(opt$ont)) {
  ont <- toupper(opt$ont)
  if (!ont %in% c("BP", "CC", "MF")) stop("--ont must be BP, CC or MF.")

  ontology_names <- c(BP = "Biological Process", CC = "Cellular Component", MF = "Molecular Function")
  csv_name <- paste0("GO_", ont, "_enrichment_results.csv")
  plot_name <- paste0("GO_dotplot_", ont, ".png")
  title_text <- paste("GO", ontology_names[[ont]], "Enrichment")

  if (length(entrez_ids) == 0) {
    write.csv(data.frame(), file.path(opt$outdir, csv_name), row.names = FALSE)
    save_blank_plot(file.path(opt$outdir, plot_name), paste(title_text, "- No target genes available"))
    quit(save = "no", status = 0)
  }

  enrichment_args <- list(
  gene = entrez_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = ont,
  pAdjustMethod = "BH",
  readable = TRUE
  )
  
  if (opt$significance == "pvalue") {
    enrichment_args$pvalueCutoff <- opt$cutoff
    enrichment_args$qvalueCutoff <- 1
  } else {
    enrichment_args$pvalueCutoff <- 1
    enrichment_args$qvalueCutoff <- opt$cutoff
  }

  ego <- do.call(enrichGO, enrichment_args)
  

  ego_df <- if (!is.null(ego)) as.data.frame(ego) else data.frame()
  write.csv(ego_df, file.path(opt$outdir, csv_name), row.names = FALSE)

  if (nrow(ego_df) > 0) {

    go_plot <- dotplot(
      ego,
      showCategory = 20,
      title = title_text
    ) +
      scale_y_discrete(
        labels = function(x) stringr::str_wrap(x, width = 50)
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold")
      )

    n_terms <- min(20, nrow(ego_df))
    plot_height <- max(4.5, min(8, nrow(ego_df) * 0.30))

    ggsave(
      file.path(opt$outdir, plot_name),
      go_plot,
      width = 9,
      height = plot_height,
      dpi = 300
    )

  } else {

    save_blank_plot(
      file.path(opt$outdir, plot_name),
      paste0(
        "GO ",
        ont,
        " Enrichment - No significant terms"
      )
    )
  }
}


# ============================================================
# KEGG ENRICHMENT
# ============================================================

if (isTRUE(opt$kegg)) {

  csv_name <- "KEGG_enrichment_results.csv"
  plot_name <- "KEGG_dotplot.png"

  if (length(entrez_ids) == 0) {

    write.csv(
      data.frame(),
      file.path(opt$outdir, csv_name),
      row.names = FALSE
    )

    save_blank_plot(
      file.path(opt$outdir, plot_name),
      "KEGG Pathway Enrichment - No target genes available"
    )

    quit(save = "no", status = 0)
  }

  enrichment_args <- list(
    gene = entrez_ids,
    organism = "hsa",
    pAdjustMethod = "BH"
  )

  if (opt$significance == "pvalue") {

    enrichment_args$pvalueCutoff <- opt$cutoff
    enrichment_args$qvalueCutoff <- 1

  } else {

    enrichment_args$pvalueCutoff <- 1
    enrichment_args$qvalueCutoff <- opt$cutoff
  }

  ekegg <- do.call(enrichKEGG, enrichment_args)

  ekegg_df <- if (!is.null(ekegg)) {
    as.data.frame(ekegg)
  } else {
    data.frame()
  }

  write.csv(
    ekegg_df,
    file.path(opt$outdir, csv_name),
    row.names = FALSE
  )

  if (nrow(ekegg_df) > 0) {

    kegg_plot <- dotplot(
      ekegg,
      showCategory = 20,
      title = "KEGG Pathway Enrichment"
    ) +
      scale_y_discrete(
        labels = function(x) stringr::str_wrap(x, width = 50)
      ) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold")
      )

    n_terms <- min(20, nrow(ekegg_df))
    plot_height <- max(4.5, min(8, nrow(ekegg_df) * 0.30))

    ggsave(
      file.path(opt$outdir, plot_name),
      kegg_plot,
      width = 9,
      height = plot_height,
      dpi = 300
    )

  } else {

    save_blank_plot(
      file.path(opt$outdir, plot_name),
      "KEGG Pathway Enrichment - No significant pathways"
    )
  }
}
