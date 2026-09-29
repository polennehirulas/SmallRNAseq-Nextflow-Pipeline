#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(multiMiR))

option_list <- list(
  make_option(c("--mirnas"), type = "character", help = "TXT file containing mature miRNA IDs"),
  make_option(c("--outdir"), type = "character", default = ".")
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(opt$mirnas)) stop("miRNA list file not found: ", opt$mirnas)

mirna_list <- readLines(opt$mirnas, warn = FALSE)
mirna_list <- trimws(mirna_list)
mirna_list <- unique(mirna_list[nzchar(mirna_list)])

empty_target_outputs <- function() {
  write.csv(data.frame(), file.path(opt$outdir, "predicted_targets.csv"), row.names = FALSE)
  write.csv(data.frame(), file.path(opt$outdir, "mirna_validated_targets_details.csv"), row.names = FALSE)
  write.csv(data.frame(target_symbol = character(0), target_entrez = character(0)),
            file.path(opt$outdir, "unique_target_genes_list.csv"), row.names = FALSE)
}

if (length(mirna_list) == 0) {
  cat("No significant mature miRNAs were supplied. Creating empty target outputs.\n")
  empty_target_outputs()
  quit(save = "no", status = 0)
}

cat("Loaded", length(mirna_list), "mature miRNAs. Querying validated multiMiR targets...\n")

multimir_results <- get_multimir(
  org     = "hsa",
  mirna   = mirna_list,
  table   = "validated",
  summary = TRUE
)

target_table <- multimir_results@data

if (is.null(target_table) || nrow(target_table) == 0) {
  cat("No validated target interactions were found.\n")
  empty_target_outputs()
  quit(save = "no", status = 0)
}

required_cols <- c("mature_mirna_id", "target_symbol", "target_entrez", "database", "experiment")
missing_cols <- setdiff(required_cols, colnames(target_table))
if (length(missing_cols) > 0) {
  stop("multiMiR result is missing expected columns: ", paste(missing_cols, collapse = ", "))
}

report_targets <- target_table[, required_cols, drop = FALSE]
unique_genes_df <- unique(target_table[, c("target_symbol", "target_entrez"), drop = FALSE])

write.csv(target_table,
          file.path(opt$outdir, "predicted_targets.csv"),
          row.names = FALSE)
write.csv(report_targets,
          file.path(opt$outdir, "mirna_validated_targets_details.csv"),
          row.names = FALSE)
write.csv(unique_genes_df,
          file.path(opt$outdir, "unique_target_genes_list.csv"),
          row.names = FALSE)

cat("Total validated interactions:", nrow(target_table), "\n")
cat("Total unique target genes:", nrow(unique_genes_df), "\n")
