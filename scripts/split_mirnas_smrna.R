#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("--dge"), type = "character",
              help = "Differential expression results CSV"),
  make_option(c("--outdir"), type = "character", default = "."),
  make_option(c("--significance"), type = "character", default = "pvalue",
            help = "Significance metric: pvalue or fdr"),
  make_option(c("--cutoff"), type = "double", default = 0.05,
            help = "Significance cutoff")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!opt$significance %in% c("pvalue", "fdr")) {
  stop("--significance must be either 'pvalue' or 'fdr'.")
}

if (opt$cutoff <= 0 || opt$cutoff >= 1) {
  stop("--cutoff must be greater than 0 and less than 1.")
}

if (!file.exists(opt$dge)) {
  stop("DE results file not found: ", opt$dge)
}

# ---------------------------------------------------------
# 1. Create output directories
# ---------------------------------------------------------

dir.create(
  file.path(opt$outdir, "Up_Regulated"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(opt$outdir, "Down_Regulated"),
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  file.path(opt$outdir, "Up_and_Down"),
  recursive = TRUE,
  showWarnings = FALSE
)

# ---------------------------------------------------------
# 2. Read differential expression results
# ---------------------------------------------------------

de <- read.csv(
  opt$dge,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!all(c("logFC") %in% colnames(de))) {
  stop("DE results must contain logFC column.")
}

if (opt$significance == "pvalue" && !"PValue" %in% colnames(de)) {
  stop("PValue column is required when --significance pvalue is selected.")
}

if (opt$significance == "fdr" && !"FDR" %in% colnames(de)) {
  stop("FDR column is required when --significance fdr is selected.")
}

# ---------------------------------------------------------
# 3. Select significant mature miRNAs
# ---------------------------------------------------------

significance_values <- if (opt$significance == "pvalue") {
  de$PValue
} else {
  de$FDR
}

sig <- de[
  !is.na(significance_values) &
  significance_values < opt$cutoff &
  !is.na(de$logFC),
  ,
  drop = FALSE
]

# ---------------------------------------------------------
# 4. Separate Up / Down
# ---------------------------------------------------------

up_mirnas <- rownames(
  sig[sig$logFC > 0, , drop = FALSE]
)

down_mirnas <- rownames(
  sig[sig$logFC < 0, , drop = FALSE]
)

# ---------------------------------------------------------
# 5. Combined significant miRNA list
# ---------------------------------------------------------

combined_mirnas <- unique(
  c(up_mirnas, down_mirnas)
)

# ---------------------------------------------------------
# 6. Write input files
# ---------------------------------------------------------

writeLines(
  up_mirnas,
  file.path(opt$outdir, "Up_Regulated", "input_mirnas.txt"),
  useBytes = TRUE
)

writeLines(
  down_mirnas,
  file.path(opt$outdir, "Down_Regulated", "input_mirnas.txt"),
  useBytes = TRUE
)

writeLines(
  combined_mirnas,
  file.path(opt$outdir, "Up_and_Down", "input_mirnas.txt"),
  useBytes = TRUE
)

# ---------------------------------------------------------
# 7. Report summary
# ---------------------------------------------------------

cat("Significance metric:", opt$significance, "\n")
cat("Significance cutoff:", opt$cutoff, "\n")
cat("Significant mature miRNAs:", nrow(sig), "\n")
cat("Up-regulated:", length(up_mirnas), "\n")
cat("Down-regulated:", length(down_mirnas), "\n")
cat("Up + Down combined:", length(combined_mirnas), "\n")
