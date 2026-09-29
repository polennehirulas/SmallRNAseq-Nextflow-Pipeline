#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(EnhancedVolcano))
suppressPackageStartupMessages(library(ggplot2))


# =========================================================
# 1. COMMAND-LINE ARGUMENTS
# =========================================================

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
    c("--controls"),
    type = "character",
    dest = "controls",
    default = "",
    help = "Comma-separated control sample IDs"
  ),
  make_option(
    c("--significance"),
    type = "character",
    dest = "significance",
    default = "pvalue",
    help = "Significance metric: pvalue or fdr"
  ),
  make_option(
    c("--cutoff"),
    type = "double",
    dest = "cutoff",
    default = 0.05,
    help = "Significance cutoff"
  ),
  make_option(
    c("--outdir"),
    type = "character",
    dest = "outdir",
    default = "."
  )
)

opt <- parse_args(OptionParser(option_list = option_list))
if (!opt$significance %in% c("pvalue", "fdr")) {
  stop("--significance must be either 'pvalue' or 'fdr'.")
}

if (opt$cutoff <= 0 || opt$cutoff >= 1) {
  stop("--cutoff must be greater than 0 and less than 1.")
}
dir.create(
  opt$outdir,
  recursive = TRUE,
  showWarnings = FALSE
)


# =========================================================
# 2. CHECK INPUT FILES
# =========================================================

if (!file.exists(opt$dge)) {
  stop("DGE results file not found: ", opt$dge)
}

if (!file.exists(opt$metadata)) {
  stop("Metadata file not found: ", opt$metadata)
}


# =========================================================
# 3. READ METADATA
# =========================================================

md <- read.csv(
  opt$metadata,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!all(c("sample_id", "condition") %in% colnames(md))) {
  stop(
    "Metadata must contain 'sample_id' and 'condition' columns."
  )
}

md$sample_id <- trimws(as.character(md$sample_id))
md$condition <- trimws(as.character(md$condition))

conditions <- unique(md$condition)

if (length(conditions) != 2) {
  stop(
    "Volcano plot requires exactly 2 conditions in metadata. Found: ",
    paste(conditions, collapse = ", ")
  )
}

# ---------------------------------------------------------
# Identify control/reference condition using --controls
# ---------------------------------------------------------

if (!nzchar(trimws(opt$controls))) {
  stop(
    "Control sample IDs must be provided with --controls so ",
    "the volcano plot can correctly define Up/Down."
  )
}

control_ids <- trimws(
  unlist(strsplit(opt$controls, ",", fixed = TRUE))
)

control_ids <- control_ids[nzchar(control_ids)]

# Match control sample IDs to metadata
control_idx <- match(control_ids, md$sample_id)

if (anyNA(control_idx)) {
  control_idx <- match(
    chartr(".", "-", control_ids),
    chartr(".", "-", md$sample_id)
  )
}

if (anyNA(control_idx)) {
  stop(
    "Control sample IDs not found in metadata: ",
    paste(control_ids[is.na(control_idx)], collapse = ", ")
  )
}

control_condition <- unique(
  md$condition[control_idx]
)

if (length(control_condition) != 1) {
  stop(
    "All control samples must belong to the same condition."
  )
}

comparison_condition <- setdiff(
  conditions,
  control_condition
)

if (length(comparison_condition) != 1) {
  stop(
    "Could not identify exactly one comparison condition."
  )
}

cat("Control/reference condition:", control_condition, "\n")
cat("Comparison/patient condition:", comparison_condition, "\n")


# =========================================================
# 4. READ DIFFERENTIAL EXPRESSION RESULTS
# =========================================================

de_data <- read.csv(
  opt$dge,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!"logFC" %in% colnames(de_data)) {
  stop("DE results must contain 'logFC' column.")
}

if (opt$significance == "pvalue" && !"PValue" %in% colnames(de_data)) {
  stop("PValue column is required when --significance pvalue is selected.")
}

if (opt$significance == "fdr" && !"FDR" %in% colnames(de_data)) {
  stop("FDR column is required when --significance fdr is selected.")
}


# =========================================================
# 5. CREATE DYNAMIC VOLCANO CATEGORIES
# =========================================================

up_label <- paste(
  "Up in",
  comparison_condition
)

down_label <- paste(
  "Down in",
  comparison_condition
)

significance_values <- if (opt$significance == "pvalue") {
  de_data$PValue
} else {
  de_data$FDR
}

keyvals <- rep(
  "Not significant",
  nrow(de_data)
)

keyvals[
  !is.na(significance_values) &
  significance_values < opt$cutoff &
  de_data$logFC > 0
] <- up_label

keyvals[
  !is.na(significance_values) &
  significance_values < opt$cutoff &
  de_data$logFC < 0
] <- down_label

names(keyvals) <- keyvals


# =========================================================
# 6. CUSTOM COLORS
# =========================================================

col_custom <- c(
  setNames("#d1493a", up_label),
  setNames("#3a9fd1", down_label),
  "Not significant" = "#bebebe"
)

# =========================================================
# 7. LABEL SIGNIFICANT miRNAs
# =========================================================

select_lab <- rownames(de_data)[
  !is.na(significance_values) &
  significance_values < opt$cutoff
]


# =========================================================
# 8. CREATE VOLCANO PLOT
# =========================================================
y_variable <- if (opt$significance == "pvalue") {
  "PValue"
} else {
  "FDR"
}

v_plot <- EnhancedVolcano(
  toptable = de_data,
  lab = rownames(de_data),
  x = "logFC",
  y = y_variable,,

  pCutoff = opt$cutoff,
  FCcutoff = 0,

  pointSize = 3.5,
  labSize = 4.5,

  colCustom = col_custom[names(keyvals)],

  colAlpha = 0.8,

  selectLab = select_lab,

  legendPosition = "right",
  legendLabSize = 11,
  legendIconSize = 4.0,

  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = "grey30",

  gridlines.major = TRUE,
  gridlines.minor = FALSE,

  border = "partial",

  title = "Volcano Plot of Differentially Expressed Mature miRNAs",

  subtitle = paste(
    comparison_condition,
    "vs",
    control_condition,
    "|",
    toupper(opt$significance),
    "<",
    opt$cutoff
),

  caption = NULL
) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 14
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 11
    )
  )


# =========================================================
# 9. SAVE PLOT
# =========================================================

ggsave(
  file.path(
    opt$outdir,
    "Volcano_mature_miRNAs.png"
  ),
  v_plot,
  width = 9,
  height = 8,
  dpi = 300
)
