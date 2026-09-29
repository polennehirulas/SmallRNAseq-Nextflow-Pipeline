
# Small RNA-seq Downstream Analysis Pipeline

A modular **Nextflow DSL2 pipeline** for downstream analysis of small RNA-seq mature miRNA expression data. The workflow performs filtering, differential expression, miRNA target analysis, functional enrichment, visualization, and result summarization.

## Workflow

```text
mature_counts.csv + mature_logtpm.csv + metadata.csv
                         │
                         ▼
                    FILTER_SMRNA
                         │
                         ▼
                  DIFF_EXP_SMRNA
                         │
                         ▼
                 SPLIT_MIRNAS_SMRNA
                    │          │
                    ▼          ▼
              Target genes   Enrichment
                    │          │
                    └────┬─────┘
                         ▼
                 VISUALIZATION
                         │
                         ▼
                  SUMMARY_COMPILE
                         │
                         ▼
                    CSV_TO_XLSX
```

##Environment

The repository includes `environment.yaml`, which specifies the Python, R, and Bioconductor dependencies used by the pipeline.

The workflow requires a compatible **Nextflow** installation and an environment containing the dependencies specified in this file.

## Inputs

### Expression matrices

| File                | Description                   |
| ------------------- | ----------------------------- |
| `mature_counts.csv` | Raw mature miRNA count matrix |
| `mature_logtpm.csv` | Log-transformed TPM matrix    |

Samples should correspond between the two matrices.

### Metadata

The pipeline uses a separate metadata file:

```text
input/metadata/metadata.csv
```

Required columns:

```text
sample_id,condition
```

Example:

```csv
sample_id,condition
SR_C1,Control
SR_C2,Control
SR_T1,Treated
SR_T2,Treated
```

Sample IDs are matched to the expression matrices, including automatic `.` → `-` normalization when necessary.

## Configuration

Main parameters are defined in `nextflow.config`:

```groovy
params {
    metadata = "${projectDir}/input/metadata/metadata.csv"
    controls = "SR_C1,SR_C2,SR_C3"

    de_significance = "fdr"
    de_cutoff = 0.05

    enrichment_significance = "qvalue"
    enrichment_cutoff = 0.05

    smrna_counts = "${projectDir}/mature_counts.csv"
    smrna_logtpm = "${projectDir}/mature_logtpm.csv"

    outdir = "${projectDir}/smRNA_results"
}
```

* `controls`: comma-separated control sample IDs.
* `de_significance`: `pvalue` or `fdr`.
* `de_cutoff`: significance threshold for differential expression and significant-miRNA splitting.
* `enrichment_significance`: enrichment significance metric.
* `enrichment_cutoff`: enrichment threshold.
* `outdir`: final results directory.

## Analysis Steps

| Process               | Script                                                        | Main output                        |
| --------------------- | ------------------------------------------------------------- | ---------------------------------- |
| `FILTER_SMRNA`        | `filter_smrna.R`                                              | Filtered counts/logTPM             |
| `DIFF_EXP_SMRNA`      | `diff_exp_smrna.R`                                            | DE results, `dge.rds`, BCV plot    |
| `SPLIT_MIRNAS_SMRNA`  | `split_mirnas_smrna.R`                                        | Up, Down, and combined miRNA lists |
| `TARGET_GENE_SMRNA`   | `target_gene_smrna.R`                                         | Validated miRNA target genes       |
| `ENRICHMENT_SMRNA`    | `enrichment_smrna.R`                                          | GO BP/CC/MF and KEGG results/plots |
| `VISUALIZATION_SMRNA` | `PCA_plot_smrna.R`, `heatmap_smrna.R`, `volcano_plot_smrna.R` | PCA, heatmap, volcano plots        |
| `SUMMARY_COMPILE`     | —                                                             | Combined `summary.txt`             |
| `CSV_TO_XLSX`         | —                                                             | Excel versions of CSV results      |

### miRNA branches

Target and enrichment analysis are performed separately for:

```text
Up_Regulated/
Down_Regulated/
Up_and_Down/
```

The significant miRNA threshold is inherited from the DE configuration.

## Main Outputs

```text
smRNA_results/
├── filtered/
├── dge/
├── split/
│   ├── Up_Regulated/
│   ├── Down_Regulated/
│   └── Up_and_Down/
├── targets/
├── enrichment/
├── visualization/
├── summary/
└── xlsx/
```

Important results include:

* filtered mature miRNA matrices
* `differential_expression_results.csv`
* `dge.rds`
* `BCV_plot.png`
* significant Up/Down/combined miRNA lists
* validated miRNA target tables
* unique target gene lists
* GO Biological Process, Cellular Component and Molecular Function results
* KEGG enrichment results
* PCA, heatmap and volcano plots
* `summary.txt`
* Excel versions of CSV results

## Run

From the pipeline directory:

```bash
nextflow run main.nf
```

Resume an interrupted or completed run:

```bash
nextflow run main.nf -resume
```

The pipeline is controlled through `nextflow.config`, so metadata, controls, statistical thresholds, input files and output location can be changed without modifying the R scripts.

