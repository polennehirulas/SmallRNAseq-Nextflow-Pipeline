// main.nf - small RNA-seq downstream analysis
// Includes stage-level summaries and a final pipeline_summary.txt

params.metadata = "${projectDir}/input/metadata/metadata.csv"

// =========================================================
// 1. FILTER MATURE miRNAs
// =========================================================
process FILTER_SMRNA {
    tag "filter mature miRNAs"
    publishDir "${params.outdir}/filtered", mode: 'copy'

    input:
    path counts
    path logtpm
    path metadata
    val controls

    output:
    path "mature_counts_filtered.csv", emit: counts_filtered
    path "mature_logtpm_filtered.csv", emit: logtpm_filtered
    path "metadata.csv", emit: metadata_copy
    path "filter_summary.txt", emit: summary

    script:
    """
    echo "============================================================" > filter_summary.txt
    echo "FILTER_SMRNA" >> filter_summary.txt
    echo "============================================================" >> filter_summary.txt
    echo "Input counts: ${counts}" >> filter_summary.txt
    echo "Input logTPM: ${logtpm}" >> filter_summary.txt
    echo "Metadata: ${metadata}" >> filter_summary.txt
    echo "Controls/reference samples: ${controls ?: '<not specified>'}" >> filter_summary.txt
    echo "" >> filter_summary.txt

    Rscript ${projectDir}/scripts/filter_smrna.R \\
        --counts ${counts} \\
        --logtpm ${logtpm} \\
        --metadata ${metadata} \\
        --controls "${controls ?: ''}" \\
        --outdir . >> filter_summary.txt
    """
}

// =========================================================
// 2. DIFFERENTIAL EXPRESSION
// =========================================================
process DIFF_EXP_SMRNA {
    tag "edgeR differential expression"
    publishDir "${params.outdir}/dge", mode: 'copy'

    input:
    path counts
    path metadata
    val controls

    output:
    path "differential_expression_results.csv", emit: dge_results
    path "dge.rds", emit: dge_rds
    path "BCV_plot.png", emit: bcv_plot
    path "dge_summary.txt", emit: summary

    script:
    """
    echo "============================================================" > dge_summary.txt
    echo "DIFF_EXP_SMRNA" >> dge_summary.txt
    echo "============================================================" >> dge_summary.txt
    echo "Filtered counts: ${counts}" >> dge_summary.txt
    echo "Metadata: ${metadata}" >> dge_summary.txt
    echo "Controls/reference samples: ${controls ?: '<not specified>'}" >> dge_summary.txt
    echo "" >> dge_summary.txt

    Rscript ${projectDir}/scripts/diff_exp_smrna.R \
        --counts ${counts} \
        --metadata ${metadata} \
        --controls "${controls ?: ''}" \
        --significance ${params.de_significance} \
        --cutoff ${params.de_cutoff} \
        --outdir . >> dge_summary.txt

    echo "" >> dge_summary.txt
    echo "Additional DGE output counts:" >> dge_summary.txt

    python3 - <<'PY' >> dge_summary.txt
import csv

path = "differential_expression_results.csv"

with open(path, newline="") as f:
    rows = list(csv.DictReader(f))

n_tested = len(rows)
n_sig = 0
n_up = 0
n_down = 0

significance = "${params.de_significance}"
cutoff = float("${params.de_cutoff}")

for row in rows:
    try:
        fc = float(row.get("logFC", "nan"))

        if significance == "pvalue":
            sig_value = float(row.get("PValue", "nan"))
        elif significance == "fdr":
            sig_value = float(row.get("FDR", "nan"))
        else:
            continue

    except (TypeError, ValueError):
        continue

    if sig_value < cutoff:
        n_sig += 1

        if fc > 0:
            n_up += 1
        elif fc < 0:
            n_down += 1

print(f"Mature miRNAs tested: {n_tested}")
print(f"Significance metric: {significance}")
print(f"Significance cutoff: {cutoff}")
print(f"Significant mature miRNAs: {n_sig}")
print(f"Significant up-regulated in comparison: {n_up}")
print(f"Significant down-regulated in comparison: {n_down}")
PY
    """
}

// =========================================================
// 3. SPLIT SIGNIFICANT miRNAs
// =========================================================
process SPLIT_MIRNAS_SMRNA {
    tag "split significant mature miRNAs"
    publishDir "${params.outdir}/targets", mode: 'copy'

    input:
    path dge_results
    val de_significance
    val de_cutoff

    output:
    path "Up_Regulated/input_mirnas.txt", emit: up_mirnas
    path "Down_Regulated/input_mirnas.txt", emit: down_mirnas
    path "Up_and_Down/input_mirnas.txt", emit: combined_mirnas
    path "split_summary.txt", emit: summary

    script:
    """
    echo "============================================================" > split_summary.txt
    echo "SPLIT_MIRNAS_SMRNA" >> split_summary.txt
    echo "============================================================" >> split_summary.txt
    echo "DGE results: ${dge_results}" >> split_summary.txt
    echo "Significance metric: ${de_significance}" >> split_summary.txt
    echo "Significance cutoff: ${de_cutoff}" >> split_summary.txt
    echo "" >> split_summary.txt

    Rscript ${projectDir}/scripts/split_mirnas_smrna.R \\
        --dge ${dge_results} \\
        --significance ${de_significance} \\
        --cutoff ${de_cutoff} \\
        --outdir . >> split_summary.txt
    """
}

// =========================================================
// 4. VALIDATED miRNA TARGET GENES
// =========================================================
process TARGET_GENE_SMRNA {
    tag "${direction} validated targets"
    publishDir "${params.outdir}/targets/${direction}", mode: 'copy'

    input:
    tuple val(direction), path(mirna_list)

    output:
    tuple val(direction), path("predicted_targets.csv"), emit: predicted_targets
    tuple val(direction), path("mirna_validated_targets_details.csv"), emit: target_details
    tuple val(direction), path("unique_target_genes_list.csv"), emit: unique_target_genes
    path "target_summary_${direction}.txt", emit: summary

    script:
    """
    echo "============================================================" > target_summary_${direction}.txt
    echo "TARGET_GENE_SMRNA - ${direction}" >> target_summary_${direction}.txt
    echo "============================================================" >> target_summary_${direction}.txt
    echo "Input miRNA list: ${mirna_list}" >> target_summary_${direction}.txt
    echo "Direction: ${direction}" >> target_summary_${direction}.txt
    echo "" >> target_summary_${direction}.txt

    Rscript ${projectDir}/scripts/target_gene_smrna.R \\
        --mirnas ${mirna_list} \\
        --outdir . >> target_summary_${direction}.txt
    """
}

// =========================================================
// 5. GO BP / CC / MF + KEGG ENRICHMENT
// =========================================================
process ENRICHMENT_SMRNA {
    tag "${direction} enrichment"
    publishDir "${params.outdir}/enrichment/${direction}", mode: 'copy'

    input:
    tuple val(direction), path(target_genes)
    val enrichment_significance
    val enrichment_cutoff

    output:
    path "GO_BP_enrichment_results.csv", emit: go_bp
    path "GO_CC_enrichment_results.csv", emit: go_cc
    path "GO_MF_enrichment_results.csv", emit: go_mf
    path "KEGG_enrichment_results.csv", emit: kegg

    path "GO_dotplot_BP.png", emit: go_bp_plot
    path "GO_dotplot_CC.png", emit: go_cc_plot
    path "GO_dotplot_MF.png", emit: go_mf_plot
    path "KEGG_dotplot.png", emit: kegg_plot

    path "enrichment_summary_${direction}.txt", emit: summary

    script:
    """
    echo "============================================================" > enrichment_summary_${direction}.txt
    echo "ENRICHMENT_SMRNA - ${direction}" >> enrichment_summary_${direction}.txt
    echo "============================================================" >> enrichment_summary_${direction}.txt
    echo "Target gene file: ${target_genes}" >> enrichment_summary_${direction}.txt
    echo "Direction: ${direction}" >> enrichment_summary_${direction}.txt
    echo "Significance metric: ${enrichment_significance}" >> enrichment_summary_${direction}.txt
    echo "Significance cutoff: ${enrichment_cutoff}" >> enrichment_summary_${direction}.txt
    echo "" >> enrichment_summary_${direction}.txt

    echo "--- GO Biological Process ---" >> enrichment_summary_${direction}.txt
    Rscript ${projectDir}/scripts/enrichment_smrna.R \
        --targets ${target_genes} \
        --outdir . \
        --direction ${direction} \
        --ont BP \
        --significance ${enrichment_significance} \
        --cutoff ${enrichment_cutoff} >> enrichment_summary_${direction}.txt

    echo "" >> enrichment_summary_${direction}.txt
    echo "--- GO Cellular Component ---" >> enrichment_summary_${direction}.txt
    Rscript ${projectDir}/scripts/enrichment_smrna.R \
        --targets ${target_genes} \
        --outdir . \
        --direction ${direction} \
        --ont CC \
        --significance ${enrichment_significance} \
        --cutoff ${enrichment_cutoff} >> enrichment_summary_${direction}.txt

    echo "" >> enrichment_summary_${direction}.txt
    echo "--- GO Molecular Function ---" >> enrichment_summary_${direction}.txt
    Rscript ${projectDir}/scripts/enrichment_smrna.R \
        --targets ${target_genes} \
        --outdir . \
        --direction ${direction} \
        --ont MF \
        --significance ${enrichment_significance} \
        --cutoff ${enrichment_cutoff} >> enrichment_summary_${direction}.txt

    echo "" >> enrichment_summary_${direction}.txt
    echo "--- KEGG ---" >> enrichment_summary_${direction}.txt
    Rscript ${projectDir}/scripts/enrichment_smrna.R \
        --targets ${target_genes} \
        --outdir . \
        --direction ${direction} \
        --kegg \
        --significance ${enrichment_significance} \
        --cutoff ${enrichment_cutoff} >> enrichment_summary_${direction}.txt

    echo "" >> enrichment_summary_${direction}.txt
    echo "Enrichment result counts:" >> enrichment_summary_${direction}.txt

    csv_count() {
        n=\$(wc -l < "\$1" | tr -d ' ')
        if [ "\$n" -gt 1 ]; then
            echo \$((n - 1))
        else
            echo 0
        fi
    }

    echo "GO BP significant terms: \$(csv_count GO_BP_enrichment_results.csv)" >> enrichment_summary_${direction}.txt
    echo "GO CC significant terms: \$(csv_count GO_CC_enrichment_results.csv)" >> enrichment_summary_${direction}.txt
    echo "GO MF significant terms: \$(csv_count GO_MF_enrichment_results.csv)" >> enrichment_summary_${direction}.txt
    echo "KEGG significant pathways: \$(csv_count KEGG_enrichment_results.csv)" >> enrichment_summary_${direction}.txt
    """
}

// =========================================================
// 6. VISUALIZATION
// =========================================================
process VISUALIZATION_SMRNA {
    tag "mature miRNA visualizations"
    publishDir "${params.outdir}/plots", mode: 'copy'

    input:
    path dge_rds
    path dge_results
    path metadata
    val controls

    output:
    path "PCA_normalized_miRNA_expression.png", emit: pca
    path "Significant_miRNAs_Heatmap.png", emit: heatmap
    path "Volcano_mature_miRNAs.png", emit: volcano
    path "visualization_summary.txt", emit: summary

    script:
    """
    echo "============================================================" > visualization_summary.txt
    echo "VISUALIZATION_SMRNA" >> visualization_summary.txt
    echo "============================================================" >> visualization_summary.txt
    echo "DGE RDS: ${dge_rds}" >> visualization_summary.txt
    echo "DGE results: ${dge_results}" >> visualization_summary.txt
    echo "Metadata: ${metadata}" >> visualization_summary.txt
    echo "Controls/reference samples: ${controls ?: '<not specified>'}" >> visualization_summary.txt
    echo "" >> visualization_summary.txt

    echo "--- PCA ---" >> visualization_summary.txt
    Rscript ${projectDir}/scripts/PCA_plot_smrna.R \\
        --dge-rds ${dge_rds} \\
        --metadata ${metadata} \\
        --outdir . >> visualization_summary.txt

    echo "" >> visualization_summary.txt
    echo "--- Heatmap ---" >> visualization_summary.txt
    Rscript ${projectDir}/scripts/heatmap_smrna.R \
        --dge-rds dge.rds \
        --dge differential_expression_results.csv \
        --metadata metadata.csv \
        --outdir . \
        --significance ${params.de_significance} \
        --cutoff ${params.de_cutoff} >> visualization_summary.txt

    echo "" >> visualization_summary.txt
    echo "--- Volcano ---" >> visualization_summary.txt
    Rscript ${projectDir}/scripts/volcano_plot_smrna.R \
        --dge-rds ${dge_rds} \
        --dge ${dge_results} \
        --metadata ${metadata} \
        --controls "${controls ?: ''}" \
        --significance ${params.de_significance} \
        --cutoff ${params.de_cutoff} \
        --outdir . >> visualization_summary.txt

    echo "" >> visualization_summary.txt
    echo "Plot files:" >> visualization_summary.txt

    for f in PCA_normalized_miRNA_expression.png Significant_miRNAs_Heatmap.png Volcano_mature_miRNAs.png; do
        if [ -s "\$f" ]; then
            echo "Generated: \$f" >> visualization_summary.txt
        else
            echo "Missing: \$f" >> visualization_summary.txt
        fi
    done
    """
}

// =========================================================
// 7. COMPILE ONE FINAL PIPELINE SUMMARY
// =========================================================
process SUMMARY_COMPILE {
    tag "compile pipeline summary"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path summaries

    output:
    path "pipeline_summary.txt"

    script:
    """
    echo "============================================================" > pipeline_summary.txt
    echo "SMALL RNA-SEQ DOWNSTREAM ANALYSIS SUMMARY" >> pipeline_summary.txt
    echo "============================================================" >> pipeline_summary.txt
    echo "" >> pipeline_summary.txt
    echo "Output directory: ${params.outdir}" >> pipeline_summary.txt
    echo "Controls/reference samples: ${params.controls ?: '<not specified>'}" >> pipeline_summary.txt
    echo "" >> pipeline_summary.txt

    echo "============================================================" >> pipeline_summary.txt
    echo "1. FILTERING" >> pipeline_summary.txt
    echo "============================================================" >> pipeline_summary.txt
    cat filter_summary.txt >> pipeline_summary.txt
    echo "" >> pipeline_summary.txt

    echo "============================================================" >> pipeline_summary.txt
    echo "2. DIFFERENTIAL EXPRESSION" >> pipeline_summary.txt
    echo "============================================================" >> pipeline_summary.txt
    cat dge_summary.txt >> pipeline_summary.txt
    echo "" >> pipeline_summary.txt

    echo "============================================================" >> pipeline_summary.txt
    echo "3. SIGNIFICANT miRNA SPLITTING" >> pipeline_summary.txt
    echo "============================================================" >> pipeline_summary.txt
    cat split_summary.txt >> pipeline_summary.txt
    echo "" >> pipeline_summary.txt

    echo "============================================================" >> pipeline_summary.txt
    echo "4. VALIDATED miRNA TARGETS" >> pipeline_summary.txt
    echo "============================================================" >> pipeline_summary.txt
    for f in target_summary_*.txt; do
        [ -e "\$f" ] || continue
        cat "\$f" >> pipeline_summary.txt
        echo "" >> pipeline_summary.txt
    done

    echo "============================================================" >> pipeline_summary.txt
    echo "5. FUNCTIONAL ENRICHMENT" >> pipeline_summary.txt
    echo "============================================================" >> pipeline_summary.txt
    for f in enrichment_summary_*.txt; do
        [ -e "\$f" ] || continue
        cat "\$f" >> pipeline_summary.txt
        echo "" >> pipeline_summary.txt
    done

    echo "============================================================" >> pipeline_summary.txt
    echo "6. VISUALIZATION" >> pipeline_summary.txt
    echo "============================================================" >> pipeline_summary.txt
    cat visualization_summary.txt >> pipeline_summary.txt
    echo "" >> pipeline_summary.txt

    echo "============================================================" >> pipeline_summary.txt
    echo "PIPELINE COMPLETED" >> pipeline_summary.txt
    echo "============================================================" >> pipeline_summary.txt
    """
}

// =========================================================
// 8. CONVERT ALL CSV RESULTS TO XLSX
// =========================================================
process CSV_TO_XLSX {
    tag "convert CSV results to XLSX"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val trigger

    output:
    path "csv_to_xlsx.done"

    script:
    """
    python3 ${projectDir}/scripts/csv_to_xlsx.py \\
        --results ${params.outdir}

    touch csv_to_xlsx.done
    """
}

// =========================================================
// WORKFLOW
// =========================================================
workflow {

    metadata_ch = Channel.value(
        file(params.metadata, checkIfExists: true)
    )

    log.info "=============================================="
    log.info "SMALL RNA-SEQ CUSTOM METADATA:"
    log.info "${params.metadata}"
    log.info "CONTROLS:"
    log.info "${params.controls ?: '<not specified>'}"
    log.info "=============================================="

    smrna_counts_ch = Channel.fromPath(
        params.smrna_counts,
        checkIfExists: true
    )

    smrna_logtpm_ch = Channel.fromPath(
        params.smrna_logtpm,
        checkIfExists: true
    )

    // -----------------------------------------------------
    // FILTER
    // -----------------------------------------------------
    filter_out = FILTER_SMRNA(
        smrna_counts_ch,
        smrna_logtpm_ch,
        metadata_ch,
        params.controls ?: ''
    )

    // -----------------------------------------------------
    // DIFFERENTIAL EXPRESSION
    // -----------------------------------------------------
    dge_out = DIFF_EXP_SMRNA(
        filter_out.counts_filtered,
        metadata_ch,
        params.controls ?: ''
    )

    // -----------------------------------------------------
    // SPLIT SIGNIFICANT miRNAs
    // -----------------------------------------------------
    split_out = SPLIT_MIRNAS_SMRNA(
        dge_out.dge_results,
        params.de_significance,
        params.de_cutoff
    )

    // -----------------------------------------------------
    // THREE TARGET ANALYSIS BRANCHES
    //   1. Up_Regulated
    //   2. Down_Regulated
    //   3. Up_and_Down (all significant)
    // -----------------------------------------------------
    up_input = split_out.up_mirnas.map { f ->
        tuple('Up_Regulated', f)
    }

    down_input = split_out.down_mirnas.map { f ->
        tuple('Down_Regulated', f)
    }

    combined_input = split_out.combined_mirnas.map { f ->
        tuple('Up_and_Down', f)
    }

    target_input = up_input
        .mix(down_input)
        .mix(combined_input)

    target_outputs = TARGET_GENE_SMRNA(target_input)

    // -----------------------------------------------------
    // ENRICHMENT FOR ALL THREE TARGET BRANCHES
    // -----------------------------------------------------
    // TARGET_GENE_SMRNA already emits tuples:
    // (direction, unique_target_genes.csv)
    enrichment_outputs = ENRICHMENT_SMRNA(
       target_outputs.unique_target_genes,
       params.enrichment_significance,
       params.enrichment_cutoff
    )

    // -----------------------------------------------------
    // VISUALIZATION
    // -----------------------------------------------------
    viz_out = VISUALIZATION_SMRNA(
        dge_out.dge_rds,
        dge_out.dge_results,
        metadata_ch,
        params.controls ?: ''
    )

    // -----------------------------------------------------
    // COMBINE STAGE SUMMARIES
    // -----------------------------------------------------
    all_summaries = filter_out.summary
        .mix(dge_out.summary)
        .mix(split_out.summary)
        .mix(target_outputs.summary)
        .mix(enrichment_outputs.summary)
        .mix(viz_out.summary)
        .collect()

    summary_out = SUMMARY_COMPILE(all_summaries)

    // -----------------------------------------------------
    // XLSX CONVERSION IS THE FINAL STEP
    // -----------------------------------------------------
    CSV_TO_XLSX(summary_out)
}
