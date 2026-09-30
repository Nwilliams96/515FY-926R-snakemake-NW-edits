INTERNAL_STANDARD_REPORT_FIGURES = []
INTERNAL_STANDARD_REPORT_TABLE = []
UNMERGED_16S_REPORT_SUMMARY = []
if USE_INTERNAL_STANDARDS:
    INTERNAL_STANDARD_REPORT_FIGURES = [
        rules.intstd_correct_data.output.recovery_plot_png,
        rules.intstd_correct_data.output.domain_plot_png,
    ]
    INTERNAL_STANDARD_REPORT_TABLE = [rules.intstd_correct_data.output.corrected]
if RUN_UNMERGED_16S_RESCUE:
    UNMERGED_16S_REPORT_SUMMARY = [UNMERGED_16S_SUMMARY]


rule generate_pipeline_report:
    input:
        samples="config/samples.tsv",
        split_summary="results/" + config["studyName"] + ".eukfrac-per-sample.tsv",
        stats16s="results/02-proks/04-DADA2d-plaintext-exports/" + config["studyName"] + ".16S.latest_stats.tsv",
        stats18s="results/02-euks/09-DADA2d-plaintext-exports/" + config["studyName"] + ".18S.latest_stats.tsv",
        correction_factors=rules.merge_prok_euk.output.correction_factors,
        species_assignment=rules.assign_SILVA_144_species.output.summary,
        chloroplast_audit=rules.resolve_PR2_plastids.output.audit,
        chloroplast_summary=rules.resolve_PR2_plastids.output.summary,
        cutadapt_qc=expand(
            "results/00-trimmed/{sample}.qc.txt", sample=samples["sample"]
        ),
        quality_16s="results/02-proks/02-quality-plots-R1-R2/",
        quality_18s_paired="results/02-euks/02-quality-plots-R1-R2/",
        quality_18s_concatenated="results/02-euks/07-quality-plots-concat/",
        long_data=RESULTS_LONG_DATA,
        internal_standard_figures=INTERNAL_STANDARD_REPORT_FIGURES,
        internal_standard_table=INTERNAL_STANDARD_REPORT_TABLE,
        unmerged_16s_summary=UNMERGED_16S_REPORT_SUMMARY
    output:
        html="results/07-report/" + config["studyName"] + ".pipeline-report.html"
    log:
        "logs/06-report/generate_pipeline_report.log"
    script:
        "../scripts/generate_pipeline_report.py"
