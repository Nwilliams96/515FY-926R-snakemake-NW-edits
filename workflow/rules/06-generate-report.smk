INTERNAL_STANDARD_REPORT_FIGURES = []
INTERNAL_STANDARD_REPORT_TABLE = []
if USE_INTERNAL_STANDARDS:
    INTERNAL_STANDARD_REPORT_FIGURES = [
        rules.intstd_correct_data.output.recovery_plot_png,
        rules.intstd_correct_data.output.domain_plot_png,
    ]
    INTERNAL_STANDARD_REPORT_TABLE = [rules.intstd_correct_data.output.corrected]


rule generate_pipeline_report:
    input:
        samples="config/samples.tsv",
        split_summary="results/" + config["studyName"] + ".eukfrac-per-sample.tsv",
        stats16s="results/02-proks/04-DADA2d-plaintext-exports/" + config["studyName"] + ".16S.latest_stats.tsv",
        stats18s="results/02-euks/09-DADA2d-plaintext-exports/" + config["studyName"] + ".18S.latest_stats.tsv",
        cutadapt_qc=expand(
            "results/00-trimmed/{sample}.qc.txt", sample=samples["sample"]
        ),
        long_data="results/04-formatted/" + config["studyName"] + ".long_data.tsv",
        internal_standard_figures=INTERNAL_STANDARD_REPORT_FIGURES,
        internal_standard_table=INTERNAL_STANDARD_REPORT_TABLE
    output:
        html="results/07-report/" + config["studyName"] + ".pipeline-report.html"
    log:
        "logs/06-report/generate_pipeline_report.log"
    script:
        "../scripts/generate_pipeline_report.py"
