INTERNAL_STANDARD_REPORT_FIGURES = []
if USE_INTERNAL_STANDARDS:
    INTERNAL_STANDARD_REPORT_FIGURES = [
        rules.intstd_correct_data.output.recovery_plot_png,
        rules.intstd_correct_data.output.domain_plot_png,
    ]


rule generate_pipeline_report:
    input:
        samples="config/samples.tsv",
        split_summary="results/" + config["studyName"] + ".eukfrac-per-sample.tsv",
        stats16s="results/02-proks/04-DADA2d-plaintext-exports/" + config["studyName"] + ".16S.latest_stats.tsv",
        stats18s="results/02-euks/09-DADA2d-plaintext-exports/" + config["studyName"] + ".18S.latest_stats.tsv",
        long_data="results/04-formatted/" + config["studyName"] + ".long_data.tsv",
        quality_directories=[
            "results/02-proks/02-quality-plots-R1-R2/",
            "results/02-euks/02-quality-plots-R1-R2/",
            "results/02-euks/07-quality-plots-concat/",
        ],
        trimming_logs=expand("logs/00-trimming/{sample}.log", sample=samples["sample"]),
        internal_standard_figures=INTERNAL_STANDARD_REPORT_FIGURES
    output:
        html="results/07-report/" + config["studyName"] + ".pipeline-report.html"
    log:
        "logs/06-report/generate_pipeline_report.log"
    script:
        "../scripts/generate_pipeline_report.py"
