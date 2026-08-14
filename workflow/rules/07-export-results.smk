rule export_results_for_download:
    input:
        data_files=RESULTS_EXPORT_INPUTS
    output:
        data_files=RESULTS_EXPORT_OUTPUTS,
        summaries=RESULTS_EXPORT_SUMMARIES,
        readme=RESULTS_EXPORT_DIR + "/README.txt"
    script:
        "../scripts/export_results.py"
