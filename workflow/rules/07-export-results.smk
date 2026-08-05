rule export_results_for_download:
    input:
        data_files=RESULTS_EXPORT_INPUTS
    output:
        data_files=RESULTS_EXPORT_OUTPUTS,
        readme="Results-Export/README.txt"
    script:
        "../scripts/export_results.py"
