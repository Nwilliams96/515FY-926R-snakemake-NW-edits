This workflow is in progress, and not yet ready for production...

Please continue to use the bash workflow until otherwise indicated on the base README.

Instructions for Nathan:

1. Install latest version of snakemake via conda, making sure to do the full snakemake install instead of "minimal"
[Snakemake install instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
2. Activate snakemake conda env
```
conda activate snakemake
```
3. Install slurm plugin
```pip install snakemake-executor-plugin-slurm```
4. Clone the latest version of the eASV pipeline:
```
git clone https://github.com/jcmcnch/eASV-pipeline-for-515Y-926R.git
```
5. Enter repo, navigate to `snakemake` folder
```
cd eASV-pipeline-for-515Y-926R/snakemake
```
6. Clone another repo with test data (yes, inside the other repo)
```
git clone https://github.com/stfx-microeco-lab/nativemicrobiota
```
7. Now you should have all the data needed to run a toy version of the pipeline. Here's how I was able to successfully run it on my laptop in case you want to test it on yours first before running on CARC:
```
snakemake --use-conda --cores 4
```
8. To run with slurm, it should be as simple as running the following command:
```
snakemake --executor slurm --default-resources --jobs N #replacing "N" with the number of jobs
#or possibly with the --use-conda flag
snakemake --use-conda --executor slurm --default-resources --jobs N
```

## Run report

The final workflow target creates a self-contained report at:

```text
results/07-report/<studyName>.pipeline-report.html
```

Open that file in a web browser for a concise summary of the configuration,
read retention through trimming, 16S/18S splitting and DADA2, sample-level QC,
domain composition, abundant taxa, and—when enabled—the internal-standard
figures. The report is intended for rapid review; retain the complete results
directory and QIIME 2 visualizations for detailed analysis.

For convenient downloading, the workflow also copies the report and formatted
data tables into the top-level `Results-Export/` folder. When internal standards
are enabled, the corrected ASV table is included there as well.

For your reference, in case the above commands are incomplete:
https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html#usage
