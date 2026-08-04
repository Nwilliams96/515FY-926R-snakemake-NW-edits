# import basic packages
import os
import pandas as pd
from snakemake.utils import validate

# Resolve optional workflow branches and the shared database location once so
# every rule uses the same paths. Database files always live outside the
# run-specific repository clone.
USE_PREEXISTING_DATABASES = config.get("use_preexisting_databases", False)
USE_INTERNAL_STANDARDS = config.get("use_internal_standards", False)

DATABASE_DIR = os.path.normpath(config["database_dir"])
DATABASE_PREFIX = DATABASE_DIR + os.sep
BBSPLIT_DB_DIR = os.path.join(
    DATABASE_DIR, "bbsplit-db", "EUK-PROK-bbsplit-db"
)
SILVA_CLASSIFIER = os.path.join(
    DATABASE_DIR,
    "classification",
    "SILVA",
    "silva-ssu-nr99-tax-dereplicated-sliced_"
    + config["fwdPrimer"]
    + "_"
    + config["revPrimer"]
    + "_dereplicated_final_classifier_USE_ME.qza",
)
PR2_CLASSIFIER = os.path.join(
    DATABASE_DIR,
    "classification",
    "PR2",
    "pr2_version_5.1.1_SSU_dada2.clean.culled.derep-sliced_"
    + config["fwdPrimer"]
    + "_"
    + config["revPrimer"]
    + "_dereplicated_final_classifier_USE_ME.qza",
)

# Small, user-facing files copied into one folder at the end of the run.
RESULTS_EXPORT_INPUTS = [
    "results/04-formatted/" + config["studyName"] + ".long_data.tsv",
    "results/04-formatted/" + config["studyName"] + ".asv_sequences.tsv",
    "results/07-report/" + config["studyName"] + ".pipeline-report.html",
]
RESULTS_EXPORT_OUTPUTS = [
    "Results-Export/" + config["studyName"] + ".long_data.tsv",
    "Results-Export/" + config["studyName"] + ".asv_sequences.tsv",
    "Results-Export/" + config["studyName"] + ".pipeline-report.html",
]
if USE_INTERNAL_STANDARDS:
    RESULTS_EXPORT_INPUTS.append(
        "results/05-internal-std-corrected/"
        + config["studyName"]
        + ".ISD_corrected_asv_table.tsv"
    )
    RESULTS_EXPORT_OUTPUTS.append(
        "Results-Export/"
        + config["studyName"]
        + ".ISD_corrected_asv_table.tsv"
    )

# read sample sheet
samples = (
    pd.read_csv(config["samplesheet"], sep="\t", dtype={"sample": str})
    .set_index("sample", drop=False)
    .sort_index()
)

# define output as function
def get_final_output():
    final_output = expand(
        "results/01-split/{sample}.{organism}.R1.fastq.gz", sample=samples["sample"], organism=["prok","euk"]
        #"results/01-split/{sample}.prok.fastq",
        #"results/00-trimmed/{sample}.{direction}.fastq",
        #sample=samples["sample"], direction=["1","2"]
    )

    final_output.extend(
        [directory(BBSPLIT_DB_DIR), SILVA_CLASSIFIER, PR2_CLASSIFIER]
    )
#    final_output.append("results/02-proks/manifest.tsv"),
    final_output.append("results/02-proks/16S.qza"),
    final_output.append(directory("results/02-proks/02-quality-plots-R1-R2/")),
    final_output.append(directory("results/02-proks/03-DADA2d/")),
    final_output.append(directory("results/02-proks/04-DADA2d-plaintext-exports")),
    final_output.append(directory("results/02-proks/05-classified")),
    final_output.append(directory("results/02-proks/07-SILVA-only-barplots/")),
    final_output.append("results/02-proks/09-subsetting/split-seqs/exclude_o__Chloroplast_subset_filtered_seqs.qza"),
    final_output.append("results/02-proks/09-subsetting/reclassified/include_o__Chloroplast_subset_reclassified_PR2.qza"),
    final_output.append("results/02-proks/10-exports/" + config["studyName"] + ".taxonomy.tsv"),
    final_output.append("results/02-proks/10-exports/" + config["studyName"] + ".all-16S-seqs.with-tax.tsv"),
    final_output.append("results/02-proks/sample-metadata.tsv"),
    final_output.append("results/02-euks/18S-viz.qza"),
    final_output.append(directory("results/02-euks/02-quality-plots-R1-R2/")),
    final_output.append(directory("results/02-euks/07-quality-plots-concat")),
    final_output.append(directory("results/02-euks/08-DADA2d")),
    final_output.append("results/02-euks/18S-concat.qza"),
    final_output.append(expand("results/02-euks/04-concatenated/{sample}.euk.concatenated.fastq", sample=samples["sample"])),
    final_output.append(directory("results/02-euks/09-DADA2d-plaintext-exports/")),
    final_output.append(directory("results/02-euks/10-classified/")),
    final_output.append("results/02-euks/sample-metadata.tsv"),
    final_output.append(directory("results/02-euks/12-SILVA-only-barplots/")),
    final_output.append("results/02-euks/14-subsetting/reclassified-PR2/fixed/taxonomy-without-spaces.qza"),
    final_output.append("results/02-euks/14-subsetting/split-tables/include_Metazoa_PR2_filtered_table.qza"),
    final_output.append("results/02-euks/15-exports/" + config["studyName"] + ".include_Metazoa_PR2_filtered_table.with-tax.tsv"),
    final_output.append("results/" + config["studyName"] + ".eukfrac-per-sample.tsv"),
    final_output.append("results/" + config["studyName"] + ".eukfrac-all.tsv"),
    final_output.append("results/03-merged/" + config["studyName"] + ".merged_uncorrected.tsv"),
    final_output.append("results/02-proks/10-exports/" + config["studyName"] + ".Synechococcales.proportal-classified.tsv"),
    final_output.append("results/04-formatted/" + config["studyName"] + ".long_data.tsv"),
    final_output.append(
        "results/07-report/" + config["studyName"] + ".pipeline-report.html"
    )
    if USE_INTERNAL_STANDARDS:
        final_output.append(
            "results/05-internal-std-corrected/"
            + config["studyName"]
            + ".ISD_corrected_asv_table.tsv"
        )

    final_output.extend(RESULTS_EXPORT_OUTPUTS)
    final_output.append("Results-Export/README.txt")

    return final_output

# validate sample sheet and config file
validate(samples, schema="../../config/schemas/samples.schema.yml")
validate(config, schema="../../config/schemas/config.schema.yml")
