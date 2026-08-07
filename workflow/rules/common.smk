# import basic packages
import os
import re
from itertools import combinations
import pandas as pd
from snakemake.exceptions import WorkflowError
from snakemake.utils import validate

# Resolve optional workflow branches and the shared database location once so
# every rule uses the same paths. Database files always live outside the
# run-specific repository clone.
USE_PREEXISTING_DATABASES = config.get("use_preexisting_databases", False)
USE_INTERNAL_STANDARDS = config.get("use_internal_standards", False)

# New configs store internal standards as an ordered YAML list. Continue to
# accept the older intstd1/intstd2/intstd3 mapping so existing studies remain
# runnable after updating the workflow.
configured_intstds = config.get("intstds", [])
if isinstance(configured_intstds, dict):
    INTERNAL_STANDARD_IDS = [
        str(value).strip() for value in configured_intstds.values()
    ]
elif isinstance(configured_intstds, list):
    INTERNAL_STANDARD_IDS = [str(value).strip() for value in configured_intstds]
else:
    raise WorkflowError("intstds must be an ordered list of internal-standard IDs")

if USE_INTERNAL_STANDARDS:
    if not INTERNAL_STANDARD_IDS:
        raise WorkflowError(
            "use_internal_standards is true, but intstds does not contain any IDs"
        )
    invalid_internal_standard_ids = [
        standard_id
        for standard_id in INTERNAL_STANDARD_IDS
        if not re.fullmatch(r"[A-Za-z0-9._-]+", standard_id)
    ]
    if invalid_internal_standard_ids:
        raise WorkflowError(
            "Internal-standard IDs may contain only letters, numbers, periods, "
            "underscores, and hyphens: "
            + ", ".join(invalid_internal_standard_ids)
        )
    if len(set(INTERNAL_STANDARD_IDS)) != len(INTERNAL_STANDARD_IDS):
        raise WorkflowError("Internal-standard IDs in intstds must be unique")

INTERNAL_STANDARD_SLOTS = [
    f"isd_{index}" for index in range(1, len(INTERNAL_STANDARD_IDS) + 1)
]
INTERNAL_STANDARD_PAIRS = list(combinations(INTERNAL_STANDARD_SLOTS, 2))
INTERNAL_STANDARD_METHOD_STEMS = (
    [f"{slot}_recovery_ratio" for slot in INTERNAL_STANDARD_SLOTS]
    + ["mean_recovery_ratio", "median_recovery_ratio"]
    + [
        f"{first}_{second}_mean_recovery_ratio"
        for first, second in INTERNAL_STANDARD_PAIRS
    ]
)

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

if USE_PREEXISTING_DATABASES:
    missing_database_resources = [
        path
        for path in [BBSPLIT_DB_DIR, SILVA_CLASSIFIER, PR2_CLASSIFIER]
        if not os.path.exists(path)
    ]
    if missing_database_resources:
        missing_list = "\n  - ".join(missing_database_resources)
        raise WorkflowError(
            "use_preexisting_databases is true, but the following configured "
            "database resources were not found:\n  - "
            + missing_list
            + "\nSet database_dir to the shared database root that contains "
            "bbsplit-db/ and classification/, or set "
            "use_preexisting_databases to false so the workflow builds them."
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

if USE_INTERNAL_STANDARDS:
    expected_internal_standard_columns = [
        f"{standard_id}_ng" for standard_id in INTERNAL_STANDARD_IDS
    ]
    missing_internal_standard_columns = [
        column
        for column in expected_internal_standard_columns
        if column not in samples.columns
    ]
    if missing_internal_standard_columns:
        raise WorkflowError(
            "The sample sheet is missing internal-standard amount column(s): "
            + ", ".join(missing_internal_standard_columns)
        )

# define output as function
def get_final_output():
    final_output = expand(
        "results/01-split/{sample}.{organism}.R1.fastq.gz", sample=samples["sample"], organism=["prok","euk"]
        #"results/01-split/{sample}.prok.fastq",
        #"results/00-trimmed/{sample}.{direction}.fastq",
        #sample=samples["sample"], direction=["1","2"]
    )

    final_output.extend([BBSPLIT_DB_DIR, SILVA_CLASSIFIER, PR2_CLASSIFIER])
#    final_output.append("results/02-proks/manifest.tsv"),
    final_output.append("results/02-proks/16S.qza"),
    final_output.append("results/02-proks/02-quality-plots-R1-R2/"),
    final_output.append("results/02-proks/03-DADA2d/"),
    final_output.append("results/02-proks/04-DADA2d-plaintext-exports"),
    final_output.append("results/02-proks/05-classified"),
    final_output.append("results/02-proks/07-SILVA-only-barplots/"),
    final_output.append("results/02-proks/09-subsetting/split-seqs/exclude_o__Chloroplast_subset_filtered_seqs.qza"),
    final_output.append("results/02-proks/09-subsetting/reclassified/include_o__Chloroplast_subset_reclassified_PR2.qza"),
    final_output.append("results/02-proks/10-exports/" + config["studyName"] + ".taxonomy.tsv"),
    final_output.append("results/02-proks/10-exports/" + config["studyName"] + ".all-16S-seqs.with-tax.tsv"),
    final_output.append("results/02-proks/sample-metadata.tsv"),
    final_output.append("results/02-euks/18S-viz.qza"),
    final_output.append("results/02-euks/02-quality-plots-R1-R2/"),
    final_output.append("results/02-euks/07-quality-plots-concat"),
    final_output.append("results/02-euks/08-DADA2d"),
    final_output.append("results/02-euks/18S-concat.qza"),
    final_output.append(expand("results/02-euks/04-concatenated/{sample}.euk.concatenated.fastq", sample=samples["sample"])),
    final_output.append("results/02-euks/09-DADA2d-plaintext-exports/"),
    final_output.append("results/02-euks/10-classified/"),
    final_output.append("results/02-euks/sample-metadata.tsv"),
    final_output.append("results/02-euks/12-SILVA-only-barplots/"),
    final_output.append("results/02-euks/14-subsetting/reclassified-PR2/fixed/taxonomy-without-spaces.qza"),
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
