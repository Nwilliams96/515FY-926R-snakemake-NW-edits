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
# The three database families can be managed independently. The legacy setting
# remains the fallback for BBsplit because that was the only database branch it
# controlled. SILVA and PR2 retain their historical create-if-missing behavior
# when the new settings are absent.
USE_PREEXISTING_BBSPLIT_DATABASE = config.get(
    "use_preexisting_bbsplit_database",
    config.get("use_preexisting_databases", False),
)
USE_PREEXISTING_SILVA_DATABASE = config.get(
    "use_preexisting_silva_database", False
)
USE_PREEXISTING_PR2_DATABASE = config.get(
    "use_preexisting_pr2_database", False
)
USE_INTERNAL_STANDARDS = config.get("use_internal_standards", False)

# SILVA 144 classifiers were built with QIIME 2 2026.7 and scikit-learn 1.7.1.
# Transparently migrate configs that still point at one of this repository's
# older bundled QIIME environments. Explicit external environment paths remain
# under the user's control.
DEFAULT_QIIME2_ENV = "../envs/rachis-qiime2-linux-64-2026.7.yml"
configured_qiime2_env = str(config.get("qiime2version", "")).strip()
if (
    not configured_qiime2_env
    or os.path.basename(configured_qiime2_env) in {
        "qiime2-amplicon-ubuntu-2024-5-conda.yml",
        "qiime2-amplicon-ubuntu-2025-7-conda.yml",
    }
):
    config["qiime2version"] = DEFAULT_QIIME2_ENV

SILVA_VERSION = "144"
config["SILVAversion"] = SILVA_VERSION

# PR2 is now used as an independent plastid screen for every 16S ASV.  A PR2
# lineage replaces SILVA only when it contains the explicit PR2 plastid marker
# ":plas" and its classifier confidence meets this threshold.  Keeping the
# resolver threshold separate from QIIME's classifier threshold gives the
# audit table enough information to record low-confidence plastid candidates.
CHLOROPLAST_PR2_MIN_CONFIDENCE = float(
    config.get("chloroplast_pr2_min_confidence", 0.7)
)
if not 0 <= CHLOROPLAST_PR2_MIN_CONFIDENCE <= 1:
    raise WorkflowError("chloroplast_pr2_min_confidence must be between 0 and 1")

# DADA2 parameters are configurable per marker-gene path. These fallbacks are
# the values used by the workflow before the controls were exposed, so older
# study configs remain reproducible after updating the pipeline.
DADA2_CONFIG = config.get("dada2", {})
DADA2_PROK_CONFIG = DADA2_CONFIG.get("prokaryotes", {})
DADA2_EUK_CONFIG = DADA2_CONFIG.get("eukaryotes", {})

DADA2_PROK_DEFAULTS = {
    "trunc_len_f": 0,
    "trunc_len_r": 0,
    "max_ee_f": 2.0,
    "max_ee_r": 2.0,
    "trunc_q": 2,
    "min_overlap": 12,
    "pooling_method": "independent",
    "chimera_method": "consensus",
    "min_fold_parent_over_abundance": 1.0,
    "n_reads_learn": 1000000,
}
DADA2_EUK_DEFAULTS = {
    "max_ee": 2.0,
    "trunc_q": 0,
    "pooling_method": "independent",
    "chimera_method": "consensus",
    "min_fold_parent_over_abundance": 1.0,
    "n_reads_learn": 1000000,
}

DADA2_PROK = {**DADA2_PROK_DEFAULTS, **DADA2_PROK_CONFIG}
DADA2_EUK = {**DADA2_EUK_DEFAULTS, **DADA2_EUK_CONFIG}

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
INTERNAL_STANDARD_ID_PAIRS = list(combinations(INTERNAL_STANDARD_IDS, 2))
ALL_INTERNAL_STANDARD_IDS_STEM = "_".join(INTERNAL_STANDARD_IDS)
ISD_CORRECTED_LONG_TABLE = (
    "results/05-internal-std-corrected/"
    + config["studyName"]
    + "."
    + ALL_INTERNAL_STANDARD_IDS_STEM
    + ".ISD_corrected_long_data.tsv"
)
ISD_FILTERED_LONG_TABLE = (
    "results/05-internal-std-corrected/"
    + config["studyName"]
    + ".ISD_removed_long_data.tsv"
)
INTERNAL_STANDARD_METHOD_STEMS = (
    [f"{standard_id}_recovery_ratio" for standard_id in INTERNAL_STANDARD_IDS]
    + [
        f"mean_{ALL_INTERNAL_STANDARD_IDS_STEM}_recovery_ratio",
        f"median_{ALL_INTERNAL_STANDARD_IDS_STEM}_recovery_ratio",
    ]
    + [
        f"mean_{first}_and_{second}_recovery_ratio"
        for first, second in INTERNAL_STANDARD_ID_PAIRS
    ]
)

NEW_AMPLICON_CONCENTRATIONS_FILE = (
    "config/prok_and_euk_SSU_amplicon_molarities.tsv"
)
PREVIOUS_AMPLICON_CONCENTRATIONS_FILE = (
    "config/prok_and_euk_SSU_amplicon_concentrations.tsv"
)
LEGACY_AMPLICON_CONCENTRATIONS_FILE = "config/bioanalyzer.tsv"
AMPLICON_CONCENTRATIONS_FILE = next(
    (
        path
        for path in [
            NEW_AMPLICON_CONCENTRATIONS_FILE,
            PREVIOUS_AMPLICON_CONCENTRATIONS_FILE,
            LEGACY_AMPLICON_CONCENTRATIONS_FILE,
        ]
        if os.path.exists(path)
    ),
    NEW_AMPLICON_CONCENTRATIONS_FILE,
)

DATABASE_DIR = os.path.normpath(config["database_dir"])
DATABASE_PREFIX = DATABASE_DIR + os.sep
BBSPLIT_DB_DIR = os.path.join(
    DATABASE_DIR, "bbsplit-db", "EUK-PROK-bbsplit-db"
)
SILVA_144_BASE_URL = (
    "https://www.arb-silva.de/fileadmin/silva_databases/current/"
    "QIIME2/2026.7/SSU/"
)
if (
    config["fwdPrimer"] == "GTGYCAGCMGCCGCGGTAA"
    and config["revPrimer"] == "CCGYCAATTYMTTTRAGTTT"
):
    SILVA_CLASSIFIER_FILENAME = (
        "SILVA_144_SSURef_NR99_uniform_classifier_V4V5-515f-926r.qza"
    )
    SILVA_CLASSIFIER_URL = (
        SILVA_144_BASE_URL
        + "V4V5-515f-926r/uniform/"
        + SILVA_CLASSIFIER_FILENAME
    )
    SILVA_CLASSIFIER_MD5 = "f7757b01eb82e0ac78bf06427e410095"
else:
    # A region-specific classifier would be invalid for other primer pairs.
    # Use SILVA's full-length classifier as the safe general SSU fallback.
    SILVA_CLASSIFIER_FILENAME = (
        "SILVA_144_SSURef_NR99_full-length_prokaryotes-"
        "rank-propagation_uniq_genus_uniform_classifier.qza"
    )
    SILVA_CLASSIFIER_URL = (
        SILVA_144_BASE_URL
        + "full-length/uniform/"
        + SILVA_CLASSIFIER_FILENAME
    )
    SILVA_CLASSIFIER_MD5 = "e16b9bc78dafa2f344563e570bbf25ba"

SILVA_CLASSIFIER = os.path.join(
    DATABASE_DIR, "classification", "SILVA", SILVA_CLASSIFIER_FILENAME
)
SILVA_SPECIES_REFERENCE_FILENAME = "silva_v144_assignSpecies.fa.gz"
SILVA_SPECIES_REFERENCE_URL = (
    "https://www.arb-silva.de/fileadmin/silva_databases/release_144/"
    "DADA2/1.36.0/SSU/" + SILVA_SPECIES_REFERENCE_FILENAME
)
SILVA_SPECIES_REFERENCE_MD5 = "444de7c0cce0b66addda7a3f8b38e012"
SILVA_SPECIES_REFERENCE = os.path.join(
    DATABASE_DIR,
    "classification",
    "SILVA",
    SILVA_SPECIES_REFERENCE_FILENAME,
)
PR2_CLASSIFIER = os.path.join(
    DATABASE_DIR,
    "classification",
    "PR2",
    "pr2_version_5.1.1_SSU_dada2.clean.culled.derep-sliced_"
    + config["fwdPrimer"]
    + "_"
    + config["revPrimer"]
    + "_dereplicated_final_classifier_qiime2-2026.7.qza",
)

required_preexisting_database_resources = []
if USE_PREEXISTING_BBSPLIT_DATABASE:
    required_preexisting_database_resources.append(
        ("BBsplit", BBSPLIT_DB_DIR)
    )
if USE_PREEXISTING_SILVA_DATABASE:
    required_preexisting_database_resources.extend(
        [
            ("SILVA classifier", SILVA_CLASSIFIER),
            ("SILVA species reference", SILVA_SPECIES_REFERENCE),
        ]
    )
if USE_PREEXISTING_PR2_DATABASE:
    required_preexisting_database_resources.append(
        ("PR2 classifier", PR2_CLASSIFIER)
    )

if required_preexisting_database_resources:
    missing_database_resources = [
        (label, path)
        for label, path in required_preexisting_database_resources
        if not os.path.exists(path)
    ]
    if missing_database_resources:
        missing_list = "\n  - ".join(
            f"{label}: {path}" for label, path in missing_database_resources
        )
        raise WorkflowError(
            "One or more databases are marked as pre-existing, but the "
            "following required resources were not found:\n  - "
            + missing_list
            + "\nCorrect database_dir or mark only the missing database "
            "family as not pre-existing so the workflow prepares it."
        )

# Small, user-facing files copied into a project-specific folder at the end of
# the run. Older configs without projectName use studyName.
PROJECT_NAME = str(config.get("projectName", config["studyName"])).strip()
if not re.fullmatch(r"[A-Za-z0-9._-]+", PROJECT_NAME):
    raise WorkflowError(
        "projectName may contain only letters, numbers, periods, underscores, "
        "and hyphens"
    )
RESULTS_EXPORT_DIR = PROJECT_NAME + "-Results-Export"
RESULTS_LONG_DATA = (
    ISD_FILTERED_LONG_TABLE
    if USE_INTERNAL_STANDARDS
    else "results/04-formatted/" + config["studyName"] + ".long_data.tsv"
)
RESULTS_EXPORT_INPUTS = [
    RESULTS_LONG_DATA,
    "results/04-formatted/" + config["studyName"] + ".asv_sequences.tsv",
    "results/07-report/" + config["studyName"] + ".pipeline-report.html",
    "results/02-proks/09-subsetting/tax-merged/"
    + config["studyName"]
    + ".PR2-plastid-routing-audit.tsv",
    "results/02-proks/09-subsetting/tax-merged/"
    + config["studyName"]
    + ".PR2-plastid-routing-summary.tsv",
]
RESULTS_EXPORT_OUTPUTS = [
    RESULTS_EXPORT_DIR + "/" + config["studyName"] + ".long_data.tsv",
    RESULTS_EXPORT_DIR + "/" + config["studyName"] + ".asv_sequences.tsv",
    RESULTS_EXPORT_DIR + "/" + config["studyName"] + ".pipeline-report.html",
    RESULTS_EXPORT_DIR + "/" + config["studyName"] + ".PR2-plastid-routing-audit.tsv",
    RESULTS_EXPORT_DIR + "/" + config["studyName"] + ".PR2-plastid-routing-summary.tsv",
]
RESULTS_EXPORT_SUMMARIES = [
    RESULTS_EXPORT_DIR + "/" + config["studyName"] + ".phylum_summary.tsv",
    RESULTS_EXPORT_DIR + "/" + config["studyName"] + ".order_summary.tsv",
]
if USE_INTERNAL_STANDARDS:
    RESULTS_EXPORT_INPUTS.append(ISD_CORRECTED_LONG_TABLE)
    RESULTS_EXPORT_OUTPUTS.append(
        RESULTS_EXPORT_DIR + "/"
        + config["studyName"]
        + "."
        + ALL_INTERNAL_STANDARD_IDS_STEM
        + ".ISD_corrected_long_data.tsv"
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
    final_output.append("results/02-proks/09-subsetting/reclassified/all_16S_ASVs_PR2.classified.qza"),
    final_output.append("results/02-proks/09-subsetting/tax-merged/chloroplasts-PR2-reclassified-merged-classification.qza"),
    final_output.append("results/02-proks/09-subsetting/tax-merged/" + config["studyName"] + ".PR2-plastid-routing-audit.tsv"),
    final_output.append("results/02-proks/09-subsetting/tax-merged/" + config["studyName"] + ".PR2-plastid-routing-summary.tsv"),
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
        final_output.append(ISD_CORRECTED_LONG_TABLE)

    final_output.extend(RESULTS_EXPORT_OUTPUTS)
    final_output.extend(RESULTS_EXPORT_SUMMARIES)
    final_output.append(RESULTS_EXPORT_DIR + "/README.txt")

    return final_output

# validate sample sheet and config file
validate(samples, schema="../../config/schemas/samples.schema.yml")
validate(config, schema="../../config/schemas/config.schema.yml")
