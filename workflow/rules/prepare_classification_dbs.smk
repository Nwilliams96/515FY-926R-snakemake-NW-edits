import zipfile


# AI-assisted change — see CODE_PROVENANCE.md, commits 86ca36d and d230025.
# Purpose: use the release-pinned SILVA 144 resources and publish downloads
# safely when independent projects share one database directory.
# Human review status: pending.

# SILVA 144 is supplied as an official, cryptographically signed QIIME 2
# classifier. The matching 515F/926R classifier is used for this pipeline's
# default primers; other primer pairs use SILVA's full-length SSU classifier.
rule download_SILVA_144_classifier:
    input:
        rules.initialize_database_directories.output.marker
    output:
        SILVA_CLASSIFIER
    params:
        url=SILVA_CLASSIFIER_URL,
        md5=SILVA_CLASSIFIER_MD5,
    log:
        "logs/SILVA_144_classifier_download.log"
    priority: 50
    script:
        "../scripts/download_verified_reference.py"


# The signed QIIME classifier above intentionally stops at Genus. Species are
# added separately only for unambiguous, exact ASV matches to SILVA's official
# DADA2 species reference.
rule download_SILVA_144_species_reference:
    input:
        rules.initialize_database_directories.output.marker
    output:
        SILVA_SPECIES_REFERENCE
    params:
        url=SILVA_SPECIES_REFERENCE_URL,
        md5=SILVA_SPECIES_REFERENCE_MD5,
    log:
        "logs/SILVA_144_species_reference_download.log"
    priority: 50
    script:
        "../scripts/download_verified_reference.py"


# AI-assisted change — see CODE_PROVENANCE.md, commits 8d63a81 and 1a25e17.
# Purpose: retain the original rule-per-step PR2 construction workflow while
# keeping all intermediate files private to this project. Only the completed,
# validated classifier is published to the shared database directory.
# Human review status: pending.

PR2_BUILD_ROOT = os.path.join(
    ".snakemake",
    "database-build",
    "PR2",
    "pr2-5.1.1-" + config["fwdPrimer"] + "-" + config["revPrimer"],
)
PR2_SOURCE_GZ = os.path.join(PR2_BUILD_ROOT, "pr2_version_5.1.1_SSU_dada2.fasta.gz")
PR2_SOURCE_FASTA = os.path.join(PR2_BUILD_ROOT, "pr2_version_5.1.1_SSU_dada2.fasta")
PR2_CLEAN_FASTA = os.path.join(PR2_BUILD_ROOT, "pr2_version_5.1.1_SSU_dada2.clean.fasta")
PR2_TAXONOMY_TSV = os.path.join(PR2_BUILD_ROOT, "pr2_version_5.1.1_SSU_dada2.taxonomy.tsv")
PR2_CLEAN_QZA = os.path.join(PR2_BUILD_ROOT, "pr2_version_5.1.1_SSU_dada2.clean.qza")
PR2_TAXONOMY_QZA = os.path.join(PR2_BUILD_ROOT, "pr2_version_5.1.1_SSU_dada2.taxonomy.qza")
PR2_CULLED_QZA = os.path.join(PR2_BUILD_ROOT, "pr2_version_5.1.1_SSU_dada2.clean.culled.qza")
PR2_DEREP_QZA = os.path.join(PR2_BUILD_ROOT, "pr2_version_5.1.1_SSU_dada2.clean.culled.derep.qza")
PR2_DEREP_TAXONOMY_QZA = os.path.join(PR2_BUILD_ROOT, "pr2_version_5.1.1_SSU_dada2.taxonomy.derep.qza")
PR2_SLICED_QZA = os.path.join(PR2_BUILD_ROOT, "pr2_version_5.1.1_SSU_dada2.sliced.qza")
PR2_EXTRACTION_STATS_QZA = os.path.join(PR2_BUILD_ROOT, "pr2_version_5.1.1_SSU_dada2.read-extraction-stats.qza")
PR2_SLICED_DEREP_QZA = os.path.join(PR2_BUILD_ROOT, "pr2_version_5.1.1_SSU_dada2.sliced.derep.qza")
PR2_SLICED_TAXONOMY_QZA = os.path.join(PR2_BUILD_ROOT, "pr2_version_5.1.1_SSU_dada2.sliced.taxonomy.qza")
PR2_STAGED_CLASSIFIER = os.path.join(PR2_BUILD_ROOT, "pr2_version_5.1.1.classifier.qza")


def reusable_pr2_classifier_archive(path):
    """Cheap DAG-time integrity check; final QIIME validation remains a rule."""
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return False
    try:
        with zipfile.ZipFile(path) as artifact:
            metadata_files = [
                name for name in artifact.namelist() if name.endswith("/metadata.yaml")
            ]
            if len(metadata_files) != 1 or artifact.testzip() is not None:
                return False
            metadata = artifact.read(metadata_files[0]).decode("utf-8", errors="replace")
            return "TaxonomicClassifier" in metadata
    except (OSError, zipfile.BadZipFile):
        return False


def existing_or_staged_pr2_classifier(wildcards):
    # A user who explicitly selected a pre-existing database receives a clear
    # validation failure rather than an unrequested rebuild. Otherwise, a
    # missing/corrupt archive schedules the visible construction rules below.
    if USE_PREEXISTING_PR2_DATABASE or reusable_pr2_classifier_archive(PR2_CLASSIFIER):
        return PR2_CLASSIFIER
    return PR2_STAGED_CLASSIFIER


rule download_pr2_reference:
    input:
        rules.initialize_database_directories.output.marker
    output:
        temp(PR2_SOURCE_GZ)
    params:
        url="https://github.com/pr2database/pr2database/releases/download/v5.1.1/pr2_version_5.1.1_SSU_dada2.fasta.gz",
        checksum="0c8728abcbb2126eed2c7e587f820cbce39c138cdfdb51239bbf18621498462d",
        checksum_algorithm="sha256",
    log:
        "logs/PR2_classification_db_prep/00-download-reference.log"
    priority: 50
    script:
        "../scripts/download_verified_reference.py"


rule decompress_pr2_reference:
    input:
        PR2_SOURCE_GZ
    output:
        temp(PR2_SOURCE_FASTA)
    log:
        "logs/PR2_classification_db_prep/01-decompress-reference.log"
    priority: 49
    shell:
        "gzip -dc {input:q} > {output:q} 2> {log:q}"


rule clean_pr2_fasta_extract_headers:
    input:
        PR2_SOURCE_FASTA
    output:
        clean=temp(PR2_CLEAN_FASTA),
        headers=temp(PR2_TAXONOMY_TSV)
    log:
        "logs/PR2_classification_db_prep/02-reformat-reference.log"
    priority: 48
    script:
        "../scripts/reformat_pr2_reference.py"


rule import_pr2_fasta:
    input:
        rules.clean_pr2_fasta_extract_headers.output.clean
    output:
        temp(PR2_CLEAN_QZA)
    log:
        "logs/PR2_classification_db_prep/03-import-sequences.log"
    priority: 47
    conda:
        config["qiime2version"]
    shell:
        "qiime tools import --type 'FeatureData[Sequence]' "
        "--input-path {input:q} --output-path {output:q} > {log:q} 2>&1"


rule import_pr2_taxonomy:
    input:
        rules.clean_pr2_fasta_extract_headers.output.headers
    output:
        temp(PR2_TAXONOMY_QZA)
    log:
        "logs/PR2_classification_db_prep/04-import-taxonomy.log"
    priority: 47
    conda:
        config["qiime2version"]
    shell:
        "qiime tools import --type 'FeatureData[Taxonomy]' "
        "--input-format HeaderlessTSVTaxonomyFormat "
        "--input-path {input:q} --output-path {output:q} > {log:q} 2>&1"


rule cull_pr2_seqs:
    input:
        rules.import_pr2_fasta.output
    output:
        temp(PR2_CULLED_QZA)
    log:
        "logs/PR2_classification_db_prep/05-cull-sequences.log"
    priority: 46
    conda:
        config["qiime2version"]
    shell:
        "qiime rescript cull-seqs --i-sequences {input:q} "
        "--o-clean-sequences {output:q} > {log:q} 2>&1"


rule derep_seqs_taxonomy:
    input:
        culled=rules.cull_pr2_seqs.output,
        taxonomy=rules.import_pr2_taxonomy.output
    output:
        derepseqs=temp(PR2_DEREP_QZA),
        dereptaxa=temp(PR2_DEREP_TAXONOMY_QZA)
    log:
        "logs/PR2_classification_db_prep/06-dereplicate-full-reference.log"
    priority: 45
    conda:
        config["qiime2version"]
    shell:
        "qiime rescript dereplicate --i-sequences {input.culled:q} "
        "--i-taxa {input.taxonomy:q} --p-mode uniq --p-rank-handles disable "
        "--o-dereplicated-sequences {output.derepseqs:q} "
        "--o-dereplicated-taxa {output.dereptaxa:q} > {log:q} 2>&1"


rule extract_primers_pr2:
    input:
        rules.derep_seqs_taxonomy.output.derepseqs
    output:
        slicedDNA=temp(PR2_SLICED_QZA),
        extractionStats=temp(PR2_EXTRACTION_STATS_QZA)
    params:
        fwdPrimer=config["fwdPrimer"],
        revPrimer=config["revPrimer"]
    log:
        "logs/PR2_classification_db_prep/07-extract-primer-region.log"
    priority: 44
    conda:
        config["qiime2version"]
    threads: 8
    shell:
        "qiime feature-classifier extract-reads --i-sequences {input:q} "
        "--p-f-primer {params.fwdPrimer:q} --p-r-primer {params.revPrimer:q} "
        "--p-n-jobs {threads} --p-read-orientation forward "
        "--o-reads {output.slicedDNA:q} "
        "--o-read-extraction-stats {output.extractionStats:q} > {log:q} 2>&1"


rule dereplicate_extracted_pr2_reads:
    input:
        slicedDNA=rules.extract_primers_pr2.output.slicedDNA,
        dereplicatedTaxa=rules.derep_seqs_taxonomy.output.dereptaxa
    output:
        slicedDNAdereplicated=temp(PR2_SLICED_DEREP_QZA),
        dereplicatedTaxaSliced=temp(PR2_SLICED_TAXONOMY_QZA)
    log:
        "logs/PR2_classification_db_prep/08-dereplicate-primer-region.log"
    priority: 43
    conda:
        config["qiime2version"]
    shell:
        "qiime rescript dereplicate --i-sequences {input.slicedDNA:q} "
        "--i-taxa {input.dereplicatedTaxa:q} --p-mode uniq "
        "--o-dereplicated-sequences {output.slicedDNAdereplicated:q} "
        "--o-dereplicated-taxa {output.dereplicatedTaxaSliced:q} > {log:q} 2>&1"


rule train_classifier_pr2:
    input:
        slicedDNAdereplicated=rules.dereplicate_extracted_pr2_reads.output.slicedDNAdereplicated,
        dereplicatedTaxaSliced=rules.dereplicate_extracted_pr2_reads.output.dereplicatedTaxaSliced
    output:
        temp(PR2_STAGED_CLASSIFIER)
    log:
        "logs/PR2_classification_db_prep/09-train-classifier.log"
    priority: 42
    conda:
        config["qiime2version"]
    threads: 1
    resources:
        mem_mb=32000,
        runtime=360,
    shell:
        "qiime feature-classifier fit-classifier-naive-bayes "
        "--i-reference-reads {input.slicedDNAdereplicated:q} "
        "--i-reference-taxonomy {input.dereplicatedTaxaSliced:q} "
        "--o-classifier {output:q} > {log:q} 2>&1"


rule ensure_pr2_classifier:
    input:
        classifier=existing_or_staged_pr2_classifier
    output:
        marker=PR2_READY_MARKER
    params:
        classifier=PR2_CLASSIFIER,
        lock=PR2_BUILD_LOCK,
    log:
        "logs/PR2_classification_db_prep/10-validate-and-publish.log"
    priority: 41
    conda:
        config["qiime2version"]
    resources:
        mem_mb=2000,
        runtime=30,
    script:
        "../scripts/publish_pr2_classifier.py"
