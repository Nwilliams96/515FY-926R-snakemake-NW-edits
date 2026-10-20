# AI-assisted change — three database choices, human review pending.
# Purpose: preserve named scientific steps while selecting local reuse,
# classifier download, or RESCRIPt training independently for each database.
# Configuration/path resolution and publication helpers contain no training.

def classifier_source(record):
    if record["mode"] == "existing" or database_config.reusable_classifier(
        record["classifier"], record["checksum"]
    ):
        return record["classifier"]
    return record["staged"]


# A missing local build must pass its environment's native version check
# before acquisition. No SILVA-144 adapter or automatic training fallback.
if SILVA_DB["mode"] == "train":
    rule check_silva_reference_capability:
        output:
            os.path.join(SILVA_BUILD_ROOT, "capability.json")
        params:
            family="SILVA",
            version=SILVA_VERSION
        conda:
            config["qiime2version"]
        script:
            "../scripts/check_reference_capability.py"

if PR2_DB["mode"] == "train":
    rule check_pr2_reference_capability:
        output:
            os.path.join(PR2_BUILD_ROOT, "capability.json")
        params:
            family="PR2",
            version=PR2_VERSION
        conda:
            config["qiime2version"]
        script:
            "../scripts/check_reference_capability.py"


if SILVA_DB["mode"] == "download":
    rule download_SILVA_classifier:
        output:
            temp(SILVA_DB["staged"])
        params:
            url=SILVA_DB["url"],
            md5=SILVA_DB["checksum"],
            validate_classifier=True
        conda:
            config["qiime2version"]
        log:
            "logs/SILVA_classifier_download.log"
        script:
            "../scripts/download_verified_reference.py"

if PR2_DB["mode"] == "download":
    rule download_PR2_classifier:
        output:
            temp(PR2_DB["staged"])
        params:
            url=PR2_DB["url"],
            md5=PR2_DB["checksum"],
            validate_classifier=True
        conda:
            config["qiime2version"]
        log:
            "logs/PR2_classifier_download.log"
        script:
            "../scripts/download_verified_reference.py"


# This supporting reference is independent of how the classifier was obtained.
# In particular, an existing classifier may live in a read-only directory.
if RUN_SILVA_SPECIES_MATCHING:
    rule download_SILVA_144_species_reference:
        input:
            rules.initialize_database_directories.output.marker
        output:
            SILVA_SPECIES_REFERENCE
        params:
            url=SILVA_SPECIES_REFERENCE_URL,
            md5=SILVA_SPECIES_REFERENCE_MD5
        log:
            "logs/SILVA_144_species_reference_download.log"
        script:
            "../scripts/download_verified_reference.py"


if SILVA_DB["mode"] == "train":
    # Original SILVA stages restored, with private intermediate paths.
    SILVA_FILES = {
        name: os.path.join(SILVA_BUILD_ROOT, name + ".qza")
        for name in ("rna", "taxonomy", "dna", "culled", "filtered", "discarded",
                     "derep", "derep-taxonomy", "sliced", "extraction-stats",
                     "sliced-derep", "sliced-taxonomy")
    }
    SILVA_RANKS = ["domain"] + (["kingdom"] if SILVA_VERSION == "144" else []) + [
        "phylum", "class", "order", "family", "genus"
    ]

    rule download_SILVA:
        input:
            rules.check_silva_reference_capability.output
        output:
            seqs=temp(SILVA_FILES["rna"]),
            taxonomy=temp(SILVA_FILES["taxonomy"])
        params:
            version=SILVA_VERSION,
            ranks=SILVA_RANKS
        conda:
            config["qiime2version"]
        log:
            "logs/SILVA_classification_db_prep/01-download-reference.log"
        shell:
            "qiime rescript get-silva-data --p-version {params.version:q} "
            "--p-target SSURef_NR99 --p-ranks {params.ranks:q} "
            "--p-rank-propagation --p-no-include-species-labels "
            "--o-silva-sequences {output.seqs:q} "
            "--o-silva-taxonomy {output.taxonomy:q} > {log:q} 2>&1"

    rule reverse_transcribe:
        input:
            rules.download_SILVA.output.seqs
        output:
            temp(SILVA_FILES["dna"])
        conda:
            config["qiime2version"]
        log:
            "logs/SILVA_classification_db_prep/02-reverse-transcribe.log"
        shell:
            "qiime rescript reverse-transcribe --i-rna-sequences {input:q} "
            "--o-dna-sequences {output:q} > {log:q} 2>&1"

    rule qc_seqs_cull:
        input:
            rules.reverse_transcribe.output
        output:
            temp(SILVA_FILES["culled"])
        conda:
            config["qiime2version"]
        log:
            "logs/SILVA_classification_db_prep/03-cull.log"
        shell:
            "qiime rescript cull-seqs --i-sequences {input:q} "
            "--p-num-degenerates 5 --p-homopolymer-length 8 "
            "--o-clean-sequences {output:q} > {log:q} 2>&1"

    # The original helper was missing from the imported repository. These
    # explicit selected cutoffs are documented, not claimed as recovered code.
    rule qc_seqs_filter:
        input:
            sequences=rules.qc_seqs_cull.output,
            taxonomy=rules.download_SILVA.output.taxonomy
        output:
            filteredDNA=temp(SILVA_FILES["filtered"]),
            discardedDNA=temp(SILVA_FILES["discarded"])
        conda:
            config["qiime2version"]
        log:
            "logs/SILVA_classification_db_prep/04-filter-length.log"
        shell:
            "qiime rescript filter-seqs-length-by-taxon "
            "--i-sequences {input.sequences:q} --i-taxonomy {input.taxonomy:q} "
            "--p-labels Bacteria Archaea Eukaryota --p-min-lens 1200 900 1400 "
            "--o-filtered-seqs {output.filteredDNA:q} "
            "--o-discarded-seqs {output.discardedDNA:q} > {log:q} 2>&1"

    rule qc_seqs_dereplicate:
        input:
            sequences=rules.qc_seqs_filter.output.filteredDNA,
            taxonomy=rules.download_SILVA.output.taxonomy
        output:
            sequences=temp(SILVA_FILES["derep"]),
            taxonomy=temp(SILVA_FILES["derep-taxonomy"])
        params:
            ranks=SILVA_RANKS
        conda:
            config["qiime2version"]
        log:
            "logs/SILVA_classification_db_prep/05-dereplicate-full-reference.log"
        shell:
            "qiime rescript dereplicate --i-sequences {input.sequences:q} "
            "--i-taxa {input.taxonomy:q} --p-mode uniq --p-rank-handles {params.ranks:q} "
            "--o-dereplicated-sequences {output.sequences:q} "
            "--o-dereplicated-taxa {output.taxonomy:q} > {log:q} 2>&1"

    rule extract_primers:
        input:
            rules.qc_seqs_dereplicate.output.sequences
        output:
            slicedDNA=temp(SILVA_FILES["sliced"]),
            extractionStats=temp(SILVA_FILES["extraction-stats"])
        params:
            fwdPrimer=config["fwdPrimer"],
            revPrimer=config["revPrimer"]
        conda:
            config["qiime2version"]
        threads: 8
        log:
            "logs/SILVA_classification_db_prep/06-extract-primer-region.log"
        shell:
            "qiime feature-classifier extract-reads --i-sequences {input:q} "
            "--p-f-primer {params.fwdPrimer:q} --p-r-primer {params.revPrimer:q} "
            "--p-n-jobs {threads} --p-read-orientation forward "
            "--o-reads {output.slicedDNA:q} "
            "--o-read-extraction-stats {output.extractionStats:q} > {log:q} 2>&1"

    rule dereplicated_sliced_data:
        input:
            sequences=rules.extract_primers.output.slicedDNA,
            taxonomy=rules.qc_seqs_dereplicate.output.taxonomy
        output:
            sequences=temp(SILVA_FILES["sliced-derep"]),
            taxonomy=temp(SILVA_FILES["sliced-taxonomy"])
        params:
            ranks=SILVA_RANKS
        conda:
            config["qiime2version"]
        log:
            "logs/SILVA_classification_db_prep/07-dereplicate-primer-region.log"
        shell:
            "qiime rescript dereplicate --i-sequences {input.sequences:q} "
            "--i-taxa {input.taxonomy:q} --p-mode uniq --p-rank-handles {params.ranks:q} "
            "--o-dereplicated-sequences {output.sequences:q} "
            "--o-dereplicated-taxa {output.taxonomy:q} > {log:q} 2>&1"

    rule train_classifier:
        input:
            sequences=(SILVA_FILES["sliced-derep"] if config["silva_reference_mode"] == "amplicon" else SILVA_FILES["derep"]),
            taxonomy=(SILVA_FILES["sliced-taxonomy"] if config["silva_reference_mode"] == "amplicon" else SILVA_FILES["derep-taxonomy"])
        output:
            temp(SILVA_DB["staged"])
        conda:
            config["qiime2version"]
        resources:
            mem_mb=110000,
            runtime=720
        log:
            "logs/SILVA_classification_db_prep/08-train-classifier.log"
        shell:
            "qiime feature-classifier fit-classifier-naive-bayes "
            "--i-reference-reads {input.sequences:q} "
            "--i-reference-taxonomy {input.taxonomy:q} "
            "--o-classifier {output:q} > {log:q} 2>&1"


if PR2_DB["mode"] == "train":
    PR2_CLEAN_QZA = os.path.join(PR2_BUILD_ROOT, "sequences.qza")
    PR2_RAW_TAXONOMY_QZA = os.path.join(PR2_BUILD_ROOT, "raw-taxonomy.qza")
    PR2_TAXONOMY_QZA = os.path.join(PR2_BUILD_ROOT, "taxonomy.qza")
    PR2_CULLED_QZA = os.path.join(PR2_BUILD_ROOT, "culled.qza")
    PR2_DEREP_QZA = os.path.join(PR2_BUILD_ROOT, "derep.qza")
    PR2_DEREP_TAXONOMY_QZA = os.path.join(PR2_BUILD_ROOT, "derep-taxonomy.qza")
    PR2_SLICED_QZA = os.path.join(PR2_BUILD_ROOT, "sliced.qza")
    PR2_EXTRACTION_STATS_QZA = os.path.join(PR2_BUILD_ROOT, "extraction-stats.qza")
    PR2_SLICED_DEREP_QZA = os.path.join(PR2_BUILD_ROOT, "sliced-derep.qza")
    PR2_SLICED_TAXONOMY_QZA = os.path.join(PR2_BUILD_ROOT, "sliced-taxonomy.qza")
    PR2_STAGED_CLASSIFIER = PR2_DB["staged"]

    rule download_pr2_reference:
        input:
            rules.check_pr2_reference_capability.output
        output:
            sequences=temp(PR2_CLEAN_QZA),
            taxonomy=temp(PR2_RAW_TAXONOMY_QZA)
        params:
            version=PR2_VERSION
        conda:
            config["qiime2version"]
        log:
            "logs/PR2_classification_db_prep/00-download-reference.log"
        shell:
            "qiime rescript get-pr2-data --p-version {params.version:q} "
            "--p-ranks domain supergroup division subdivision class order family genus species "
            "--o-pr2-sequences {output.sequences:q} "
            "--o-pr2-taxonomy {output.taxonomy:q} > {log:q} 2>&1"

    # RESCRIPt prefixes PR2 ranks. Strip only those prefixes, preserving all
    # nine ranks and :plas. Otherwise downstream code would mistake d__ for SILVA.
    rule normalize_pr2_reference_taxonomy:
        input:
            rules.download_pr2_reference.output.taxonomy
        output:
            temp(PR2_TAXONOMY_QZA)
        conda:
            config["qiime2version"]
        log:
            "logs/PR2_classification_db_prep/01-normalize-taxonomy.log"
        shell:
            r"""
            qiime rescript edit-taxonomy --i-taxonomy {input:q} \
              --p-search-strings '(^|;)(?:d__|sgr__|dv__|dvs__|c__|o__|f__|g__|s__)' \
              --p-replacement-strings '\1' --p-use-regex \
              --o-edited-taxonomy {output:q} > {log:q} 2>&1
            """

    rule cull_pr2_seqs:
        input:
            PR2_CLEAN_QZA
        output:
            temp(PR2_CULLED_QZA)
        log:
            "logs/PR2_classification_db_prep/05-cull-sequences.log"
        priority: 46
        conda:
            config["qiime2version"]
        shell:
            "qiime rescript cull-seqs --i-sequences {input:q} "
            "--p-num-degenerates 5 --p-homopolymer-length 8 "
            "--o-clean-sequences {output:q} > {log:q} 2>&1"


    rule derep_seqs_taxonomy:
        input:
            culled=rules.cull_pr2_seqs.output,
            taxonomy=PR2_TAXONOMY_QZA
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
            "--i-taxa {input.dereplicatedTaxa:q} --p-mode uniq --p-rank-handles disable "
            "--o-dereplicated-sequences {output.slicedDNAdereplicated:q} "
            "--o-dereplicated-taxa {output.dereplicatedTaxaSliced:q} > {log:q} 2>&1"


    rule train_classifier_pr2:
        input:
            slicedDNAdereplicated=(PR2_SLICED_DEREP_QZA if config["pr2_reference_mode"] == "amplicon" else PR2_DEREP_QZA),
            dereplicatedTaxaSliced=(PR2_SLICED_TAXONOMY_QZA if config["pr2_reference_mode"] == "amplicon" else PR2_DEREP_TAXONOMY_QZA)
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




# Publication is infrastructure only: scientific commands remain above.
# Existing classifiers are validated without writing to their source directory.
rule ensure_silva_classifier:
    input:
        classifier=lambda wc: classifier_source(SILVA_DB)
    output:
        marker=SILVA_READY_MARKER
    params:
        classifier=SILVA_CLASSIFIER,
        lock=SILVA_DB["lock"],
        read_only=SILVA_DB["mode"] == "existing",
        md5=SILVA_DB["checksum"],
        label="SILVA"
    conda:
        config["qiime2version"]
    log:
        "logs/SILVA_classification_db_prep/validate-and-publish.log"
    script:
        "../scripts/publish_pr2_classifier.py"


rule ensure_pr2_classifier:
    input:
        classifier=lambda wc: classifier_source(PR2_DB)
    output:
        marker=PR2_READY_MARKER
    params:
        classifier=PR2_CLASSIFIER,
        lock=PR2_BUILD_LOCK,
        read_only=PR2_DB["mode"] == "existing",
        md5=PR2_DB["checksum"],
        label="PR2"
    conda:
        config["qiime2version"]
    log:
        "logs/PR2_classification_db_prep/10-validate-and-publish.log"
    resources:
        mem_mb=2000,
        runtime=30
    script:
        "../scripts/publish_pr2_classifier.py"
