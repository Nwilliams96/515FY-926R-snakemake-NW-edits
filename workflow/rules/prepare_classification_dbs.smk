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

rule ensure_pr2_classifier:
    input:
        rules.initialize_database_directories.output.marker
    output:
        marker=PR2_READY_MARKER
    params:
        classifier=PR2_CLASSIFIER,
        lock=PR2_BUILD_LOCK,
        source_url="https://github.com/pr2database/pr2database/releases/download/v5.1.1/pr2_version_5.1.1_SSU_dada2.fasta.gz",
        source_checksum="0c8728abcbb2126eed2c7e587f820cbce39c138cdfdb51239bbf18621498462d",
        fwd_primer=config["fwdPrimer"],
        rev_primer=config["revPrimer"],
        allow_build=not USE_PREEXISTING_PR2_DATABASE,
    log:
        "logs/PR2_classification_db_prep.log"
    priority: 50
    conda:
        config["qiime2version"]
    threads: 8
    resources:
        mem_mb=32000,
        runtime=360,
    script:
        "../scripts/ensure_pr2_classifier.py"
 
