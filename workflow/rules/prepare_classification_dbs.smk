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

rule clean_pr2_fasta_extract_headers:
    input:
        DATABASE_PREFIX + "classification/PR2/pr2_version_5.1.1_SSU_dada2.fasta"
    output:
        clean=temp(DATABASE_PREFIX + "classification/PR2/pr2_version_5.1.1_SSU_dada2.clean.fasta"),
        headers=DATABASE_PREFIX + "classification/PR2/pr2_version_5.1.1_SSU_dada2.headers.txt"
    priority: 50
    script:
        "../scripts/reformat_pr2_reference.py"

rule import_pr2_fasta:
    input:
        DATABASE_PREFIX + "classification/PR2/pr2_version_5.1.1_SSU_dada2.clean.fasta"
    output:
        temp(DATABASE_PREFIX + "classification/PR2/pr2_version_5.1.1_SSU_dada2.clean.qza")
    priority: 49
    conda:
        config["qiime2version"]
    shell:
        "qiime tools import --type 'FeatureData[Sequence]' "
        "--input-path {input:q} --output-path {output:q}"

rule import_pr2_taxonomy:
    input:
        DATABASE_PREFIX + "classification/PR2/pr2_version_5.1.1_SSU_dada2.headers.txt"
    output:
        temp(DATABASE_PREFIX + "classification/PR2/pr2_version_5.1.1_SSU_dada2.headers.qza")
    priority: 49
    conda:
        config["qiime2version"]
    shell:
        "qiime tools import --type 'FeatureData[Taxonomy]' "
        "--input-format HeaderlessTSVTaxonomyFormat "
        "--input-path {input:q} --output-path {output:q}"

rule cull_pr2_seqs:
    input:
        DATABASE_PREFIX + "classification/PR2/pr2_version_5.1.1_SSU_dada2.clean.qza"
    output:
        temp(DATABASE_PREFIX + "classification/PR2/pr2_version_5.1.1_SSU_dada2.clean.culled.qza")
    priority: 48
    conda:
        config["qiime2version"]
    shell:
        "qiime rescript cull-seqs --i-sequences {input:q} "
        "--o-clean-sequences {output:q}"

rule derep_seqs_taxonomy:
    input:
        culled=DATABASE_PREFIX + "classification/PR2/pr2_version_5.1.1_SSU_dada2.clean.culled.qza",
        taxonomy=DATABASE_PREFIX + "classification/PR2/pr2_version_5.1.1_SSU_dada2.headers.qza"
    output:
        derepseqs=temp(DATABASE_PREFIX + "classification/PR2/pr2_version_5.1.1_SSU_dada2.clean.culled.derep.qza"),
        dereptaxa=temp(DATABASE_PREFIX + "classification/PR2/pr2_version_5.1.1_SSU_dada2.headers.derep.qza")
    priority: 47
    conda:
        config["qiime2version"]
    shell:
        "qiime rescript dereplicate --i-sequences {input.culled:q} "
        "--i-taxa {input.taxonomy:q} --p-mode uniq --p-rank-handles disable "
        "--o-dereplicated-sequences {output.derepseqs:q} "
        "--o-dereplicated-taxa {output.dereptaxa:q}"

rule extract_primers_pr2:
    input:
        DATABASE_PREFIX + "classification/PR2/pr2_version_5.1.1_SSU_dada2.clean.culled.derep.qza"
    params:
        fwdPrimer=config["fwdPrimer"],
        revPrimer=config["revPrimer"]
    output:
        slicedDNA=temp(DATABASE_PREFIX + "classification/PR2/pr2_version_5.1.1_SSU_dada2.clean.culled.derep-sliced_" + config["fwdPrimer"] + "_" + config["revPrimer"] + ".qza")
    log:
        "logs/PR2_classification_db_prep_qc_PR2_seqs_extract_primers.log"
    priority: 45
    conda:
        config["qiime2version"]
    threads: 8
    shell:
        "qiime feature-classifier extract-reads --i-sequences {input:q} "
        "--p-f-primer {params.fwdPrimer:q} --p-r-primer {params.revPrimer:q} "
        "--p-n-jobs {threads} --p-read-orientation forward "
        "--o-reads {output.slicedDNA:q} 2> {log:q}"

rule dereplicate_extracted_pr2_reads:
    input:
        slicedDNA=rules.extract_primers_pr2.output.slicedDNA,
        dereplicatedTaxa=DATABASE_PREFIX + "classification/PR2/pr2_version_5.1.1_SSU_dada2.headers.derep.qza"
    output:
        slicedDNAdereplicated=temp(DATABASE_PREFIX + "classification/PR2/pr2_version_5.1.1_SSU_dada2.clean.culled.derep-sliced_" + config["fwdPrimer"] + "_" + config["revPrimer"] + "_derep.qza"),
        dereplicatedTaxaSliced=temp(DATABASE_PREFIX + "classification/PR2/pr2_version_5.1.1_SSU_dada2.headers.derep_" + config["fwdPrimer"] + "_" + config["revPrimer"] + "_derep.qza")
    log:
        "logs/PR2_classification_db_prep_qc_PR2_seqs_dereplicate_sliced_data.log"
    priority: 44
    conda:
        config["qiime2version"]
    shell:
        "qiime rescript dereplicate --i-sequences {input.slicedDNA:q} "
        "--i-taxa {input.dereplicatedTaxa:q} --p-mode uniq "
        "--o-dereplicated-sequences {output.slicedDNAdereplicated:q} "
        "--o-dereplicated-taxa {output.dereplicatedTaxaSliced:q} 2> {log:q}"

rule train_classifier_pr2:
    input:
        slicedDNAdereplicated=rules.dereplicate_extracted_pr2_reads.output.slicedDNAdereplicated,
        dereplicatedTaxaSliced=rules.dereplicate_extracted_pr2_reads.output.dereplicatedTaxaSliced   
    output:
        PR2_CLASSIFIER
    log:
        "logs/PR2_classification_db_prep_qc_PR2_seqs_train_sliced_classifier.log"
    priority: 43
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
        "--o-classifier {output:q} 2> {log:q}"
 
