#set up rules with different scripts for rescript pipeline as implemented in bash, may need to split out
#make variables for everything that can have variables - db version, primer sequence, etc

rule download_SILVA:
    output:
        seqs=temp(DATABASE_PREFIX + "classification/SILVA/silva-ssu-nr99-rna-seqs.qza"),
        taxonomy=temp(DATABASE_PREFIX + "classification/SILVA/silva-ssu-nr99-tax.qza")
    params:
        SILVAversion=config["SILVAversion"]
    log:
        "logs/SILVA_classification_db_prep_download.log"
    priority: 50
    conda:
        config["qiime2version"]
    script:
        "../scripts/tax-classifier-construction/SILVA/01-rescript-dl-database-file.sh"

rule reverse_transcribe:
    input:
        DATABASE_PREFIX + "classification/SILVA/silva-ssu-nr99-rna-seqs.qza"
    output:
        temp(DATABASE_PREFIX + "classification/SILVA/silva-ssu-nr99-dna-seqs.qza")
    log:
        "logs/SILVA_classification_db_prep_reverse_transcribe.log"
    priority: 49
    conda:
        config["qiime2version"]
    script:
        "../scripts/tax-classifier-construction/SILVA/02-reverse-transcribe.sh"

rule qc_seqs_cull:
    input:
        rawDNA=DATABASE_PREFIX + "classification/SILVA/silva-ssu-nr99-dna-seqs.qza"
    output:
        cleanDNA=temp(DATABASE_PREFIX + "classification/SILVA/silva-ssu-nr99-dna-seqs-culled.qza")
    log:
        "logs/SILVA_classification_db_prep_qc_SILVA_seqs_cull.log"
    priority: 48
    conda:
        config["qiime2version"]
    script:
        "../scripts/tax-classifier-construction/SILVA/03-qc-seqs-cull.sh"

rule qc_seqs_filter:
    input:
        cleanDNA=DATABASE_PREFIX + "classification/SILVA/silva-ssu-nr99-dna-seqs-culled.qza",
        taxonomy=DATABASE_PREFIX + "classification/SILVA/silva-ssu-nr99-tax.qza"
    output:
        filteredDNA=temp(DATABASE_PREFIX + "classification/SILVA/silva-ssu-nr99-dna-seqs-culled-filtered.qza"),
        discardedDNA=temp(DATABASE_PREFIX + "classification/SILVA/silva-ssu-nr99-dna-seqs-culled-discarded.qza")
    log:
        "logs/SILVA_classification_db_prep_qc_SILVA_seqs_filter.log"
    priority: 47
    conda:
        config["qiime2version"]
    script:
        "../scripts/tax-classifier-construction/SILVA/04-qc-seqs-length-filter.sh"

rule qc_seqs_dereplicate:
    input:
        filteredDNA=DATABASE_PREFIX + "classification/SILVA/silva-ssu-nr99-dna-seqs-culled-filtered.qza",
        taxonomy=DATABASE_PREFIX + "classification/SILVA/silva-ssu-nr99-tax.qza"
    output:
        dereplicatedDNA=temp(DATABASE_PREFIX + "classification/SILVA/silva-ssu-nr99-dna-seqs-culled-filtered-dereplicated.qza"),
        dereplicatedTaxa=temp(DATABASE_PREFIX + "classification/SILVA/silva-ssu-nr99-tax-dereplicated.qza")
    log:
        "logs/SILVA_classification_db_prep_qc_SILVA_seqs_dereplicate.log"
    priority: 46
    conda:
        config["qiime2version"]
    script:
        "../scripts/tax-classifier-construction/SILVA/05-qc-seqs-dereplicate.sh"

rule extract_primers:
    input:
        dereplicatedDNA=DATABASE_PREFIX + "classification/SILVA/silva-ssu-nr99-dna-seqs-culled-filtered-dereplicated.qza"
    params:
        fwdPrimer=config["fwdPrimer"],
        revPrimer=config["revPrimer"]
    output:
        slicedDNA=temp(DATABASE_PREFIX + "classification/SILVA/silva-ssu-nr99-tax-dereplicated-sliced_" + config["fwdPrimer"] + "_" + config["revPrimer"] + ".qza")
    log:
        "logs/SILVA_classification_db_prep_qc_SILVA_seqs_extract_primers.log"
    priority: 45
    conda:
        config["qiime2version"]
    script:
        "../scripts/tax-classifier-construction/SILVA/06-extract-primers.sh"

rule dereplicated_sliced_data:
    input:
        slicedDNA=rules.extract_primers.output.slicedDNA,
        dereplicatedTaxa=DATABASE_PREFIX + "classification/SILVA/silva-ssu-nr99-tax-dereplicated.qza"
    output:
        slicedDNAdereplicated=temp(DATABASE_PREFIX + "classification/SILVA/silva-ssu-nr99-tax-dereplicated-sliced_" + config["fwdPrimer"] + "_" + config["revPrimer"] + "_dereplicated.qza"),
        dereplicatedTaxaSliced=temp(DATABASE_PREFIX + "classification/SILVA/silva-ssu-nr99-tax-dereplicated_" + config["fwdPrimer"] + "_" + config["revPrimer"] + "_dereplicated.qza")
    log:
        "logs/SILVA_classification_db_prep_qc_SILVA_seqs_dereplicate_sliced_data.log"
    priority: 44
    conda:
        config["qiime2version"]
    script:
        "../scripts/tax-classifier-construction/SILVA/07-deduplicate-sliced-data.sh"

rule train_classifier:
    input:
        slicedDNAdereplicated=rules.dereplicated_sliced_data.output.slicedDNAdereplicated,
        dereplicatedTaxaSliced=rules.dereplicated_sliced_data.output.dereplicatedTaxaSliced
    output:
        DATABASE_PREFIX + "classification/SILVA/silva-ssu-nr99-tax-dereplicated-sliced_" + config["fwdPrimer"] + "_" + config["revPrimer"] + "_dereplicated_final_classifier_USE_ME.qza"
    log:
        "logs/SILVA_classification_db_prep_qc_SILVA_seqs_train_sliced_classifier.log"
    priority: 43
    conda:
        config["qiime2version"]
    script:
        "../scripts/tax-classifier-construction/SILVA/08-train-classifier.sh"

rule clean_pr2_fasta_extract_headers:
    input:
        DATABASE_PREFIX + "classification/PR2/pr2_version_5.1.1_SSU_dada2.fasta"
    output:
        clean=temp(DATABASE_PREFIX + "classification/PR2/pr2_version_5.1.1_SSU_dada2.clean.fasta"),
        headers=DATABASE_PREFIX + "classification/PR2/pr2_version_5.1.1_SSU_dada2.headers.txt"
    priority: 50
    script:
        "../scripts/tax-classifier-construction/PR2/script_to_reformat_PR.sh"

rule import_pr2_fasta:
    input:
        DATABASE_PREFIX + "classification/PR2/pr2_version_5.1.1_SSU_dada2.clean.fasta"
    output:
        temp(DATABASE_PREFIX + "classification/PR2/pr2_version_5.1.1_SSU_dada2.clean.qza")
    priority: 49
    conda:
        config["qiime2version"]
    script:
        "../scripts/tax-classifier-construction/PR2/00-import-fasta.sh"

rule import_pr2_taxonomy:
    input:
        DATABASE_PREFIX + "classification/PR2/pr2_version_5.1.1_SSU_dada2.headers.txt"
    output:
        temp(DATABASE_PREFIX + "classification/PR2/pr2_version_5.1.1_SSU_dada2.headers.qza")
    priority: 49
    conda:
        config["qiime2version"]
    script:
        "../scripts/tax-classifier-construction/PR2/01-import-headers.sh"

rule cull_pr2_seqs:
    input:
        DATABASE_PREFIX + "classification/PR2/pr2_version_5.1.1_SSU_dada2.clean.qza"
    output:
        temp(DATABASE_PREFIX + "classification/PR2/pr2_version_5.1.1_SSU_dada2.clean.culled.qza")
    priority: 48
    conda:
        config["qiime2version"]
    script:
        "../scripts/tax-classifier-construction/PR2/02-qc-seqs-cull.sh"

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
    script:
        "../scripts/tax-classifier-construction/PR2/03-qc-seqs-dereplicate.sh"

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
    script:
        "../scripts/tax-classifier-construction/PR2/04-extract-primers.sh"

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
    script:
        "../scripts/tax-classifier-construction/PR2/05-deduplicate-sliced-data.sh"

rule train_classifier_pr2:
    input:
        slicedDNAdereplicated=rules.dereplicate_extracted_pr2_reads.output.slicedDNAdereplicated,
        dereplicatedTaxaSliced=rules.dereplicate_extracted_pr2_reads.output.dereplicatedTaxaSliced   
    output:
        DATABASE_PREFIX + "classification/PR2/pr2_version_5.1.1_SSU_dada2.clean.culled.derep-sliced_" + config["fwdPrimer"] + "_" + config["revPrimer"] + "_dereplicated_final_classifier_USE_ME.qza"
    log:
        "logs/PR2_classification_db_prep_qc_PR2_seqs_train_sliced_classifier.log"
    priority: 43
    conda:
        config["qiime2version"]
    script:
        "../scripts/tax-classifier-construction/PR2/06-train-classifier.sh"
 
