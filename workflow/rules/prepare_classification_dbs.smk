#set up rules with different scripts for rescript pipeline as implemented in bash, may need to split out
#make variables for everything that can have variables - db version, primer sequence, etc

rule download_SILVA:
    input:
        rules.initialize_database_directories.output.marker
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
    shell:
        "qiime rescript get-silva-data --p-version {params.SILVAversion:q} "
        "--p-target SSURef_NR99 --p-include-species-labels "
        "--o-silva-sequences {output.seqs:q} --o-silva-taxonomy {output.taxonomy:q} "
        "2> {log:q}"

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
    shell:
        "qiime rescript reverse-transcribe --i-rna-sequences {input:q} "
        "--o-dna-sequences {output:q} 2> {log:q}"

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
    shell:
        "qiime rescript cull-seqs --i-sequences {input.rawDNA:q} "
        "--o-clean-sequences {output.cleanDNA:q} 2> {log:q}"

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
    shell:
        "qiime rescript filter-seqs-length-by-taxon "
        "--i-sequences {input.cleanDNA:q} --i-taxonomy {input.taxonomy:q} "
        "--p-labels Archaea Bacteria Eukaryota --p-min-lens 900 1200 1400 "
        "--o-filtered-seqs {output.filteredDNA:q} "
        "--o-discarded-seqs {output.discardedDNA:q} 2> {log:q}"

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
    shell:
        "qiime rescript dereplicate --i-sequences {input.filteredDNA:q} "
        "--i-taxa {input.taxonomy:q} --p-mode uniq "
        "--o-dereplicated-sequences {output.dereplicatedDNA:q} "
        "--o-dereplicated-taxa {output.dereplicatedTaxa:q} 2> {log:q}"

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
    threads: 8
    shell:
        "qiime feature-classifier extract-reads "
        "--i-sequences {input.dereplicatedDNA:q} "
        "--p-f-primer {params.fwdPrimer:q} --p-r-primer {params.revPrimer:q} "
        "--p-n-jobs {threads} --p-read-orientation forward "
        "--o-reads {output.slicedDNA:q} 2> {log:q}"

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
    shell:
        "qiime rescript dereplicate --i-sequences {input.slicedDNA:q} "
        "--i-taxa {input.dereplicatedTaxa:q} --p-mode uniq "
        "--o-dereplicated-sequences {output.slicedDNAdereplicated:q} "
        "--o-dereplicated-taxa {output.dereplicatedTaxaSliced:q} 2> {log:q}"

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
    threads: 1
    resources:
        mem_mb=32000,
        runtime=360,
    shell:
        "qiime feature-classifier fit-classifier-naive-bayes "
        "--i-reference-reads {input.slicedDNAdereplicated:q} "
        "--i-reference-taxonomy {input.dereplicatedTaxaSliced:q} "
        "--o-classifier {output:q} 2> {log:q}"

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
        DATABASE_PREFIX + "classification/PR2/pr2_version_5.1.1_SSU_dada2.clean.culled.derep-sliced_" + config["fwdPrimer"] + "_" + config["revPrimer"] + "_dereplicated_final_classifier_USE_ME.qza"
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
 
