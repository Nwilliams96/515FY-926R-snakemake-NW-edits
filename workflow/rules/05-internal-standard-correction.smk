if config["use_internal_standards"]:

    rule prepare_16S_BLASTdb:
        input:
            intstd1="config/intstd_fastas/" + config["intstds"]["intstd1"] + ".fasta",
            intstd2="config/intstd_fastas/" + config["intstds"]["intstd2"] + ".fasta",
            intstd3="config/intstd_fastas/" + config["intstds"]["intstd3"] + ".fasta"
        output:
            intstd1="config/intstd_fastas/" + config["intstds"]["intstd1"] + ".fasta.nhr",
            intstd2="config/intstd_fastas/" + config["intstds"]["intstd2"] + ".fasta.nhr",
            intstd3="config/intstd_fastas/" + config["intstds"]["intstd3"] + ".fasta.nhr"
        conda:
            "../envs/blast-env.yaml"
        script:
            "../scripts/B00-makeblastdb.sh"

    rule identify_intsd_ASVS:
        input:
            latestseqs="results/02-proks/04-DADA2d-plaintext-exports/" + config["studyName"] + ".16S.latest_seqs.fasta",
            intstd1="config/intstd_fastas/" + config["intstds"]["intstd1"] + ".fasta.nhr",
            intstd2="config/intstd_fastas/" + config["intstds"]["intstd2"] + ".fasta.nhr",
            intstd3="config/intstd_fastas/" + config["intstds"]["intstd3"] + ".fasta.nhr"
        output:
            intstd1tsv="config/intstd_fastas/" + config["intstds"]["intstd1"] + ".asvs.outfmt6.tsv",
            intstd1asvs="config/intstd_fastas/" + config["intstds"]["intstd1"] + ".asvs.txt",
            intstd2tsv="config/intstd_fastas/" + config["intstds"]["intstd2"] + ".asvs.outfmt6.tsv",
            intstd2asvs="config/intstd_fastas/" + config["intstds"]["intstd2"] + ".asvs.txt",
            intstd3tsv="config/intstd_fastas/" + config["intstds"]["intstd3"] + ".asvs.outfmt6.tsv",
            intstd3asvs="config/intstd_fastas/" + config["intstds"]["intstd3"] + ".asvs.txt"
        conda:
            "../envs/blast-env.yaml"
        script:
            "../scripts/B01-blast-ASVs.sh"

    rule intstd_correct_data:
        input:
            asv_table="results/04-formatted/" + config["studyName"] + ".long_data.tsv",
            intstd1asvs="config/intstd_fastas/" + config["intstds"]["intstd1"] + ".asvs.txt",
            intstd2asvs="config/intstd_fastas/" + config["intstds"]["intstd2"] + ".asvs.txt",
            intstd3asvs="config/intstd_fastas/" + config["intstds"]["intstd3"] + ".asvs.txt",
            isd="config/internal_stds.tsv",
            isd_added="config/samples.tsv"
        params:
            intstd1name=config["intstds"]["intstd1"],
            intstd2name=config["intstds"]["intstd2"],
            intstd3name=config["intstds"]["intstd3"],
        output:
            corrected="results/05-internal-std-corrected/" + config["studyName"] + ".ISD_corrected_asv_table.tsv",
            isd_1_recovery_ratio="results/05-internal-std-corrected/" + config["studyName"] + ".asv_table_isd_1_recovery_ratio.tsv",
            isd_2_recovery_ratio="results/05-internal-std-corrected/" + config["studyName"] + ".asv_table_isd_2_recovery_ratio.tsv",
            isd_3_recovery_ratio="results/05-internal-std-corrected/" + config["studyName"] + ".asv_table_isd_3_recovery_ratio.tsv",
            mean_recovery_ratio="results/05-internal-std-corrected/" + config["studyName"] + ".asv_table_mean_recovery_ratio.tsv",
            median_recovery_ratio="results/05-internal-std-corrected/" + config["studyName"] + ".asv_table_median_recovery_ratio.tsv",
            isd_1_isd_2_mean_recovery_ratio="results/05-internal-std-corrected/" + config["studyName"] + ".asv_table_isd_1_isd_2_mean_recovery_ratio.tsv",
            isd_1_isd_3_mean_recovery_ratio="results/05-internal-std-corrected/" + config["studyName"] + ".asv_table_isd_1_isd_3_mean_recovery_ratio.tsv",
            isd_2_isd_3_mean_recovery_ratio="results/05-internal-std-corrected/" + config["studyName"] + ".asv_table_isd_2_isd_3_mean_recovery_ratio.tsv",
            recovery_plot="results/06-figures/" + config["studyName"] + ".recovery_ratios.pdf",
            domain_plot="results/06-figures/" + config["studyName"] + ".Domain_by_sampleID.pdf"
        conda:
            "../envs/r-tidyverse-2.0.0.yml"
        log:
            "logs/05-internal-std-correction/prepare_long_data_corrected.log"
        script:
            "../scripts/correct-w-isds-snakemake-v5.R"
