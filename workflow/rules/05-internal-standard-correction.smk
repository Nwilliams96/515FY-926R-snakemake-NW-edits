if USE_INTERNAL_STANDARDS:

    INTERNAL_STANDARD_FASTAS = expand(
        "config/intstd_fastas/{standard}.fasta",
        standard=INTERNAL_STANDARD_IDS,
    )
    INTERNAL_STANDARD_ASV_LISTS = expand(
        "config/intstd_fastas/{standard}.asvs.txt",
        standard=INTERNAL_STANDARD_IDS,
    )
    INTERNAL_STANDARD_METHOD_TABLES = expand(
        "results/05-internal-std-corrected/" + config["studyName"] + ".asv_table_{method}.tsv",
        method=INTERNAL_STANDARD_METHOD_STEMS,
    )


    rule prepare_internal_standard_fastas:
        input:
            table="config/internal_stds.tsv"
        output:
            fastas=INTERNAL_STANDARD_FASTAS
        params:
            standard_ids=INTERNAL_STANDARD_IDS
        script:
            "../scripts/prepare_internal_standard_fastas.py"


    rule prepare_16S_BLASTdb:
        input:
            fasta="config/intstd_fastas/{standard}.fasta"
        output:
            nhr="config/intstd_fastas/{standard}.fasta.nhr"
        wildcard_constraints:
            standard="|".join(re.escape(value) for value in INTERNAL_STANDARD_IDS)
        conda:
            "../envs/blast-env.yaml"
        shell:
            "makeblastdb -dbtype nucl -in {input.fasta:q} && test -s {output.nhr:q}"


    rule identify_intsd_ASVS:
        input:
            latestseqs="results/02-proks/04-DADA2d-plaintext-exports/" + config["studyName"] + ".16S.latest_seqs.fasta",
            database="config/intstd_fastas/{standard}.fasta.nhr"
        output:
            matches="config/intstd_fastas/{standard}.asvs.outfmt6.tsv",
            asvs="config/intstd_fastas/{standard}.asvs.txt"
        wildcard_constraints:
            standard="|".join(re.escape(value) for value in INTERNAL_STANDARD_IDS)
        params:
            database=lambda wildcards: f"config/intstd_fastas/{wildcards.standard}.fasta"
        conda:
            "../envs/blast-env.yaml"
        shell:
            "blastn -query {input.latestseqs:q} -db {params.database:q} "
            "-outfmt 6 -perc_identity 99 -qcov_hsp_perc 100 > {output.matches:q} "
            "&& cut -f1 {output.matches:q} > {output.asvs:q}"


    rule intstd_correct_data:
        input:
            asv_table="results/04-formatted/" + config["studyName"] + ".long_data.tsv",
            standard_asvs=INTERNAL_STANDARD_ASV_LISTS,
            isd="config/internal_stds.tsv",
            isd_added="config/samples.tsv"
        params:
            standard_ids=INTERNAL_STANDARD_IDS,
            standard_slots=INTERNAL_STANDARD_SLOTS,
            method_stems=INTERNAL_STANDARD_METHOD_STEMS
        output:
            corrected=ISD_CORRECTED_LONG_TABLE,
            filtered=ISD_FILTERED_LONG_TABLE,
            method_tables=INTERNAL_STANDARD_METHOD_TABLES,
            recovery_plot="results/06-figures/" + config["studyName"] + ".recovery_ratios.pdf",
            domain_plot="results/06-figures/" + config["studyName"] + ".Domain_by_sampleID.pdf",
            recovery_plot_png="results/06-figures/" + config["studyName"] + ".recovery_ratios.png",
            domain_plot_png="results/06-figures/" + config["studyName"] + ".Domain_by_sampleID.png"
        conda:
            "../envs/r-tidyverse-2.0.0.yml"
        log:
            "logs/05-internal-std-correction/prepare_long_data_corrected.log"
        script:
            "../scripts/correct-w-isds-snakemake-v5.R"
