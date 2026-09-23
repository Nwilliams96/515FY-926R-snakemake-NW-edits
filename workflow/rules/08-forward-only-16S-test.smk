"""Optional forward-read-only 16S sensitivity analysis.

Run with the dedicated launcher or request the ``forward_only_16S_test``
target explicitly.  This branch reuses the primer-trimmed, BBsplit-assigned R1
FASTQs and canonical paired-end outputs; it is deliberately excluded from the
main workflow's default target.
"""


FORWARD_ONLY_ROOT = "results/08-forward-only-16S-test"
FORWARD_ONLY_DADA2 = FORWARD_ONLY_ROOT + "/02-DADA2d"
FORWARD_ONLY_TAXONOMY = FORWARD_ONLY_ROOT + "/03-taxonomy"
FORWARD_ONLY_ANALYSIS = FORWARD_ONLY_ROOT + "/04-comparison"
FORWARD_ONLY_REFERENCE = FORWARD_ONLY_ROOT + "/05-reference-lengths"
FORWARD_ONLY_EXPORT = PROJECT_NAME + "-Forward-Only-16S-Test"


rule forward_only_16S_test:
    input:
        FORWARD_ONLY_EXPORT + "/README.txt",
        FORWARD_ONLY_EXPORT + "/forward_only_summary.tsv",
        FORWARD_ONLY_EXPORT + "/forward_only_ASV_comparison.tsv",
        FORWARD_ONLY_EXPORT + "/taxon_presence_comparison.tsv",
        FORWARD_ONLY_EXPORT + "/insert_length_estimates.tsv",
        FORWARD_ONLY_EXPORT + "/insert_length_summary.tsv",
    default_target: False


rule create_forward_only_16S_manifest:
    input:
        r1=expand(
            "results/01-split/{sample}.prok.R1.fastq.gz",
            sample=samples["sample"],
        )
    output:
        FORWARD_ONLY_ROOT + "/01-import/manifest.tsv"
    params:
        sample_ids=list(samples["sample"])
    conda:
        config["qiime2version"]
    script:
        "../scripts/create_forward_only_manifest.py"


rule import_forward_only_16S:
    input:
        manifest=rules.create_forward_only_16S_manifest.output
    output:
        sequences=FORWARD_ONLY_ROOT + "/01-import/forward_reads.qza"
    conda:
        config["qiime2version"]
    log:
        "logs/08-forward-only-16S-test/01-import.log"
    script:
        "../scripts/import_forward_only_16S.sh"


rule denoise_forward_only_16S:
    input:
        sequences=rules.import_forward_only_16S.output.sequences
    output:
        directory(FORWARD_ONLY_DADA2),
        sequences=FORWARD_ONLY_DADA2 + "/representative_sequences.qza",
        stats=FORWARD_ONLY_DADA2 + "/denoising_stats.qza",
        table=FORWARD_ONLY_DADA2 + "/table.qza",
        transitions=FORWARD_ONLY_DADA2 + "/base_transition_stats.qza",
    params:
        trunc_len=DADA2_PROK["trunc_len_f"],
        max_ee=DADA2_PROK["max_ee_f"],
        trunc_q=DADA2_PROK["trunc_q"],
        pooling_method=DADA2_PROK["pooling_method"],
        chimera_method=DADA2_PROK["chimera_method"],
        min_fold_parent_over_abundance=DADA2_PROK[
            "min_fold_parent_over_abundance"
        ],
        n_reads_learn=DADA2_PROK["n_reads_learn"],
    threads: 8
    resources:
        mem_mb=96000,
    conda:
        config["qiime2version"]
    log:
        "logs/08-forward-only-16S-test/02-DADA2.log"
    script:
        "../scripts/denoise_forward_only_16S.sh"


rule classify_forward_only_16S_with_SILVA:
    input:
        sequences=rules.denoise_forward_only_16S.output.sequences,
        classDB=SILVA_CLASSIFIER,
    output:
        classified=FORWARD_ONLY_TAXONOMY + "/SILVA.classified.qza"
    params:
        classDB=""
    threads: 8
    resources:
        mem_mb=64000,
    conda:
        config["qiime2version"]
    log:
        "logs/08-forward-only-16S-test/03-SILVA.log"
    script:
        "../scripts/P05-classify-eASVs.sh"


rule classify_forward_only_16S_with_PR2:
    input:
        sequences=rules.denoise_forward_only_16S.output.sequences,
        database_ready=rules.ensure_pr2_classifier.output.marker,
    output:
        classified=FORWARD_ONLY_TAXONOMY + "/PR2.classified.qza"
    params:
        classDB=PR2_CLASSIFIER
    threads: 8
    resources:
        mem_mb=64000,
    conda:
        config["qiime2version"]
    log:
        "logs/08-forward-only-16S-test/03-PR2.log"
    script:
        "../scripts/P05-classify-eASVs.sh"


rule resolve_forward_only_PR2_plastids:
    input:
        silva_taxonomy=rules.classify_forward_only_16S_with_SILVA.output.classified,
        pr2_taxonomy=rules.classify_forward_only_16S_with_PR2.output.classified,
    output:
        taxonomy=FORWARD_ONLY_TAXONOMY + "/resolved_taxonomy.qza",
        audit=FORWARD_ONLY_TAXONOMY + "/PR2-plastid-routing-audit.tsv",
        summary=FORWARD_ONLY_TAXONOMY + "/PR2-plastid-routing-summary.tsv",
    params:
        min_confidence=CHLOROPLAST_PR2_MIN_CONFIDENCE
    conda:
        config["qiime2version"]
    log:
        "logs/08-forward-only-16S-test/03-resolve-PR2.log"
    script:
        "../scripts/resolve_pr2_plastids.py"


rule compare_forward_only_16S_to_paired:
    input:
        forward_table=rules.denoise_forward_only_16S.output.table,
        forward_sequences=rules.denoise_forward_only_16S.output.sequences,
        forward_stats=rules.denoise_forward_only_16S.output.stats,
        forward_taxonomy=rules.resolve_forward_only_PR2_plastids.output.taxonomy,
        forward_silva_taxonomy=rules.classify_forward_only_16S_with_SILVA.output.classified,
        paired_table=rules.denoise_prok_dada2.output.proktable,
        paired_sequences=rules.denoise_prok_dada2.output.prokrepseqs,
        paired_taxonomy=rules.resolve_PR2_plastids.output.taxonomy,
        paired_silva_taxonomy=rules.classify_ASVs.output.classified,
    output:
        asv_comparison=FORWARD_ONLY_ANALYSIS + "/forward_only_ASV_comparison.tsv",
        taxon_comparison=FORWARD_ONLY_ANALYSIS + "/taxon_presence_comparison.tsv",
        unique_fasta=FORWARD_ONLY_ANALYSIS + "/forward_only_unique_ASVs.fasta",
        summary=FORWARD_ONLY_ANALYSIS + "/forward_only_summary.tsv",
        forward_table=FORWARD_ONLY_ANALYSIS + "/forward_only_feature_table.tsv",
        forward_sequences=FORWARD_ONLY_ANALYSIS + "/forward_only_sequences.fasta",
        forward_taxonomy=FORWARD_ONLY_ANALYSIS + "/forward_only_taxonomy.tsv",
        forward_stats=FORWARD_ONLY_ANALYSIS + "/forward_only_denoising_stats.tsv",
    conda:
        config["qiime2version"]
    log:
        "logs/08-forward-only-16S-test/04-compare.log"
    script:
        "../scripts/compare_forward_only_16S.py"


rule blast_forward_only_unique_ASVs:
    input:
        query=rules.compare_forward_only_16S_to_paired.output.unique_fasta,
        reference=SILVA_SPECIES_REFERENCE,
    output:
        hits=FORWARD_ONLY_REFERENCE + "/SILVA_species_blast_hits.tsv"
    threads: 8
    resources:
        mem_mb=32000,
    conda:
        config["qiime2version"]
    log:
        "logs/08-forward-only-16S-test/05-reference-BLAST.log"
    script:
        "../scripts/blast_forward_only_16S.sh"


rule estimate_forward_only_insert_lengths:
    input:
        hits=rules.blast_forward_only_unique_ASVs.output.hits,
        reference=SILVA_SPECIES_REFERENCE,
        comparison=rules.compare_forward_only_16S_to_paired.output.asv_comparison,
    output:
        estimates=FORWARD_ONLY_REFERENCE + "/insert_length_estimates.tsv",
        summary=FORWARD_ONLY_REFERENCE + "/insert_length_summary.tsv",
    params:
        forward_primer=config["fwdPrimer"],
        reverse_primer=config["revPrimer"],
        trunc_len_f=DADA2_PROK["trunc_len_f"],
        trunc_len_r=DADA2_PROK["trunc_len_r"],
        min_overlap=DADA2_PROK["min_overlap"],
    conda:
        config["qiime2version"]
    script:
        "../scripts/estimate_forward_only_insert_lengths.py"


rule export_forward_only_16S_test:
    input:
        asv_comparison=rules.compare_forward_only_16S_to_paired.output.asv_comparison,
        taxon_comparison=rules.compare_forward_only_16S_to_paired.output.taxon_comparison,
        unique_fasta=rules.compare_forward_only_16S_to_paired.output.unique_fasta,
        summary=rules.compare_forward_only_16S_to_paired.output.summary,
        forward_table=rules.compare_forward_only_16S_to_paired.output.forward_table,
        forward_sequences=rules.compare_forward_only_16S_to_paired.output.forward_sequences,
        forward_taxonomy=rules.compare_forward_only_16S_to_paired.output.forward_taxonomy,
        forward_stats=rules.compare_forward_only_16S_to_paired.output.forward_stats,
        pr2_audit=rules.resolve_forward_only_PR2_plastids.output.audit,
        pr2_summary=rules.resolve_forward_only_PR2_plastids.output.summary,
        blast_hits=rules.blast_forward_only_unique_ASVs.output.hits,
        estimates=rules.estimate_forward_only_insert_lengths.output.estimates,
        length_summary=rules.estimate_forward_only_insert_lengths.output.summary,
    output:
        directory(FORWARD_ONLY_EXPORT),
        readme=FORWARD_ONLY_EXPORT + "/README.txt",
        asv_comparison=FORWARD_ONLY_EXPORT + "/forward_only_ASV_comparison.tsv",
        taxon_comparison=FORWARD_ONLY_EXPORT + "/taxon_presence_comparison.tsv",
        estimates=FORWARD_ONLY_EXPORT + "/insert_length_estimates.tsv",
        summary=FORWARD_ONLY_EXPORT + "/forward_only_summary.tsv",
        length_summary=FORWARD_ONLY_EXPORT + "/insert_length_summary.tsv",
    params:
        study_name=config["studyName"],
        trunc_len_f=DADA2_PROK["trunc_len_f"],
        trunc_len_r=DADA2_PROK["trunc_len_r"],
        min_overlap=DADA2_PROK["min_overlap"],
    script:
        "../scripts/export_forward_only_16S_test.py"
