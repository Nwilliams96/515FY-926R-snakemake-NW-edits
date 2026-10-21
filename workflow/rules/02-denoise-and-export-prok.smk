rule create_manifest_prok:
    input:
        r1=expand("results/01-split/{sample}.prok.R1.fastq.gz", sample=samples["sample"]),
        r2=expand("results/01-split/{sample}.prok.R2.fastq.gz", sample=samples["sample"])
    output:
        "results/02-proks/manifest.tsv"
    conda:
        config["qiime2version"]
    script:
        "../scripts/P00-create-manifest.sh"

rule import_prok:
    input:
        "results/02-proks/manifest.tsv"
    output:
        "results/02-proks/16S.qza"
    conda:
        config["qiime2version"]
    script:
        "../scripts/P01-import.sh"

rule visualize_prok_seq_quality:
    input:
        "results/02-proks/16S.qza"
    output:
        directory("results/02-proks/02-quality-plots-R1-R2/")
    conda:
        config["qiime2version"]
    script:
        "../scripts/P02-visualize-quality_R1-R2.sh"

rule denoise_prok_dada2:
    input:
        "results/02-proks/16S.qza"
    params:
        truncR1=DADA2_PROK["trunc_len_f"],
        truncR2=DADA2_PROK["trunc_len_r"],
        max_ee_f=DADA2_PROK["max_ee_f"],
        max_ee_r=DADA2_PROK["max_ee_r"],
        trunc_q=DADA2_PROK["trunc_q"],
        min_overlap=DADA2_PROK["min_overlap"],
        max_merge_mismatch=DADA2_PROK["max_merge_mismatch"],
        trim_overhang=DADA2_PROK["trim_overhang"],
        retain_unmerged=False,
        pooling_method=DADA2_PROK["pooling_method"],
        chimera_method=DADA2_PROK["chimera_method"],
        min_fold_parent_over_abundance=DADA2_PROK["min_fold_parent_over_abundance"],
        n_reads_learn=DADA2_PROK["n_reads_learn"]
    output:
        directory("results/02-proks/03-DADA2d/"),
        prokrepseqs="results/02-proks/03-DADA2d/representative_sequences.qza",
        prokstats="results/02-proks/03-DADA2d/denoising_stats.qza",
        proktable="results/02-proks/03-DADA2d/table.qza",
        prokbasetransitions="results/02-proks/03-DADA2d/base_transition_stats.qza"
    conda:
        config["qiime2version"]
    threads: 8
    resources:
        mem_mb=96000,
    log:
        "logs/02-denoise-and-export-prok/03-DADA2/DADA2.stderrout"
    script:
        "../scripts/P03-DADA2.sh"

rule denoise_prok_dada2_with_unmerged:
    input:
        "results/02-proks/16S.qza"
    params:
        truncR1=DADA2_PROK["trunc_len_f"],
        truncR2=DADA2_PROK["trunc_len_r"],
        max_ee_f=DADA2_PROK["max_ee_f"],
        max_ee_r=DADA2_PROK["max_ee_r"],
        trunc_q=DADA2_PROK["trunc_q"],
        min_overlap=DADA2_PROK["min_overlap"],
        max_merge_mismatch=DADA2_PROK["max_merge_mismatch"],
        trim_overhang=DADA2_PROK["trim_overhang"],
        retain_unmerged=True,
        pooling_method=DADA2_PROK["pooling_method"],
        chimera_method=DADA2_PROK["chimera_method"],
        min_fold_parent_over_abundance=DADA2_PROK["min_fold_parent_over_abundance"],
        n_reads_learn=DADA2_PROK["n_reads_learn"]
    output:
        directory("results/02-proks/03-DADA2d-with-unmerged/"),
        prokrepseqs="results/02-proks/03-DADA2d-with-unmerged/representative_sequences.qza",
        prokstats="results/02-proks/03-DADA2d-with-unmerged/denoising_stats.qza",
        proktable="results/02-proks/03-DADA2d-with-unmerged/table.qza",
        prokbasetransitions="results/02-proks/03-DADA2d-with-unmerged/base_transition_stats.qza"
    conda:
        config["qiime2version"]
    threads: 8
    resources:
        mem_mb=96000,
    log:
        "logs/02-denoise-and-export-prok/03_DADA2_with_unmerged/DADA2.stderrout"
    script:
        "../scripts/P03-DADA2.sh"

rule classify_unmerged_16S_with_SILVA:
    input:
        sequences=rules.denoise_prok_dada2_with_unmerged.output.prokrepseqs,
        # A local readiness marker schedules download/training/validation while
        # keeping the shared classifier outside Snakemake's output cleanup.
        database_ready=rules.ensure_silva_classifier.output.marker,
    params:
        classDB=SILVA_CLASSIFIER,
    output:
        classified="results/02-proks/11-unmerged-16S-rescue/SILVA.classified.qza"
    conda:
        config["qiime2version"]
    threads: 8
    resources:
        mem_mb=64000,
    log:
        "logs/02-denoise-and-export-prok/11-unmerged-SILVA.log"
    script:
        "../scripts/P05-classify-eASVs.sh"

rule classify_unmerged_16S_with_PR2:
    input:
        sequences=rules.denoise_prok_dada2_with_unmerged.output.prokrepseqs,
        database_ready=rules.ensure_pr2_classifier.output.marker,
    params:
        classDB=PR2_CLASSIFIER,
    output:
        classified="results/02-proks/11-unmerged-16S-rescue/PR2.classified.qza"
    conda:
        config["qiime2version"]
    # Each sklearn worker loads its own classifier state. Two workers keep the
    # large unmerged-ASV classification below the CARC job memory limit.
    threads: 2
    resources:
        mem_mb=110000,
    log:
        "logs/02-denoise-and-export-prok/11-unmerged-PR2.log"
    script:
        "../scripts/P05-classify-eASVs.sh"


# AI-assisted change — see provenance ledger linked from README.md, commit f13649a.
# Purpose: expose every QIIME export and BIOM conversion used by the optional
# linked/unmerged 16S supplement as a separately inspectable DAG rule.
# Human review status: pending.
rule export_unmerged_16S_retained_table:
    input:
        rules.denoise_prok_dada2_with_unmerged.output.proktable
    output:
        temp(".snakemake/classification-exports/unmerged-16S/retained-table/feature-table.biom")
    params:
        export_dir=".snakemake/classification-exports/unmerged-16S/retained-table"
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-prok/11-export-retained-table.log"
    shell:
        "mkdir -p {params.export_dir:q} && "
        "qiime tools export --input-path {input:q} --output-path {params.export_dir:q} "
        "> {log:q} 2>&1"


rule export_unmerged_16S_standard_table:
    input:
        rules.denoise_prok_dada2.output.proktable
    output:
        temp(".snakemake/classification-exports/unmerged-16S/standard-table/feature-table.biom")
    params:
        export_dir=".snakemake/classification-exports/unmerged-16S/standard-table"
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-prok/11-export-standard-table.log"
    shell:
        "mkdir -p {params.export_dir:q} && "
        "qiime tools export --input-path {input:q} --output-path {params.export_dir:q} "
        "> {log:q} 2>&1"


rule export_unmerged_16S_sequences:
    input:
        rules.denoise_prok_dada2_with_unmerged.output.prokrepseqs
    output:
        temp(".snakemake/classification-exports/unmerged-16S/sequences/linked-dna-sequences.fasta")
    params:
        export_dir=".snakemake/classification-exports/unmerged-16S/sequences"
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-prok/11-export-sequences.log"
    shell:
        "mkdir -p {params.export_dir:q} && "
        "qiime tools export --input-path {input:q} --output-path {params.export_dir:q} "
        "> {log:q} 2>&1"


rule export_unmerged_16S_SILVA_taxonomy:
    input:
        rules.classify_unmerged_16S_with_SILVA.output.classified
    output:
        temp(".snakemake/classification-exports/unmerged-16S/SILVA/taxonomy.tsv")
    params:
        export_dir=".snakemake/classification-exports/unmerged-16S/SILVA"
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-prok/11-export-SILVA-taxonomy.log"
    shell:
        "mkdir -p {params.export_dir:q} && "
        "qiime tools export --input-path {input:q} --output-path {params.export_dir:q} "
        "> {log:q} 2>&1"


rule export_unmerged_16S_PR2_taxonomy:
    input:
        rules.classify_unmerged_16S_with_PR2.output.classified
    output:
        temp(".snakemake/classification-exports/unmerged-16S/PR2-raw/taxonomy.tsv")
    params:
        export_dir=".snakemake/classification-exports/unmerged-16S/PR2-raw"
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-prok/11-export-PR2-taxonomy.log"
    shell:
        "mkdir -p {params.export_dir:q} && "
        "qiime tools export --input-path {input:q} --output-path {params.export_dir:q} "
        "> {log:q} 2>&1"


rule normalize_unmerged_16S_PR2_rank_prefixes:
    input:
        taxonomy=rules.export_unmerged_16S_PR2_taxonomy.output
    output:
        normalized=temp(".snakemake/classification-exports/unmerged-16S/PR2/taxonomy.tsv")
    log:
        "logs/02-denoise-and-export-prok/11-normalize-PR2-rank-prefixes.log"
    # Normalize the Taxon column only. RESCRIPt's edit-taxonomy action accepts
    # a taxonomy Series and discards Confidence; the supplement needs the
    # original confidence values for its unchanged plastid threshold.
    script:
        "../scripts/normalize_pr2_taxonomy.py"


rule convert_unmerged_16S_retained_biom:
    input:
        biom=rules.export_unmerged_16S_retained_table.output
    output:
        temp(".snakemake/classification-exports/unmerged-16S/retained-table.tsv")
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-prok/11-convert-retained-biom.log"
    shell:
        "biom convert -i {input.biom:q} -o {output:q} --to-tsv "
        "> {log:q} 2>&1"


rule convert_unmerged_16S_standard_biom:
    input:
        biom=rules.export_unmerged_16S_standard_table.output
    output:
        temp(".snakemake/classification-exports/unmerged-16S/standard-table.tsv")
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-prok/11-convert-standard-biom.log"
    shell:
        "biom convert -i {input.biom:q} -o {output:q} --to-tsv "
        "> {log:q} 2>&1"


rule export_unmerged_16S_supplement:
    input:
        retained_table_tsv=rules.convert_unmerged_16S_retained_biom.output,
        standard_table_tsv=rules.convert_unmerged_16S_standard_biom.output,
        sequence_fasta=rules.export_unmerged_16S_sequences.output,
        silva_export=rules.export_unmerged_16S_SILVA_taxonomy.output,
        pr2_export=rules.normalize_unmerged_16S_PR2_rank_prefixes.output.normalized,
    output:
        table=UNMERGED_16S_TABLE,
        summary=UNMERGED_16S_SUMMARY,
    params:
        min_confidence=CHLOROPLAST_PR2_MIN_CONFIDENCE,
    script:
        "../scripts/export_unmerged_16S.py"

rule export_DADA2_results:
    input:
        "results/02-proks/03-DADA2d/"
    params:
        studyName=config["studyName"]
    output:
        directory("results/02-proks/04-DADA2d-plaintext-exports/"),
        lateststats="results/02-proks/04-DADA2d-plaintext-exports/" + config["studyName"] + ".16S.latest_stats.tsv",
        latestseqs="results/02-proks/04-DADA2d-plaintext-exports/" + config["studyName"] + ".16S.latest_seqs.fasta"
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-prok/04_export_DADA2_results/DADA2_export.stderrout"
    script:
        "../scripts/P04-export-DADA2-results.sh"

rule classify_ASVs:
    input:
        sequences="results/02-proks/03-DADA2d/representative_sequences.qza",
        database_ready=rules.ensure_silva_classifier.output.marker,
    params:
        classDB=SILVA_CLASSIFIER,
    output:
        directory("results/02-proks/05-classified/"),
        classified="results/02-proks/05-classified/" + config["studyName"] + "_SILVA.classified.qza"
    conda:
        config["qiime2version"]
    threads: 8
    resources:
        mem_mb=64000,
    log:
        "logs/02-denoise-and-export-prok/05-classify-ASVs.log"
    script:
        "../scripts/P05-classify-eASVs.sh"

rule create_sample_metadata_file:
    input:
        manifest="results/02-proks/manifest.tsv",
        samplesdottsv="config/samples.tsv",
        prokstats=rules.export_DADA2_results.output.lateststats,
        eukfracpersample="results/" + config["studyName"] + ".eukfrac-per-sample.tsv"
    output:
        "results/02-proks/sample-metadata.tsv"
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-prok/06-make-sample-metadata-file.log"
    script:
        "../scripts/P06-make-sample-metadata-file.py"

rule make_SILVA_only_prok_barplots:
    input:
        proktable="results/02-proks/03-DADA2d/table.qza",
        proktax=rules.classify_ASVs.output.classified,
        prokmetadata="results/02-proks/sample-metadata.tsv"
    params:
        studyName=config["studyName"]
    output:
        directory("results/02-proks/07-SILVA-only-barplots/")
    conda:
        config["qiime2version"]
    script:
        "../scripts/P07-make-barplot.sh"

rule classify_all_16S_with_PR2:
    input:
        database_ready=rules.ensure_pr2_classifier.output.marker,
        prokseqs=rules.denoise_prok_dada2.output.prokrepseqs,
    params:
        PR2classifier=PR2_CLASSIFIER,
    output:
        classified="results/02-proks/09-subsetting/reclassified/all_16S_ASVs_PR2.classified.qza",
    # Limit classifier copies in memory; this rule can be large even when the
    # standard SILVA classification completed comfortably with more workers.
    threads: 2
    resources:
        mem_mb=110000,
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-prok/09-classify-all-16S-with-PR2.log"
    script:
        "../scripts/P09a-classify-all-16S-with-PR2.sh"


# AI-assisted change — see provenance ledger linked from README.md, commit 92f67ee.
# Purpose: make the SILVA-first, PR2-plastid-resolution order explicit. PR2
# replaces SILVA only for :plas lineages with confidence >= 0.7.
# Human review status: pending.
rule export_SILVA_taxonomy_for_plastid_resolution:
    input:
        rules.classify_ASVs.output.classified
    output:
        temp(".snakemake/classification-exports/canonical-16S/SILVA/taxonomy.tsv")
    params:
        export_dir=".snakemake/classification-exports/canonical-16S/SILVA"
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-prok/09-export-SILVA-taxonomy.log"
    shell:
        "mkdir -p {params.export_dir:q} && "
        "qiime tools export --input-path {input:q} --output-path {params.export_dir:q} "
        "> {log:q} 2>&1"


rule export_PR2_taxonomy_for_plastid_resolution:
    input:
        rules.classify_all_16S_with_PR2.output.classified
    output:
        temp(".snakemake/classification-exports/canonical-16S/PR2-raw/taxonomy.tsv")
    params:
        export_dir=".snakemake/classification-exports/canonical-16S/PR2-raw"
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-prok/09-export-PR2-taxonomy.log"
    shell:
        "mkdir -p {params.export_dir:q} && "
        "qiime tools export --input-path {input:q} --output-path {params.export_dir:q} "
        "> {log:q} 2>&1"


rule normalize_canonical_16S_PR2_rank_prefixes:
    input:
        taxonomy=rules.export_PR2_taxonomy_for_plastid_resolution.output
    output:
        normalized=temp(".snakemake/classification-exports/canonical-16S/PR2/taxonomy.tsv")
    log:
        "logs/02-denoise-and-export-prok/09-normalize-PR2-rank-prefixes.log"
    # Keep classifier confidence and feature IDs intact. Only prefixed rank
    # notation changes, so PR2 lineages cannot masquerade as SILVA taxonomy.
    # Reconciliation still applies the same >=0.7 plastid threshold afterward.
    script:
        "../scripts/normalize_pr2_taxonomy.py"


rule resolve_PR2_plastids:
    input:
        silva_taxonomy=rules.export_SILVA_taxonomy_for_plastid_resolution.output,
        pr2_taxonomy=rules.normalize_canonical_16S_PR2_rank_prefixes.output.normalized,
    output:
        taxonomy_tsv=temp("results/02-proks/09-subsetting/tax-merged/" + config["studyName"] + ".PR2-plastid-resolved-taxonomy.tsv"),
        audit="results/02-proks/09-subsetting/tax-merged/" + config["studyName"] + ".PR2-plastid-routing-audit.tsv",
        summary="results/02-proks/09-subsetting/tax-merged/" + config["studyName"] + ".PR2-plastid-routing-summary.tsv",
    params:
        min_confidence=CHLOROPLAST_PR2_MIN_CONFIDENCE,
    script:
        "../scripts/resolve_pr2_plastids.py"


rule import_resolved_PR2_plastid_taxonomy:
    input:
        rules.resolve_PR2_plastids.output.taxonomy_tsv
    output:
        taxonomy="results/02-proks/09-subsetting/tax-merged/chloroplasts-PR2-reclassified-merged-classification.qza"
    conda:
        config["qiime2version"]
    log:
        "logs/02-denoise-and-export-prok/09-import-resolved-taxonomy.log"
    shell:
        "qiime tools import --type 'FeatureData[Taxonomy]' "
        "--input-format TSVTaxonomyFormat --input-path {input:q} "
        "--output-path {output.taxonomy:q} > {log:q} 2>&1"


rule split_resolved_prok_categories:
    input:
        proktable=rules.denoise_prok_dada2.output.proktable,
        resolved_taxonomy=rules.import_resolved_PR2_plastid_taxonomy.output.taxonomy,
    output:
        includechlorotable="results/02-proks/09-subsetting/split-tables/include_o__Chloroplast_filtered_table.qza",
        excludechlorotable="results/02-proks/09-subsetting/split-tables/exclude_o__Chloroplast_filtered_table.qza",
        onlymitotable="results/02-proks/09-subsetting/split-tables/include_f__Mitochondria_filtered_table.qza",
        onlyalgaetable="results/02-proks/09-subsetting/split-tables/include_p__Cyanobacteria_NOTE_includes_chloroplasts_filtered_table.qza",
        onlycyanotable="results/02-proks/09-subsetting/split-tables/include_p__Cyanobacteria_exclude_o__Chloroplast_filtered_table.qza",
        nomitotable="results/02-proks/09-subsetting/split-tables/exclude_f__Mitochondria_filtered_table.qza",
        nomitonochlorotable="results/02-proks/09-subsetting/split-tables/exclude_o__Chloroplast_exclude_f__Mitochondria_filtered_table.qza",
        nomitonochloronocyanotable="results/02-proks/09-subsetting/split-tables/exclude_p__Cyanobacteria_exclude_f__Mitochondria_NOTE_excludes_chloroplasts_filtered_table.qza",
        onlyarchaeatable="results/02-proks/09-subsetting/split-tables/include_d__Archaea_filtered_table.qza",
        noarchaeatable="results/02-proks/09-subsetting/split-tables/exclude_d__Archaea_filtered_table.qza" 
    conda:
        config["qiime2version"]
    script:
        "../scripts/P09b-split-resolved-prok-categories.sh"

rule export_tax_convert_biom:
    input:
        mergedtax=rules.import_resolved_PR2_plastid_taxonomy.output.taxonomy,
        all16Stable="results/02-proks/03-DADA2d/table.qza",
        noarch="results/02-proks/09-subsetting/split-tables/exclude_d__Archaea_filtered_table.qza",
        nomito="results/02-proks/09-subsetting/split-tables/exclude_f__Mitochondria_filtered_table.qza",
        nochloronomito="results/02-proks/09-subsetting/split-tables/exclude_o__Chloroplast_exclude_f__Mitochondria_filtered_table.qza",
        nochloro="results/02-proks/09-subsetting/split-tables/exclude_o__Chloroplast_filtered_table.qza",
        nochloronocyanonomito="results/02-proks/09-subsetting/split-tables/exclude_p__Cyanobacteria_exclude_f__Mitochondria_NOTE_excludes_chloroplasts_filtered_table.qza",
        onlyarch="results/02-proks/09-subsetting/split-tables/include_d__Archaea_filtered_table.qza",
        onlymito="results/02-proks/09-subsetting/split-tables/include_f__Mitochondria_filtered_table.qza",
        onlychloro="results/02-proks/09-subsetting/split-tables/include_o__Chloroplast_filtered_table.qza",
        onlycyano="results/02-proks/09-subsetting/split-tables/include_p__Cyanobacteria_exclude_o__Chloroplast_filtered_table.qza",
        onlyalgae="results/02-proks/09-subsetting/split-tables/include_p__Cyanobacteria_NOTE_includes_chloroplasts_filtered_table.qza",
    output:
        genus_mergedtaxtsv=temp("results/02-proks/10-exports/" + config["studyName"] + ".taxonomy-through-genus.tsv"),
        all16Stable_biom=temp("results/02-proks/10-exports/" + config["studyName"] + ".all-16S-seqs.biom"),
        noarch_biom=temp("results/02-proks/10-exports/" + config["studyName"] + ".exclude_d__Archaea_filtered_table.biom"),
        nomito_biom=temp("results/02-proks/10-exports/" + config["studyName"] + ".exclude_f__Mitochondria_filtered_table.biom"),
        nochloronomito_biom=temp("results/02-proks/10-exports/" + config["studyName"] + ".exclude_o__Chloroplast_exclude_f__Mitochondria_filtered_table.biom"),
        nochloro_biom=temp("results/02-proks/10-exports/" + config["studyName"] + ".exclude_o__Chloroplast_filtered_table.biom"),
        nochloronocyanonomito_biom=temp("results/02-proks/10-exports/" + config["studyName"] + ".exclude_p__Cyanobacteria_exclude_f__Mitochondria_NOTE_excludes_chloroplasts_filtered_table.biom"),
        onlyarch_biom=temp("results/02-proks/10-exports/" + config["studyName"] + ".include_d__Archaea_filtered_table.biom"),
        onlymito_biom=temp("results/02-proks/10-exports/" + config["studyName"] + ".include_f__Mitochondria_filtered_table.biom"),
        onlychloro_biom=temp("results/02-proks/10-exports/" + config["studyName"] + ".include_o__Chloroplast_filtered_table.biom"),
        onlycyano_biom=temp("results/02-proks/10-exports/" + config["studyName"] + ".include_p__Cyanobacteria_exclude_o__Chloroplast_filtered_table.biom"),
        onlyalgae_biom=temp("results/02-proks/10-exports/" + config["studyName"] + ".include_p__Cyanobacteria_NOTE_includes_chloroplasts_filtered_table.qza.biom"),
    conda:
        config["qiime2version"]
    script:
        "../scripts/P10a-generate-biom-tables.sh"

rule assign_SILVA_species:
    input:
        taxonomy=rules.export_tax_convert_biom.output.genus_mergedtaxtsv,
        sequences=rules.export_DADA2_results.output.latestseqs,
        # The exact-match reference belongs to release 144 only. Older
        # classifiers pass their taxonomy through without loading this file.
        species_reference=([SILVA_SPECIES_REFERENCE] if RUN_SILVA_SPECIES_MATCHING else []),
    output:
        taxonomy="results/02-proks/10-exports/" + config["studyName"] + ".taxonomy.tsv",
        summary="results/02-proks/10-exports/" + config["studyName"] + ".SILVA-" + SILVA_VERSION + ".species-assignment-summary.tsv",
    params:
        chunk_size=500,
        enabled=RUN_SILVA_SPECIES_MATCHING,
        silva_version=SILVA_VERSION,
    conda:
        config["qiime2version"]
    threads: 1
    resources:
        # Reserve the node for this step. assignSpecies is also batched below,
        # so its actual peak memory no longer grows with the full ASV count.
        mem_mb=110000,
        runtime=360,
    log:
        "logs/02-denoise-and-export-prok/10-assign-SILVA-species.log"
    script:
        "../scripts/assign_silva_species.R"

rule add_tax_to_biom:
    input:
        mergedtaxtsv=rules.assign_SILVA_species.output.taxonomy,
        all16Stable_biom="results/02-proks/10-exports/" + config["studyName"] + ".all-16S-seqs.biom",
        noarch_biom="results/02-proks/10-exports/" + config["studyName"] + ".exclude_d__Archaea_filtered_table.biom",
        nomito_biom="results/02-proks/10-exports/" + config["studyName"] + ".exclude_f__Mitochondria_filtered_table.biom",
        nochloronomito_biom="results/02-proks/10-exports/" + config["studyName"] + ".exclude_o__Chloroplast_exclude_f__Mitochondria_filtered_table.biom",
        nochloro_biom="results/02-proks/10-exports/" + config["studyName"] + ".exclude_o__Chloroplast_filtered_table.biom",
        nochloronocyanonomito_biom="results/02-proks/10-exports/" + config["studyName"] + ".exclude_p__Cyanobacteria_exclude_f__Mitochondria_NOTE_excludes_chloroplasts_filtered_table.biom",
        onlyarch_biom="results/02-proks/10-exports/" + config["studyName"] + ".include_d__Archaea_filtered_table.biom",
        onlymito_biom="results/02-proks/10-exports/" + config["studyName"] + ".include_f__Mitochondria_filtered_table.biom",
        onlychloro_biom="results/02-proks/10-exports/" + config["studyName"] + ".include_o__Chloroplast_filtered_table.biom",
        onlycyano_biom="results/02-proks/10-exports/" + config["studyName"] + ".include_p__Cyanobacteria_exclude_o__Chloroplast_filtered_table.biom",
        onlyalgae_biom="results/02-proks/10-exports/" + config["studyName"] + ".include_p__Cyanobacteria_NOTE_includes_chloroplasts_filtered_table.qza.biom",
    output:
        all16Stable_biomtax=temp("results/02-proks/10-exports/" + config["studyName"] + ".all-16S-seqs.with-tax.biom"),
        noarch_biomtax=temp("results/02-proks/10-exports/" + config["studyName"] + ".exclude_d__Archaea_filtered_table.with-tax.biom"),
        nomito_biomtax=temp("results/02-proks/10-exports/" + config["studyName"] + ".exclude_f__Mitochondria_filtered_table.with-tax.biom"),
        nochloronomito_biomtax=temp("results/02-proks/10-exports/" + config["studyName"] + ".exclude_o__Chloroplast_exclude_f__Mitochondria_filtered_table.with-tax.biom"),
        nochloro_biomtax=temp("results/02-proks/10-exports/" + config["studyName"] + ".exclude_o__Chloroplast_filtered_table.with-tax.biom"),
        nochloronocyanonomito_biomtax=temp("results/02-proks/10-exports/" + config["studyName"] + ".exclude_p__Cyanobacteria_exclude_f__Mitochondria_NOTE_excludes_chloroplasts_filtered_table.with-tax.biom"),
        onlyarch_biomtax=temp("results/02-proks/10-exports/" + config["studyName"] + ".include_d__Archaea_filtered_table.with-tax.biom"),
        onlymito_biomtax=temp("results/02-proks/10-exports/" + config["studyName"] + ".include_f__Mitochondria_filtered_table.with-tax.biom"),
        onlychloro_biomtax=temp("results/02-proks/10-exports/" + config["studyName"] + ".include_o__Chloroplast_filtered_table.with-tax.biom"),
        onlycyano_biomtax=temp("results/02-proks/10-exports/" + config["studyName"] + ".include_p__Cyanobacteria_exclude_o__Chloroplast_filtered_table.with-tax.biom"),
        onlyalgae_biomtax=temp("results/02-proks/10-exports/" + config["studyName"] + ".include_p__Cyanobacteria_NOTE_includes_chloroplasts_filtered_table.with-tax.biom")
    conda:
        config["qiime2version"]
    script:
        "../scripts/P10b-add-tax-to-biom.sh"

rule export_biom_tsv:
    input:
        all16Stable_biomtax="results/02-proks/10-exports/" + config["studyName"] + ".all-16S-seqs.with-tax.biom",
        noarch_biomtax="results/02-proks/10-exports/" + config["studyName"] + ".exclude_d__Archaea_filtered_table.with-tax.biom",
        nomito_biomtax="results/02-proks/10-exports/" + config["studyName"] + ".exclude_f__Mitochondria_filtered_table.with-tax.biom",
        nochloronomito_biomtax="results/02-proks/10-exports/" + config["studyName"] + ".exclude_o__Chloroplast_exclude_f__Mitochondria_filtered_table.with-tax.biom",
        nochloro_biomtax="results/02-proks/10-exports/" + config["studyName"] + ".exclude_o__Chloroplast_filtered_table.with-tax.biom",
        nochloronocyanonomito_biomtax="results/02-proks/10-exports/" + config["studyName"] + ".exclude_p__Cyanobacteria_exclude_f__Mitochondria_NOTE_excludes_chloroplasts_filtered_table.with-tax.biom",
        onlyarch_biomtax="results/02-proks/10-exports/" + config["studyName"] + ".include_d__Archaea_filtered_table.with-tax.biom",
        onlymito_biomtax="results/02-proks/10-exports/" + config["studyName"] + ".include_f__Mitochondria_filtered_table.with-tax.biom",
        onlychloro_biomtax="results/02-proks/10-exports/" + config["studyName"] + ".include_o__Chloroplast_filtered_table.with-tax.biom",
        onlycyano_biomtax="results/02-proks/10-exports/" + config["studyName"] + ".include_p__Cyanobacteria_exclude_o__Chloroplast_filtered_table.with-tax.biom",
        onlyalgae_biomtax="results/02-proks/10-exports/" + config["studyName"] + ".include_p__Cyanobacteria_NOTE_includes_chloroplasts_filtered_table.with-tax.biom",
    output:
        all16Stable_biomtaxtsv="results/02-proks/10-exports/" + config["studyName"] + ".all-16S-seqs.with-tax.tsv",
        noarch_biomtaxtsv="results/02-proks/10-exports/" + config["studyName"] + ".exclude_d__Archaea_filtered_table.with-tax.tsv",
        nomito_biomtaxtsv="results/02-proks/10-exports/" + config["studyName"] + ".exclude_f__Mitochondria_filtered_table.with-tax.tsv",
        nochloronomito_biomtaxtsv="results/02-proks/10-exports/" + config["studyName"] + ".exclude_o__Chloroplast_exclude_f__Mitochondria_filtered_table.with-tax.tsv",
        nochloro_biomtaxtsv="results/02-proks/10-exports/" + config["studyName"] + ".exclude_o__Chloroplast_filtered_table.with-tax.tsv",
        nochloronocyanonomito_biomtaxtsv="results/02-proks/10-exports/" + config["studyName"] + ".exclude_p__Cyanobacteria_exclude_f__Mitochondria_NOTE_excludes_chloroplasts_filtered_table.with-tax.tsv",
        onlyarch_biomtaxtsv="results/02-proks/10-exports/" + config["studyName"] + ".include_d__Archaea_filtered_table.with-tax.tsv",
        onlymito_biomtaxtsv="results/02-proks/10-exports/" + config["studyName"] + ".include_f__Mitochondria_filtered_table.with-tax.tsv",
        onlychloro_biomtaxtsv="results/02-proks/10-exports/" + config["studyName"] + ".include_o__Chloroplast_filtered_table.with-tax.tsv",
        onlycyano_biomtaxtsv="results/02-proks/10-exports/" + config["studyName"] + ".include_p__Cyanobacteria_exclude_o__Chloroplast_filtered_table.with-tax.tsv",
        onlyalgae_biomtaxtsv="results/02-proks/10-exports/" + config["studyName"] + ".include_p__Cyanobacteria_NOTE_includes_chloroplasts_filtered_table.with-tax.tsv",
    conda:
        config["qiime2version"]
    script:
        "../scripts/P10c-export-tax-tsvs.sh"

rule proportal_classify:
    input:
        proktable="results/02-proks/10-exports/" + config["studyName"] + ".all-16S-seqs.with-tax.tsv",
        dnaseqs="results/02-proks/04-DADA2d-plaintext-exports/" + config["studyName"] + ".16S.latest_seqs.fasta"
    output:
        taxonomy="results/02-proks/10-exports/" + config["studyName"] + ".Synechococcales.proportal-classified.tsv",
    conda:
        "../envs/proportal.yml"
    script:
        "../scripts/add-proportal-taxonomy.sh"
