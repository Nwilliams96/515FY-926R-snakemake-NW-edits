import importlib.util
import hashlib
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REFORMAT_SCRIPT = ROOT / "workflow/scripts/reformat_pr2_reference.py"
SPEC = importlib.util.spec_from_file_location("reformat_pr2_reference", REFORMAT_SCRIPT)
REFORMAT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REFORMAT)

DOWNLOAD_SCRIPT = ROOT / "workflow/scripts/download_verified_reference.py"
DOWNLOAD_SPEC = importlib.util.spec_from_file_location(
    "download_verified_reference", DOWNLOAD_SCRIPT
)
DOWNLOAD = importlib.util.module_from_spec(DOWNLOAD_SPEC)
DOWNLOAD_SPEC.loader.exec_module(DOWNLOAD)


class DatabaseGenerationTests(unittest.TestCase):
    def test_database_rules_are_self_contained_and_use_supported_parallelism(self):
        classification = (ROOT / "workflow/rules/prepare_classification_dbs.smk").read_text(
            encoding="utf-8"
        )
        common = (ROOT / "workflow/rules/common.smk").read_text(encoding="utf-8")
        downloads = (ROOT / "workflow/rules/download_databases.smk").read_text(
            encoding="utf-8"
        )
        bbsplit = (ROOT / "workflow/rules/prepare_bbsplit_db.smk").read_text(
            encoding="utf-8"
        )
        pr2_builder = (ROOT / "workflow/scripts/ensure_pr2_classifier.py").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("tax-classifier-construction", classification)
        self.assertNotIn("rule import_pr2_fasta:", classification)
        self.assertNotIn("rule derep_seqs_taxonomy:", classification)
        self.assertIn("rule ensure_pr2_classifier:", classification)
        self.assertGreaterEqual(classification.count("threads: 8"), 1)
        self.assertEqual(pr2_builder.count('"fit-classifier-naive-bayes"'), 1)
        self.assertIn('"--o-read-extraction-stats", extraction_stats_qza', pr2_builder)
        self.assertIn("fcntl.flock", pr2_builder)
        self.assertIn('prefix=".pr2-build-"', pr2_builder)
        self.assertIn("os.replace(staged_classifier, classifier)", pr2_builder)
        self.assertIn('["qiime", "tools", "validate", str(path)]', pr2_builder)
        self.assertIn("rule download_SILVA_144_classifier:", classification)
        self.assertIn("rule download_SILVA_144_species_reference:", classification)
        self.assertNotIn("get-silva-data", classification)
        self.assertIn('SILVA_VERSION = "144"', common)
        self.assertIn("SILVA_144_SSURef_NR99_uniform_classifier_V4V5-515f-926r.qza", common)
        self.assertIn("f7757b01eb82e0ac78bf06427e410095", common)
        self.assertIn("silva_v144_assignSpecies.fa.gz", common)
        self.assertIn("444de7c0cce0b66addda7a3f8b38e012", common)
        self.assertIn("rachis-qiime2-linux-64-2026.7.yml", common)
        self.assertIn("_dereplicated_final_classifier_qiime2-2026.7.qza", common)
        self.assertIn("marker=PR2_READY_MARKER", classification)
        self.assertIn("threads: 8", bbsplit)
        self.assertIn("bbsplit.sh build=1 threads={threads}", bbsplit)
        self.assertIn("rule initialize_database_directories:", downloads)
        self.assertIn('DATABASE_PREFIX + "classification/SILVA"', downloads)
        self.assertNotIn("rule download_pr2:", downloads)
        self.assertIn(
            'source_checksum="0c8728abcbb2126eed2c7e587f820cbce39c138cdfdb51239bbf18621498462d"',
            classification,
        )
        self.assertIn(
            'script:\n        "../scripts/download_verified_reference.py"',
            classification,
        )
        self.assertIn("USE_PREEXISTING_BBSPLIT_DATABASE", common)
        self.assertIn("USE_PREEXISTING_SILVA_DATABASE", common)
        self.assertIn("USE_PREEXISTING_PR2_DATABASE", common)
        snakefile = (ROOT / "workflow/Snakefile").read_text(encoding="utf-8")
        self.assertIn("if not USE_PREEXISTING_BBSPLIT_DATABASE:", snakefile)

        species_script = (ROOT / "workflow/scripts/assign_silva_species.R").read_text(
            encoding="utf-8"
        )
        self.assertIn("allowMultiple = FALSE", species_script)
        self.assertIn("tryRC = TRUE", species_script)
        self.assertIn("Genus_Agrees", species_script)
        self.assertIn("assign_species_in_chunks", species_script)
        self.assertIn("stage_reference_locally", species_script)
        self.assertIn("Using node-local SILVA species reference", species_script)
        self.assertIn("gc(full = TRUE)", species_script)
        species_rule = (
            ROOT / "workflow/rules/02-denoise-and-export-prok.smk"
        ).read_text(encoding="utf-8")
        species_rule = species_rule.split("rule assign_SILVA_144_species:", 1)[1].split(
            "\nrule ", 1
        )[0]
        self.assertIn("chunk_size=500", species_rule)
        self.assertIn("mem_mb=110000", species_rule)
        qiime_env = (
            ROOT / "workflow/envs/rachis-qiime2-linux-64-2026.7.yml"
        ).read_text(encoding="utf-8")
        self.assertIn("bioconductor-dada2=", qiime_env)
        self.assertIn("q2-dada2=2026.7.0", qiime_env)
        self.assertIn("q2-feature-classifier=2026.7.0", qiime_env)
        self.assertIn("rescript=2026.7.0", qiime_env)
        self.assertIn("scikit-learn=1.7.1", qiime_env)
        self.assertNotIn("\n- deblur=", qiime_env)
        self.assertNotIn("\n- q2-deblur=", qiime_env)
        self.assertNotIn("\n- sortmerna=", qiime_env)

        self.assertNotIn(
            "input:\n        temp(",
            downloads,
        )

    def test_reference_download_is_atomic_and_checksum_verified(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source.qza"
            destination = root / "classification" / "classifier.qza"
            source.write_bytes(b"signed SILVA classifier fixture")
            checksum = hashlib.md5(source.read_bytes()).hexdigest()

            DOWNLOAD.download_verified(source.as_uri(), destination, checksum)
            self.assertEqual(destination.read_bytes(), source.read_bytes())
            self.assertFalse(destination.with_name(destination.name + ".part").exists())

            destination.unlink()
            with self.assertRaises(ValueError):
                DOWNLOAD.download_verified(source.as_uri(), destination, "0" * 32)
            self.assertFalse(destination.exists())
            self.assertFalse(destination.with_name(destination.name + ".part").exists())

            sha_destination = root / "absolute-style" / "reference.fasta.gz"
            sha_checksum = hashlib.sha256(source.read_bytes()).hexdigest()
            DOWNLOAD.download_verified(
                source.as_uri(), sha_destination, sha_checksum, "sha256"
            )
            self.assertEqual(sha_destination.read_bytes(), source.read_bytes())

    def test_shared_reference_download_uses_lock_and_unique_staging(self):
        source = DOWNLOAD_SCRIPT.read_text(encoding="utf-8")
        self.assertIn("fcntl.flock", source)
        self.assertIn("destination_lock(destination)", source)
        self.assertIn("tempfile.NamedTemporaryFile", source)
        self.assertIn("os.replace(partial, destination)", source)
        self.assertNotIn('destination.with_name(destination.name + ".part")', source)

    def test_long_data_parser_is_rank_named_and_silva_144_safe(self):
        source = (ROOT / "workflow/scripts/long-data-preparation.R").read_text(
            encoding="utf-8"
        )
        self.assertIn('Kingdom = extract_prefixed_rank(Taxonomy, "k__")', source)
        self.assertIn("parse_prefixed_taxonomy", source)
        self.assertNotIn('separate(Taxonomy, c("Domain","Phylum"', source)
        self.assertIn('"Domain", "Kingdom", "Supergroup"', source)
        self.assertIn(
            '"Taxonomy", "Domain", "Kingdom", "Supergroup"', source
        )

        subset_script = (
            ROOT / "workflow/scripts/P09b-split-resolved-prok-categories.sh"
        ).read_text(encoding="utf-8")
        self.assertIn("p__Cyanobacteria,p__Cyanobacteriota", subset_script)
        self.assertIn(":plas,o__Chloroplast", subset_script)

        prok_rules = (
            ROOT / "workflow/rules/02-denoise-and-export-prok.smk"
        ).read_text(encoding="utf-8")
        self.assertIn("rule classify_all_16S_with_PR2:", prok_rules)
        self.assertIn("rule resolve_PR2_plastids:", prok_rules)
        self.assertIn("rule split_resolved_prok_categories:", prok_rules)
        self.assertIn("all_16S_ASVs_PR2.classified.qza", prok_rules)
        self.assertIn("resolved_taxonomy=rules.resolve_PR2_plastids.output.taxonomy", prok_rules)
        self.assertNotIn("rule splitchloroplasts:", prok_rules)

        self.assertIn('str_detect(Taxonomy, regex(":plas|', source)
        self.assertIn('plastid_16S_rRNA == "yes" ~ "Chloroplast_16S"', source)

    def test_new_amplicon_concentration_name_and_isd_ids_reach_outputs(self):
        common = (ROOT / "workflow/rules/common.smk").read_text(encoding="utf-8")
        merge_rule = (ROOT / "workflow/rules/03-merge-16S-18S.smk").read_text(
            encoding="utf-8"
        )
        correction_rule = (
            ROOT / "workflow/rules/05-internal-standard-correction.smk"
        ).read_text(encoding="utf-8")
        self.assertIn("prok_and_euk_SSU_amplicon_molarities.tsv", common)
        self.assertIn("LEGACY_AMPLICON_CONCENTRATIONS_FILE", common)
        self.assertIn("amplicon_concentrations=AMPLICON_CONCENTRATIONS_FILE", merge_rule)
        self.assertIn("correction_factors=", merge_rule)
        self.assertIn('f"{standard_id}_recovery_ratio"', common)
        self.assertIn('f"mean_{first}_and_{second}_recovery_ratio"', common)
        self.assertIn("corrected=ISD_CORRECTED_LONG_TABLE", correction_rule)
        self.assertIn("filtered=ISD_FILTERED_LONG_TABLE", correction_rule)
        self.assertIn("RESULTS_LONG_DATA = (", common)
        self.assertIn("ISD_FILTERED_LONG_TABLE", common)

    def test_reformats_pr2_fasta_without_appending_stale_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "pr2.fasta"
            clean = root / "clean.fasta"
            taxonomy = root / "taxonomy.tsv"
            source.write_text(">taxon one\nACGT\n>taxon two\nTGCA\n", encoding="utf-8")
            clean.write_text("stale\n", encoding="utf-8")
            taxonomy.write_text("stale\n", encoding="utf-8")

            REFORMAT.reformat_pr2_reference(source, clean, taxonomy)

            self.assertEqual(clean.read_text(encoding="utf-8"), ">feature_1\nACGT\n>feature_2\nTGCA\n")
            self.assertEqual(
                taxonomy.read_text(encoding="utf-8"),
                "feature_1\ttaxon one\nfeature_2\ttaxon two\n",
            )


if __name__ == "__main__":
    unittest.main()
