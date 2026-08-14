import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REFORMAT_SCRIPT = ROOT / "workflow/scripts/reformat_pr2_reference.py"
SPEC = importlib.util.spec_from_file_location("reformat_pr2_reference", REFORMAT_SCRIPT)
REFORMAT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REFORMAT)


class DatabaseGenerationTests(unittest.TestCase):
    def test_database_rules_are_self_contained_and_use_supported_parallelism(self):
        classification = (ROOT / "workflow/rules/prepare_classification_dbs.smk").read_text(
            encoding="utf-8"
        )
        downloads = (ROOT / "workflow/rules/download_databases.smk").read_text(
            encoding="utf-8"
        )
        bbsplit = (ROOT / "workflow/rules/prepare_bbsplit_db.smk").read_text(
            encoding="utf-8"
        )
        self.assertNotIn("tax-classifier-construction", classification)
        self.assertEqual(classification.count("--p-n-jobs {threads}"), 2)
        self.assertGreaterEqual(classification.count("threads: 8"), 2)
        self.assertEqual(classification.count("fit-classifier-naive-bayes"), 2)
        self.assertIn("threads: 8", bbsplit)
        self.assertIn("bbsplit.sh build=1 threads={threads}", bbsplit)
        self.assertIn("rule initialize_database_directories:", downloads)
        self.assertIn('DATABASE_PREFIX + "classification/SILVA"', downloads)

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
        self.assertIn('f"{standard_id}_recovery_ratio"', common)
        self.assertIn('f"mean_{first}_and_{second}_recovery_ratio"', common)
        self.assertIn("corrected=ISD_CORRECTED_LONG_TABLE", correction_rule)

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
