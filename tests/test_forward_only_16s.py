import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "workflow/scripts/estimate_forward_only_insert_lengths.py"
SPEC = importlib.util.spec_from_file_location("forward_lengths", SCRIPT)
LENGTHS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(LENGTHS)


class ForwardOnly16STests(unittest.TestCase):
    def test_iupac_primer_matching(self):
        self.assertEqual(
            LENGTHS.mismatch_count(
                "GTGTCAGCAGCCGCGGTAA", "GTGYCAGCMGCCGCGGTAA"
            ),
            0,
        )
        self.assertEqual(
            LENGTHS.mismatch_count(
                "GTGACAGCAGCCGCGGTAA", "GTGYCAGCMGCCGCGGTAA"
            ),
            1,
        )

    def test_reference_geometry_flags_long_insert(self):
        forward_primer = "GTGYCAGCMGCCGCGGTAA"
        reverse_primer = "CCGYCAATTYMTTTRAGTTT"
        reverse_binding = LENGTHS.reverse_complement(reverse_primer)
        concrete_forward = "GTGTCAGCAGCCGCGGTAA"
        insert_length = 390
        reference = concrete_forward + ("A" * insert_length) + reverse_binding
        hit = {
            "subject_start": len(concrete_forward) + 1,
            "subject_end": len(concrete_forward) + 220,
            "query_start": 1,
            "query_length": 220,
            "alignment_length": 220,
            "identity": 100.0,
            "evalue": 0.0,
            "bitscore": 400.0,
            "subject": "synthetic",
        }
        result = LENGTHS.evaluate_hit(
            hit,
            ("synthetic reference", reference),
            forward_primer,
            reverse_binding,
            220,
            180,
            12,
        )
        self.assertIsNotNone(result)
        self.assertEqual(
            result["Estimated_post_primer_insert_length"], insert_length
        )
        self.assertEqual(
            result["Predicted_overlap_at_configured_truncation"], 10
        )
        self.assertEqual(
            result["Length_mergeable_at_configured_min_overlap"], "false"
        )

    def test_optional_target_is_wired_without_changing_default_outputs(self):
        snakefile = (ROOT / "workflow/Snakefile").read_text(encoding="utf-8")
        rules = (
            ROOT / "workflow/rules/08-forward-only-16S-test.smk"
        ).read_text(encoding="utf-8")
        denoise = (
            ROOT / "workflow/scripts/denoise_forward_only_16S.sh"
        ).read_text(encoding="utf-8")
        self.assertIn(
            'include: "rules/08-forward-only-16S-test.smk"', snakefile
        )
        self.assertIn("rule forward_only_16S_test:", rules)
        self.assertIn("default_target: False", rules)
        self.assertIn("qiime dada2 denoise-single", denoise)
        self.assertIn('results/01-split/{sample}.prok.R1.fastq.gz', rules)


if __name__ == "__main__":
    unittest.main()
