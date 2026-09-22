import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))
SCRIPT = SCRIPTS / "export_unmerged_16S.py"
SPEC = importlib.util.spec_from_file_location("export_unmerged_16S", SCRIPT)
RESCUE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RESCUE)


class Unmerged16SRescueTests(unittest.TestCase):
    def test_biom_reader_skips_header_comment(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "table.tsv"
            path.write_text(
                "# Constructed from biom file\n#OTU ID\tB\tA\nmerged\t2\t1\nlinked\t4\t3\n",
                encoding="utf-8",
            )
            self.assertEqual(
                RESCUE.read_biom_tsv(path),
                {"merged": {"B": 2.0, "A": 1.0}, "linked": {"B": 4.0, "A": 3.0}},
            )

    def test_builds_only_features_absent_from_standard_run(self):
        retained = {
            "merged": {"B": 2.0, "A": 1.0},
            "linked": {"B": 4.0, "A": 3.0},
        }
        standard = {"merged": {"B": 2.0, "A": 1.0}}
        silva = [
            {
                "feature_id": "linked",
                "taxonomy": "d__Bacteria; p__Cyanobacteriota",
                "confidence_text": "0.8",
                "confidence": 0.8,
            }
        ]
        pr2 = [
            {
                "feature_id": "linked",
                "taxonomy": "Eukaryota:plas; Chloroplastida:plas",
                "confidence_text": "0.9",
                "confidence": 0.9,
            }
        ]
        rows, summary, samples = RESCUE.build_supplement(
            retained, standard, {"linked": "ACGT...TGCA"}, silva, pr2, 0.7
        )
        self.assertEqual(samples, ["A", "B"])
        self.assertEqual([row["Feature_ID"] for row in rows], ["linked"])
        self.assertEqual(rows[0]["A"], "3")
        self.assertEqual(rows[0]["PR2_override_accepted"], "true")
        self.assertEqual(summary["supplemental_unmerged_16S_features"], 1)
        self.assertEqual(summary["supplemental_unmerged_16S_reads"], "7")
        self.assertEqual(summary["supplemental_chloroplast_reads"], "7")


if __name__ == "__main__":
    unittest.main()
