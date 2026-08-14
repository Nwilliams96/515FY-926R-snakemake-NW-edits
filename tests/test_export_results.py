import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "workflow" / "scripts" / "export_results.py"
SPEC = importlib.util.spec_from_file_location("export_results", SCRIPT)
EXPORT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(EXPORT)


class ExportResultsTests(unittest.TestCase):
    def test_copies_files_and_writes_readme(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source.tsv"
            source.write_text(
                "SampleID\tPhylum\tOrder\tCorrected_Sequence_Counts\n"
                "S1\tP1\tO1\t3\nS1\tP1\tO2\t1\n",
                encoding="utf-8",
            )
            destination = root / "Results-Export" / "study.long_data.tsv"
            readme = root / "Results-Export" / "README.txt"
            summaries = [
                root / "Results-Export" / "study.phylum_summary.tsv",
                root / "Results-Export" / "study.order_summary.tsv",
            ]

            EXPORT.export_results(
                [source], [destination], readme, "study", summaries
            )

            self.assertEqual(destination.read_text(encoding="utf-8"), source.read_text(encoding="utf-8"))
            self.assertIn("study.long_data.tsv", readme.read_text(encoding="utf-8"))
            self.assertIn("S1\tP1\t4\t1", summaries[0].read_text(encoding="utf-8"))
            self.assertIn("study.order_summary.tsv", readme.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
