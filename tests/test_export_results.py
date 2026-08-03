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
            source.write_text("sample\tvalue\nS1\t1\n", encoding="utf-8")
            destination = root / "Results-Export" / "study.long_data.tsv"
            readme = root / "Results-Export" / "README.txt"

            EXPORT.export_results([source], [destination], readme, "study")

            self.assertEqual(destination.read_text(encoding="utf-8"), source.read_text(encoding="utf-8"))
            self.assertIn("study.long_data.tsv", readme.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
