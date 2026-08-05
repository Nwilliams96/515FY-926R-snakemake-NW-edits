import importlib.util
import tempfile
import unittest
import zipfile
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "workflow" / "scripts" / "generate_pipeline_report.py"
SPEC = importlib.util.spec_from_file_location("generate_pipeline_report", SCRIPT)
REPORT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REPORT)


class GeneratePipelineReportTests(unittest.TestCase):
    def test_compiles_after_snakemake_generated_preamble(self):
        source = "snakemake = None\n" + SCRIPT.read_text(encoding="utf-8")
        compile(source, str(SCRIPT), "exec")

    def test_renders_summary_and_optional_figure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "samples.tsv").write_text("sample\tcondition\nS_1\ttest\n", encoding="utf-8")
            (root / "split.tsv").write_text(
                "sample\tprok_seqs_split\teuk_seqs_split\ttotal_seqs_split\teukfrac_split\n"
                "S_1\t700\t300\t1000\t0.3\n",
                encoding="utf-8",
            )
            stats = (
                "sample-id\tinput\tfiltered\tnon-chimeric\n"
                "S-1\t{input}\t{filtered}\t{final}\n"
            )
            (root / "16.tsv").write_text(stats.format(input=700, filtered=650, final=600), encoding="utf-8")
            (root / "18.tsv").write_text(stats.format(input=300, filtered=250, final=200), encoding="utf-8")
            (root / "long.tsv").write_text(
                "SampleID\tDomain\tPhylum\tDivision\tASV_hash\tCorrected_Sequence_Counts\n"
                "S-1\tBacteria\tProteobacteria\t\tA1\t600\n"
                "S-1\tEukaryota\t\tDinoflagellata\tA2\t200\n",
                encoding="utf-8",
            )
            (root / "trim.log").write_text(
                "Total read pairs processed: 1,200\nPairs written (passing filters): 1,000\n",
                encoding="utf-8",
            )
            quality_dir = root / "quality"
            quality_dir.mkdir()
            with zipfile.ZipFile(quality_dir / "visualization.qzv", "w") as archive:
                archive.writestr(
                    "example/data/forward-seven-number-summaries.tsv",
                    "position\t50%\n1\t35\n2\t34\n",
                )
            # A tiny valid PNG is sufficient to verify base64 embedding.
            png = root / "figure.png"
            png.write_bytes(
                b"\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x01\x00\x00\x00\x01"
                b"\x08\x06\x00\x00\x00\x1f\x15\xc4\x89"
            )
            output = root / "report.html"
            REPORT.render_report(
                {"studyName": "test-run", "trunclens": {"truncR1": 220, "truncR2": 180}},
                {
                    "samples": root / "samples.tsv",
                    "split_summary": root / "split.tsv",
                    "stats16s": root / "16.tsv",
                    "stats18s": root / "18.tsv",
                    "long_data": root / "long.tsv",
                    "quality_directories": [quality_dir],
                    "trimming_logs": [root / "trim.log"],
                    "internal_standard_figures": [png],
                },
                output,
            )
            rendered = output.read_text(encoding="utf-8")
            self.assertIn("test-run", rendered)
            self.assertIn("Raw paired reads", rendered)
            self.assertIn("S-1", rendered)
            self.assertIn("Median Phred score", rendered)
            self.assertIn("Proteobacteria", rendered)
            self.assertIn("data:image/png;base64,", rendered)
            self.assertIn("trunclens.truncR1", rendered)


if __name__ == "__main__":
    unittest.main()
