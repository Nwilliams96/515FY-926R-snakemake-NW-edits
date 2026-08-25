import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "workflow" / "scripts" / "generate_pipeline_report.py"
SPEC = importlib.util.spec_from_file_location("generate_pipeline_report", SCRIPT)
REPORT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REPORT)


class GeneratePipelineReportTests(unittest.TestCase):
    def test_compiles_after_snakemake_generated_preamble(self):
        source = "snakemake = None\n" + SCRIPT.read_text(encoding="utf-8")
        compile(source, str(SCRIPT), "exec")

    def test_non_finite_values_are_treated_as_missing(self):
        for value in ("NaN", "nan", "Inf", "-Inf", float("nan"), float("inf")):
            with self.subTest(value=value):
                self.assertIsNone(REPORT.number(value))
                self.assertEqual(REPORT.fmt_count(value), "—")
                self.assertEqual(REPORT.fmt_percent(value), "—")

        chart = REPORT.svg_log_series(
            ["S-1", "S-2"],
            {
                "S-1": {"Bacteria": float("nan"), "Eukaryota": 100},
                "S-2": {"Bacteria": float("inf"), "Eukaryota": 200},
            },
            ["Bacteria", "Eukaryota"],
            "Non-finite regression test",
        )
        self.assertIn("Non-finite regression test", chart)
        self.assertNotIn(">nan<", chart.lower())
        self.assertNotIn(">inf<", chart.lower())

    def test_renders_summary_and_optional_figure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "samples.tsv").write_text(
                "sample\tcondition\tLatitude [degrees_north]\tLongitude [degrees_east]\tDepth (m)\n"
                "S_1\ttest\t-39.49\t99.39\t2\n",
                encoding="utf-8",
            )
            (root / "split.tsv").write_text(
                "sample\tprok_seqs_split\teuk_seqs_split\ttotal_seqs_split\teukfrac_split\n"
                "S_1\t700\t300\t1000\t0.3\n",
                encoding="utf-8",
            )
            stats = (
                "sample-id\tinput\tfiltered\tdenoised\tmerged\tnon-chimeric\n"
                "S-1\t{input}\t{filtered}\t{denoised}\t{merged}\t{final}\n"
            )
            (root / "16.tsv").write_text(
                stats.format(input=700, filtered=650, denoised=640, merged=610, final=600),
                encoding="utf-8",
            )
            (root / "18.tsv").write_text(
                stats.format(input=300, filtered=250, denoised=220, merged="", final=200),
                encoding="utf-8",
            )
            (root / "S_1.qc.txt").write_text(
                "Total read pairs processed:              1,250\n"
                "Pairs written (passing filters):          1,100 (88.0%)\n",
                encoding="utf-8",
            )
            (root / "long.tsv").write_text(
                "SampleID\tDomain\tPhylum\tDivision\tClass\tSequence_Type\t"
                "plastid_16S_rRNA\tASV_hash\tCorrected_Sequence_Counts\n"
                "S-1\tBacteria\tProteobacteria\t\t\tProkaryotic_16S\tno\tA1\t600\n"
                "S-1\tEukaryota\t\tDinoflagellata\t\tEukaryote_18S\tno\tA2\t200\n"
                "S-1\tEukaryota\t\tChloroplastida\t\tChloroplast_16S\tyes\tA3\t100\n"
                "S-1\tEukaryota\t\t\tMitochondria\tEukaryote_18S\tno\tA4\t50\n"
                "S-1\tUnassigned\t\t\t\tUnassigned\tno\tA5\t25\n",
                encoding="utf-8",
            )
            # A tiny valid PNG is sufficient to verify base64 embedding.
            png = root / "figure.png"
            png.write_bytes(
                b"\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR\x00\x00\x00\x01\x00\x00\x00\x01"
                b"\x08\x06\x00\x00\x00\x1f\x15\xc4\x89"
            )
            (root / "corrected.tsv").write_text(
                "SampleID\tDomain\tisd_1_recovery_ratio\trecovery_mean\t"
                "Copies_BP_recovery_ratio\n"
                "S-1\tBacteria\t0.25\t0.25\t2400\n"
                "S-1\tEukaryota\t0.25\t0.25\t800\n",
                encoding="utf-8",
            )
            output = root / "report.html"
            REPORT.render_report(
                {
                    "studyName": "test-run",
                    "intstds": ["BP"],
                    "trunclens": {"truncR1": 220, "truncR2": 180},
                    "dada2": {
                        "prokaryotes": {"max_ee_f": 3.5},
                        "eukaryotes": {"max_ee": 4.5},
                    },
                },
                {
                    "samples": root / "samples.tsv",
                    "split_summary": root / "split.tsv",
                    "stats16s": root / "16.tsv",
                    "stats18s": root / "18.tsv",
                    "cutadapt_qc": [root / "S_1.qc.txt"],
                    "long_data": root / "long.tsv",
                    "internal_standard_figures": [png],
                    "internal_standard_table": [root / "corrected.tsv"],
                },
                output,
            )
            rendered = output.read_text(encoding="utf-8")
            self.assertIn("test-run", rendered)
            self.assertIn("S-1", rendered)
            self.assertIn("Filtering loss", rendered)
            self.assertIn("Denoising loss", rendered)
            self.assertIn("Pair-merging loss", rendered)
            self.assertIn("Chimera-removal loss", rendered)
            self.assertIn("Reads before filtering and quality control", rendered)
            self.assertIn("S-1: 1,250", rendered)
            self.assertIn("Primer trimming and BBsplit assignment", rendered)
            self.assertIn("BBsplit unassigned", rendered)
            self.assertIn("150 (12.0%)", rendered)
            self.assertIn("100 (9.1%)", rendered)
            self.assertIn("Reads retained after DADA2", rendered)
            self.assertIn("S-1: 800", rendered)
            self.assertIn("Not applicable", rendered)
            self.assertIn("Proteobacteria", rendered)
            self.assertIn("Interactive taxonomy bar plot", rendered)
            self.assertIn('id="taxonomy-explorer-data"', rendered)
            self.assertIn('id="taxonomy-plot-field"', rendered)
            self.assertIn('id="taxonomy-rank"', rendered)
            self.assertIn("justify-content:flex-start", rendered)
            self.assertIn("top:308px", rendered)
            self.assertIn('class="frozen-y-axis"', rendered)
            self.assertIn('axis.className = "taxonomy-y-axis"', rendered)
            self.assertIn("#332288", rendered)
            self.assertIn(
                '"metadataFields":["SampleID","Condition","Latitude","Longitude","Depth"]',
                rendered,
            )
            self.assertIn('"Phylum":{"Proteobacteria":600.0,', rendered)
            self.assertIn("Sequence assignments", rendered)
            self.assertIn("total 16S", rendered)
            self.assertIn("750</strong><span>total 16S", rendered)
            self.assertIn("100</strong><span>chloroplast 16S", rendered)
            self.assertIn("50</strong><span>mitochondrial 16S", rendered)
            self.assertIn("25</strong><span>unassigned", rendered)
            self.assertIn("Sequence-assignment counts", rendered)
            self.assertIn("Recovery ratios by sample (log10 scale)", rendered)
            self.assertIn("Domain abundance by correction method", rendered)
            self.assertIn("S-1 — BP: 0.25", rendered)
            self.assertIn("trunclens.truncR1", rendered)
            self.assertIn("Effective DADA2 settings used", rendered)
            self.assertIn("max_ee_f", rendered)
            self.assertIn("Maximum expected errors in a forward read", rendered)
            self.assertIn("<code>3.5</code>", rendered)
            self.assertIn("<code>4.5</code>", rendered)
            self.assertIn("min_overlap", rendered)
            self.assertIn("n_reads_learn", rendered)
            self.assertNotIn("Median base-quality profiles", rendered)
            self.assertNotIn("Where reads were retained or lost", rendered)

            explorer = REPORT.build_taxonomy_explorer_data(
                REPORT.read_tsv(root / "samples.tsv"),
                REPORT.read_tsv(root / "long.tsv"),
                "Corrected_Sequence_Counts",
            )
            self.assertEqual(
                explorer["metadataFields"],
                ["SampleID", "Condition", "Latitude", "Longitude", "Depth"],
            )
            self.assertIn("Phylum", explorer["ranks"])
            self.assertEqual(explorer["samples"]["S-1"]["metadata"]["SampleID"], "S-1")
            self.assertEqual(explorer["samples"]["S-1"]["metadata"]["Condition"], "test")
            self.assertEqual(explorer["samples"]["S-1"]["metadata"]["Latitude"], "-39.49")
            self.assertEqual(
                explorer["samples"]["S-1"]["taxonomy"]["Phylum"]["Proteobacteria"],
                600,
            )

            domain_chart = REPORT.svg_composition(
                {"S-1": {"Bacteria": 60, "Eukaryota": 40}},
                "Domain composition by sample",
            )
            self.assertGreater(
                domain_chart.index('class="chart-legend"'),
                domain_chart.index("</svg>"),
            )


if __name__ == "__main__":
    unittest.main()
