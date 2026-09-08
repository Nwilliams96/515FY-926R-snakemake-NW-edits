import runpy
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

try:
    import pandas as pd
except ImportError:  # The workflow's QIIME environment supplies pandas.
    pd = None


ROOT = Path(__file__).resolve().parents[1]


class SampleMetadataIdentifierTests(unittest.TestCase):
    def test_export_scripts_use_exported_biom_and_propagate_pipe_failures(self):
        euk_export = (ROOT / "workflow/scripts/E09-export-DADA2-results.sh").read_text()
        prok_export = (ROOT / "workflow/scripts/P04-export-DADA2-results.sh").read_text()

        self.assertIn("set -euo pipefail", euk_export)
        self.assertIn("set -euo pipefail", prok_export)
        self.assertIn(
            "biom convert -i ${snakemake_output[0]}/feature-table.biom",
            euk_export,
        )
        self.assertNotIn(
            "biom convert -i ${snakemake_input[0]}/feature-table.biom",
            euk_export,
        )

    @unittest.skipUnless(pd is not None, "pandas is supplied by the workflow environment")
    def test_numeric_sample_ids_remain_text_in_both_metadata_scripts(self):
        for script_name, stats_key, manifest_header in (
            (
                "E11-make-sample-metadata-file.py",
                "eukstats",
                "sample-id\tabsolute-filepath\n001\t/a/001.fastq\n002\t/a/002.fastq\n",
            ),
            (
                "P06-make-sample-metadata-file.py",
                "prokstats",
                "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n"
                "001\t/a/001.R1.fastq\t/a/001.R2.fastq\n"
                "002\t/a/002.R1.fastq\t/a/002.R2.fastq\n",
            ),
        ):
            with self.subTest(script=script_name), tempfile.TemporaryDirectory() as temp_dir:
                temp = Path(temp_dir)
                manifest = temp / "manifest.tsv"
                samples = temp / "samples.tsv"
                stats = temp / "stats.tsv"
                fractions = temp / "fractions.tsv"
                output = temp / "metadata.tsv"

                manifest.write_text(manifest_header, encoding="utf-8")
                samples.write_text("sample\tCondition\n001\tA\n002\tB\n", encoding="utf-8")
                stats.write_text(
                    "sample-id\tnon-chimeric\n001\t10\n002\t20\n",
                    encoding="utf-8",
                )
                fractions.write_text(
                    "sample\tPROK_reads\tEUK_reads\n001\t10\t5\n002\t20\t10\n",
                    encoding="utf-8",
                )

                inputs = {
                    "manifest": str(manifest),
                    "samplesdottsv": str(samples),
                    stats_key: str(stats),
                    "eukfracpersample": str(fractions),
                }
                runpy.run_path(
                    ROOT / "workflow/scripts" / script_name,
                    init_globals={"snakemake": SimpleNamespace(input=inputs, output=[str(output)])},
                )

                result = pd.read_csv(output, sep="\t", dtype={"sample-id": str})
                self.assertEqual(result["sample-id"].tolist(), ["001", "002"])


if __name__ == "__main__":
    unittest.main()
