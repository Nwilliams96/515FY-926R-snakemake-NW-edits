import runpy
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


SCRIPT = (
    Path(__file__).parents[1]
    / "workflow"
    / "scripts"
    / "prepare_internal_standard_fastas.py"
)
RULES = (
    Path(__file__).parents[1]
    / "workflow"
    / "rules"
    / "05-internal-standard-correction.smk"
)


class PrepareInternalStandardFastasTest(unittest.TestCase):
    def run_script(
        self, root, rows, configured_ids, sequence_column="full_16S_sequence"
    ):
        table = root / "internal_stds.tsv"
        table.write_text(
            "internal_std_ID\trRNA_copy_number\tgenome_len_bp\t"
            + sequence_column
            + "\n"
            + "".join(
                f"{standard_id}\t{copies}\t{genome}\t{sequence}\n"
                for standard_id, copies, genome, sequence in rows
            ),
            encoding="utf-8",
        )
        output_paths = [root / f"{standard_id}.fasta" for standard_id in configured_ids]
        snakemake = SimpleNamespace(
            input=[str(table)],
            output=SimpleNamespace(fastas=[str(path) for path in output_paths]),
            params=SimpleNamespace(
                standard_ids=configured_ids,
                forward_primer="GTGYCAGCMGCCGCGGTAA",
                reverse_primer="CCGYCAATTYMTTTRAGTTT",
                trunc_r1=220,
                trunc_r2=180,
            ),
        )

        runpy.run_path(str(SCRIPT), init_globals={"snakemake": snakemake})
        return output_paths

    def test_creates_one_fasta(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            outputs = self.run_script(
                root,
                [("Only-Standard", 5, 100, "ACGT")],
                ["Only-Standard"],
            )
            self.assertEqual(
                outputs[0].read_text(encoding="utf-8"),
                ">Only-Standard\nACGT\n",
            )

    def test_creates_more_than_three_fastas_in_configured_order(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            rows = [
                ("Custom-A", 5, 100, "ACGT"),
                ("Custom_B", 2, 200, "NNAA"),
                ("SpikeIn3", 1, 300, "RYGC"),
                ("Fourth.std", 3, 400, "BDHV"),
            ]
            configured_ids = [row[0] for row in rows]
            outputs = self.run_script(root, rows, configured_ids)

            self.assertEqual(len(outputs), 4)
            self.assertEqual(
                outputs[-1].read_text(encoding="utf-8"),
                ">Fourth.std\nBDHV\n",
            )

    def test_adds_pipeline_shaped_reference_for_concatenated_18s(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            forward_primer = "GTGCCAGCAGCCGCGGTAA"
            reverse_site = "AAACTTAAAGGAATTGACGG"
            insert = "A" * 220 + "G" * 50 + "C" * 180
            full_sequence = "TTTT" + forward_primer + insert + reverse_site + "TTTT"

            outputs = self.run_script(
                root,
                [("Euk-Standard", 2, 1000, full_sequence)],
                ["Euk-Standard"],
                sequence_column="full_SSU_sequence",
            )

            self.assertEqual(
                outputs[0].read_text(encoding="utf-8"),
                ">Euk-Standard\n"
                + full_sequence
                + "\n>Euk-Standard__18S_concatenated_220_180\n"
                + "A" * 220
                + "C" * 180
                + "\n",
            )

    def test_legacy_full_16s_sequence_column_remains_supported(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            outputs = self.run_script(
                root,
                [("Legacy", 1, 100, "ACGT")],
                ["Legacy"],
            )
            self.assertEqual(outputs[0].read_text(), ">Legacy\nACGT\n")

    def test_identification_searches_both_marker_gene_outputs(self):
        rules = RULES.read_text(encoding="utf-8")
        self.assertIn('latestseqs_16s="results/02-proks/', rules)
        self.assertIn('latestseqs_18s="results/02-euks/', rules)
        self.assertIn("sort -u > {output.asvs:q}", rules)


if __name__ == "__main__":
    unittest.main()
