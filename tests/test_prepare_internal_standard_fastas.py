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


class PrepareInternalStandardFastasTest(unittest.TestCase):
    def run_script(self, root, rows, configured_ids):
        table = root / "internal_stds.tsv"
        table.write_text(
            "internal_std_ID\trRNA_copy_number\tgenome_len_bp\tfull_16S_sequence\n"
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
            params=SimpleNamespace(standard_ids=configured_ids),
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


if __name__ == "__main__":
    unittest.main()
