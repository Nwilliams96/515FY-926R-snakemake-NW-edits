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
    def test_creates_one_valid_fasta_per_configured_standard(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            table = root / "internal_stds.tsv"
            table.write_text(
                "internal_std_ID\trRNA_copy_number\tgenome_len_bp\tfull_16S_sequence\n"
                "Custom-A\t5\t100\tACGT\n"
                "Custom_B\t2\t200\tNNAA\n"
                "SpikeIn3\t1\t300\tRYGC\n",
                encoding="utf-8",
            )
            outputs = {
                slot: str(root / f"{slot}.fasta")
                for slot in ("intstd1", "intstd2", "intstd3")
            }
            snakemake = SimpleNamespace(
                input=[str(table)],
                output=outputs,
                params={
                    "intstd1name": "Custom-A",
                    "intstd2name": "Custom_B",
                    "intstd3name": "SpikeIn3",
                },
            )

            runpy.run_path(str(SCRIPT), init_globals={"snakemake": snakemake})

            self.assertEqual(
                (root / "intstd1.fasta").read_text(encoding="utf-8"),
                ">Custom-A\nACGT\n",
            )
            self.assertEqual(
                (root / "intstd2.fasta").read_text(encoding="utf-8"),
                ">Custom_B\nNNAA\n",
            )


if __name__ == "__main__":
    unittest.main()
