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
                "BP\t5\t100\tACGT\n"
                "DR\t2\t200\tNNAA\n"
                "TT\t1\t300\tRYGC\n",
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
                    "intstd1name": "BP",
                    "intstd2name": "DR",
                    "intstd3name": "TT",
                },
            )

            runpy.run_path(str(SCRIPT), init_globals={"snakemake": snakemake})

            self.assertEqual(
                (root / "intstd1.fasta").read_text(encoding="utf-8"),
                ">BP\nACGT\n",
            )
            self.assertEqual(
                (root / "intstd2.fasta").read_text(encoding="utf-8"),
                ">DR\nNNAA\n",
            )


if __name__ == "__main__":
    unittest.main()
