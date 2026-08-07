import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class Dada2ConfigurationTests(unittest.TestCase):
    def test_historical_defaults_are_declared(self):
        common = (ROOT / "workflow/rules/common.smk").read_text(encoding="utf-8")
        for setting in (
            '"max_ee_f": 2.0',
            '"max_ee_r": 2.0',
            '"trunc_q": 2',
            '"min_overlap": 12',
            '"max_ee": 2.0',
            '"trunc_q": 0',
            '"pooling_method": "independent"',
            '"chimera_method": "consensus"',
            '"n_reads_learn": 1000000',
        ):
            self.assertIn(setting, common)

    def test_configurable_values_reach_both_dada2_commands(self):
        prok = (ROOT / "workflow/scripts/P03-DADA2.sh").read_text(encoding="utf-8")
        euk = (ROOT / "workflow/scripts/E08-DADA2.sh").read_text(encoding="utf-8")

        for option in (
            "--p-max-ee-f",
            "--p-max-ee-r",
            "--p-trunc-q",
            "--p-min-overlap",
            "--p-pooling-method",
            "--p-chimera-method",
            "--p-min-fold-parent-over-abundance",
            "--p-n-reads-learn",
        ):
            self.assertIn(option, prok)

        for option in (
            "--p-max-ee",
            "--p-trunc-q",
            "--p-pooling-method",
            "--p-chimera-method",
            "--p-min-fold-parent-over-abundance",
            "--p-n-reads-learn",
        ):
            self.assertIn(option, euk)


if __name__ == "__main__":
    unittest.main()
