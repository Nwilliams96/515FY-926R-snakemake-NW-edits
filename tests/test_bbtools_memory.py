import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class BBToolsMemoryTests(unittest.TestCase):
    def test_small_per_sample_bbtools_jobs_have_fixed_heap_caps(self):
        trimming = (ROOT / "workflow/scripts/E03-bbduk-cut-reads.sh").read_text()
        fusing = (ROOT / "workflow/scripts/E04-fuse-EUKs-withoutNs.sh").read_text()

        self.assertEqual(trimming.count("bbduk.sh -Xmx2g"), 2)
        self.assertIn("repair.sh -Xmx2g", trimming)
        self.assertIn("fuse.sh -Xmx2g", fusing)

    def test_snakemake_accounts_for_each_bbtools_job_memory(self):
        rules = (ROOT / "workflow/rules/02-denoise-and-export-euk.smk").read_text()

        self.assertGreaterEqual(rules.count("mem_mb=2500"), 2)


if __name__ == "__main__":
    unittest.main()
