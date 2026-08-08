import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class CarcEfficiencyTests(unittest.TestCase):
    def test_runner_uses_slurm_cpu_and_memory_allocations(self):
        runner = (ROOT / "run_snakemake_USC_CARC_only.sh").read_text(encoding="utf-8")
        self.assertIn("#SBATCH --cpus-per-task=8", runner)
        self.assertIn('PIPELINE_CORES="${SLURM_CPUS_PER_TASK:-1}"', runner)
        self.assertIn('--cores "${PIPELINE_CORES}"', runner)
        self.assertIn("--resources mem_mb=120000", runner)

    def test_runner_labels_jobs_from_study_name(self):
        runner = (ROOT / "run_snakemake_USC_CARC_only.sh").read_text(encoding="utf-8")
        self.assertIn("read_study_name()", runner)
        self.assertIn('JOB_LABEL="eASV-${STUDY_NAME}"', runner)
        self.assertIn('sbatch --job-name="${JOB_LABEL}" "$0" "$@"', runner)
        self.assertIn("#SBATCH --output=logs/%x_%j.out", runner)
        self.assertIn("#SBATCH --error=logs/%x_%j.err", runner)

    def test_bbsplit_receives_threads_and_realistic_memory(self):
        rules = (ROOT / "workflow/rules/01-split-16S-18S.smk").read_text(
            encoding="utf-8"
        )
        self.assertIn("rule bbsplit_prok_euk:", rules)
        self.assertIn("threads: 8", rules)
        self.assertIn("mem_mb=64000", rules)
        self.assertIn("bbsplit.sh threads={threads}", rules)

    def test_dada2_and_classifiers_receive_rule_threads(self):
        for relative_path in (
            "workflow/scripts/P03-DADA2.sh",
            "workflow/scripts/E08-DADA2.sh",
        ):
            script = (ROOT / relative_path).read_text(encoding="utf-8")
            self.assertIn("--p-n-threads ${snakemake[threads]}", script)

        for relative_path in (
            "workflow/scripts/P05-classify-eASVs.sh",
            "workflow/scripts/E10-classify-seqs.sh",
        ):
            script = (ROOT / relative_path).read_text(encoding="utf-8")
            self.assertIn("--p-n-jobs ${snakemake[threads]}", script)


if __name__ == "__main__":
    unittest.main()
