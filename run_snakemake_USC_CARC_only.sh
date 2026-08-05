#!/bin/bash
#
#SBATCH --account=fuhrman_1138
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=logs/slurm_snakemake_%j.out
#SBATCH --error=logs/slurm_snakemake_%j.err

set -euo pipefail
echo "Started at: $(date) on $(hostname)"
mkdir -p logs

# Prevent OpenBLAS / MKL from oversubscribing CPUs on CARC
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Ensure conda is available in this shell and activate your snakemake env
# (this uses the conda shims so it works with modern conda installs)
eval "$(conda shell.bash hook)"
conda activate snakemake

echo "Using snakemake: $(which snakemake) ; version: $(snakemake --version || true)"

# Optional target passed to the script becomes the Snakemake target.
# Example: sbatch run_snakemake_USC_CARC_only.sh results/02-proks/sample-metadata.tsv
TARGET="${1:-}"

# Snakemake invocation
# - --cores: number of cores for Snakemake controller (matches SBATCH cpus-per-task)
# - --use-conda: create/use conda envs declared in workflow
# - --rerun-incomplete: pick up partial outputs
# - --latency-wait: wait for N seconds when a file appears missing on shared filesystems
snakemake \
  --snakefile workflow/Snakefile \
  --cores 1 \
  --use-conda \
  --rerun-incomplete \
  --latency-wait 60 \
  ${TARGET}
rc=$?

echo "Snakemake finished with exit code ${rc} at: $(date)"
exit ${rc}
