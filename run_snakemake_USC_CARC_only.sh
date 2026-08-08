#!/bin/bash
#
#SBATCH --account=fuhrman_1138
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --output=logs/slurm_snakemake_%j.out
#SBATCH --error=logs/slurm_snakemake_%j.err

set -euo pipefail
echo "Started at: $(date) on $(hostname)"
mkdir -p logs

# Prevent libraries from taking cores outside the per-rule limits set by
# Snakemake. Multithreaded tools receive their thread counts explicitly.
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
PIPELINE_CORES="${SLURM_CPUS_PER_TASK:-1}"

SNAKEMAKE_ARGS=(
  --snakefile workflow/Snakefile
  --cores "${PIPELINE_CORES}"
  --resources mem_mb=120000
  --use-conda
  --rerun-incomplete
  --latency-wait 60
  --printshellcmds
)

if [[ -n "${TARGET}" ]]; then
  SNAKEMAKE_ARGS+=("${TARGET}")
fi

# Snakemake invocation
# - --cores: use every CPU allocated by Slurm
# - --resources: prevent concurrent rules from exceeding the node's memory
# - --use-conda: create/use conda envs declared in workflow
# - --rerun-incomplete: pick up partial outputs
# - --latency-wait: wait for N seconds when a file appears missing on shared filesystems
if snakemake "${SNAKEMAKE_ARGS[@]}"; then
  rc=0
else
  rc=$?
fi

echo "Snakemake finished with exit code ${rc} at: $(date)"
exit ${rc}
