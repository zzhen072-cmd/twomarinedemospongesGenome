#!/bin/bash
#SBATCH --job-name=finalann_gbusco
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=108G
#SBATCH --time=1-00:00:00
#SBATCH --array=0-1
#SBATCH --output=/work1/sponge/lht-njlan-final-annotation-20261007/logs/genome_busco_%A_%a.out
#SBATCH --error=/work1/sponge/lht-njlan-final-annotation-20261007/logs/genome_busco_%A_%a.err
set -euo pipefail

ROOT=/work1/sponge/lht-njlan-final-annotation-20261007
SAMPLES=(lht04 njlan)
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
FASTA="$ROOT/$SAMPLE/00_input/${SAMPLE}_final.shortids.fasta"
OUT="$ROOT/$SAMPLE/01_genome_busco"
RUN="${SAMPLE}_final_genome"
LINEAGE=/public/home/sponge/metazoa_odb10
AUG_SOURCE=/public/home/sponge/.conda/envs/busco_env/config
AUG_PRIVATE="$OUT/augustus_config"

test -s "$ROOT/.inputs_prepared"
test -s "$FASTA"
test -s "$LINEAGE/dataset.cfg"
test -d "$AUG_SOURCE"
test ! -e "$OUT/$RUN"
mkdir -p "$OUT"
if [ ! -d "$AUG_PRIVATE" ]; then
  cp -a "$AUG_SOURCE" "$AUG_PRIVATE"
fi

export PATH=/public/home/sponge/.conda/envs/gzx_qc/bin:/public/home/sponge/.conda/envs/busco_env/bin:$PATH
export AUGUSTUS_CONFIG_PATH="$AUG_PRIVATE"
busco --version
augustus --version 2>&1 | head -3
busco -i "$FASTA" -o "$RUN" -m genome -l "$LINEAGE" -c "$SLURM_CPUS_PER_TASK" \
  --offline --augustus --out_path "$OUT"
test -s "$OUT/$RUN/short_summary.txt" || test -n "$(find "$OUT/$RUN" -maxdepth 1 -name 'short_summary*.txt' -print -quit)"
