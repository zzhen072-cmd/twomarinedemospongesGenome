#!/bin/bash
#SBATCH --job-name=finalann_eggnog
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=108G
#SBATCH --time=1-00:00:00
#SBATCH --array=0-1
#SBATCH --output=/work1/sponge/lht-njlan-final-annotation-20261007/logs/eggnog_%A_%a.out
#SBATCH --error=/work1/sponge/lht-njlan-final-annotation-20261007/logs/eggnog_%A_%a.err
set -euo pipefail

ROOT=/work1/sponge/lht-njlan-final-annotation-20261007
SAMPLES=(lht04 njlan)
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
INPUT="$ROOT/$SAMPLE/04_protein_busco/${SAMPLE}_final.longest_isoform.aa"
OUT="$ROOT/$SAMPLE/05_eggnog"
DB=/work1/sponge/sci_data/databases/eggnog_v5.0.2
ENV=/work1/sponge/conda_envs/eggnog_mapper_2.1.13
PREFIX="${SAMPLE}_final.longest"

test -s "$INPUT"
test -x "$ENV/bin/emapper.py"
test -s "$DB/eggnog.db"
test -s "$DB/eggnog_proteins.dmnd"
test -s "$DB/eggnog.taxa.db"
test ! -e "$OUT/$PREFIX.emapper.annotations"
mkdir -p "$OUT/tmp"
export PATH="$ENV/bin:$PATH"
export TMPDIR="$OUT/tmp"
emapper.py -i "$INPUT" --itype proteins -m diamond \
  --cpu "$SLURM_CPUS_PER_TASK" --data_dir "$DB" \
  --dmnd_db "$DB/eggnog_proteins.dmnd" --tax_scope 33208 \
  --output "$PREFIX" --output_dir "$OUT" --temp_dir "$TMPDIR"
test -s "$OUT/$PREFIX.emapper.annotations"
