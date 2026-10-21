#!/bin/bash
#SBATCH --job-name=finalann_sprot
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=0-1
#SBATCH --output=/work1/sponge/lht-njlan-final-annotation-20261007/logs/swissprot_%A_%a.out
#SBATCH --error=/work1/sponge/lht-njlan-final-annotation-20261007/logs/swissprot_%A_%a.err
set -euo pipefail

ROOT=/work1/sponge/lht-njlan-final-annotation-20261007
SAMPLES=(lht04 njlan)
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
INPUT="$ROOT/$SAMPLE/04_protein_busco/${SAMPLE}_final.longest_isoform.aa"
OUT="$ROOT/$SAMPLE/06_swissprot"
DB=/work1/sponge/sci_data/databases/uniprot_sprot.dmnd
ENV=/work1/sponge/conda_envs/eggnog_mapper_2.1.13
ALL="$OUT/${SAMPLE}_final_vs_swissprot.tsv"
BEST="$OUT/${SAMPLE}_final_vs_swissprot.best.tsv"

test -s "$INPUT"
test -s "$DB"
test -x "$ENV/bin/diamond"
test ! -e "$ALL"
mkdir -p "$OUT"
"$ENV/bin/diamond" blastp --query "$INPUT" --db "$DB" --out "$ALL" \
  --outfmt 6 qseqid sseqid pident length qlen slen qcovhsp scovhsp evalue bitscore stitle \
  --evalue 1e-5 --query-cover 30 --max-target-seqs 5 --more-sensitive \
  --threads "$SLURM_CPUS_PER_TASK"
test -s "$ALL"
awk '!seen[$1]++' "$ALL" > "$BEST"
test -s "$BEST"
gzip "$ALL"
