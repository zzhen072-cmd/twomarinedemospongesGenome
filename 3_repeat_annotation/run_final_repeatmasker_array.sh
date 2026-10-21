#!/bin/bash
#SBATCH --job-name=finalann_rm
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=108G
#SBATCH --time=2-00:00:00
#SBATCH --array=0-1
#SBATCH --output=/work1/sponge/lht-njlan-final-annotation-20261007/logs/repeatmasker_%A_%a.out
#SBATCH --error=/work1/sponge/lht-njlan-final-annotation-20261007/logs/repeatmasker_%A_%a.err
set -euo pipefail

ROOT=/work1/sponge/lht-njlan-final-annotation-20261007
SAMPLES=(lht04 njlan)
LIBRARIES=(
  /work1/sponge/genome_annotation/lht_flye118_manual/repeatmodeler/lht_db-families.fa
  /work1/sponge/genome_annotation/lan_flye118_manual/repeatmodeler/lan_db-families.fa
)
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
LIB=${LIBRARIES[$SLURM_ARRAY_TASK_ID]}
FASTA="$ROOT/$SAMPLE/00_input/${SAMPLE}_final.shortids.fasta"
OUT="$ROOT/$SAMPLE/02_repeatmasker"
MASKED="$OUT/${SAMPLE}_final.shortids.fasta.masked"

test -s "$ROOT/.inputs_prepared"
test -s "$FASTA"
test -s "$LIB"
test ! -e "$MASKED"
mkdir -p "$OUT"
export PATH=/work1/sponge/conda_envs/repeatmodeler_ltr/bin:$PATH
RepeatMasker -pa 15 -xsmall -gff -lib "$LIB" -dir "$OUT" "$FASTA"
test -s "$MASKED"
test -s "$OUT/${SAMPLE}_final.shortids.fasta.out"
test -s "$OUT/${SAMPLE}_final.shortids.fasta.out.gff"
