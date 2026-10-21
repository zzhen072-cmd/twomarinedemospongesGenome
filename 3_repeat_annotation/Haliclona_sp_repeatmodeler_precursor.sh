#!/bin/bash
# RepeatModeler was run on the 856-contig curation precursor, not the final
# 844-contig release. The resulting library was reused for final RepeatMasker.
#SBATCH --job-name=lht_RM2LTR
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=108G
#SBATCH --time=7-00:00:00
#SBATCH --output=/work1/sponge/genome_annotation/lht_flye118_manual/logs/repeatmodeler_%j.out
#SBATCH --error=/work1/sponge/genome_annotation/lht_flye118_manual/logs/repeatmodeler_%j.err
set -euo pipefail

BASE=/work1/sponge/genome_annotation/lht_flye118_manual
GENOME=$BASE/input/lht_shortids.fasta
OUT=$BASE/repeatmodeler
ENV=/work1/sponge/conda_envs/repeatmodeler_ltr
test -s "$GENOME"
test -x "$ENV/bin/BuildDatabase"
test -x "$ENV/bin/RepeatModeler"
export PATH="$ENV/bin:$PATH"
cd "$OUT"
if [ ! -s lht_db.nsq ]; then
  BuildDatabase -name lht_db "$GENOME"
fi
RepeatModeler -database lht_db -threads "$SLURM_CPUS_PER_TASK" -LTRStruct
test -s lht_db-families.fa
echo REPEATMODELER_DONE
