#!/bin/bash
# RepeatModeler was run on the 7,366-contig curation precursor, not the final
# 7,117-contig release. The resulting library was reused for final RepeatMasker.
#SBATCH --job-name=lan_RM2LTR
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=108G
#SBATCH --time=7-00:00:00
#SBATCH --output=/work1/sponge/genome_annotation/lan_flye118_manual/logs/repeatmodeler_%j.out
#SBATCH --error=/work1/sponge/genome_annotation/lan_flye118_manual/logs/repeatmodeler_%j.err
set -euo pipefail

BASE=/work1/sponge/genome_annotation/lan_flye118_manual
GENOME=$BASE/input/lan_shortids.fasta
OUT=$BASE/repeatmodeler
ENV=/work1/sponge/conda_envs/repeatmodeler_ltr
test -s "$GENOME"
test -x "$ENV/bin/BuildDatabase"
test -x "$ENV/bin/RepeatModeler"
export PATH="$ENV/bin:$PATH"
cd "$OUT"
if [ ! -s lan_db.nsq ]; then
  BuildDatabase -name lan_db "$GENOME"
fi
RepeatModeler -database lan_db -threads "$SLURM_CPUS_PER_TASK" -LTRStruct
test -s lan_db-families.fa
echo REPEATMODELER_DONE
