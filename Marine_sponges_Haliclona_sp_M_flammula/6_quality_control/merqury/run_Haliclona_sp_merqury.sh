#!/bin/bash
#SBATCH --job-name=LHT_merqury
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/work1/sponge/lht_njlan_merqury_20261019/LHT-04/slurm-%j.out
#SBATCH --error=/work1/sponge/lht_njlan_merqury_20261019/LHT-04/slurm-%j.err
set -euo pipefail

export PATH=/public/home/sponge/.conda/envs/lht_merqury/bin:/usr/bin:/bin
export MERQURY=/public/home/sponge/.conda/envs/lht_merqury/share/merqury
OUT=/work1/sponge/lht_njlan_merqury_20261019/LHT-04
ASSEMBLY=/work1/sponge/lht_reassembly_compare/15_read_validated_curated/LHT_curated_primary.fasta
R1=/public/home/sponge/sci_data/lht/00_input/LHT_R1.fq.gz
R2=/public/home/sponge/sci_data/lht/00_input/LHT_R2.fq.gz
cd "$OUT"
sha256sum "$ASSEMBLY" "$R1" "$R2" > input_sha256.txt
meryl --version > meryl_version.txt 2>&1
meryl k=19 count threads=16 memory=56 output LHT_reads.meryl "$R1" "$R2"
merqury.sh LHT_reads.meryl "$ASSEMBLY" LHT-04
