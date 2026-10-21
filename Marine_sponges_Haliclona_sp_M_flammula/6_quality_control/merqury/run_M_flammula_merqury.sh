#!/bin/bash
#SBATCH --job-name=NJlan_merqury
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/work1/sponge/lht_njlan_merqury_20261019/NJ-lan/slurm-%j.out
#SBATCH --error=/work1/sponge/lht_njlan_merqury_20261019/NJ-lan/slurm-%j.err
set -euo pipefail

export PATH=/public/home/sponge/.conda/envs/lht_merqury/bin:/usr/bin:/bin
export MERQURY=/public/home/sponge/.conda/envs/lht_merqury/share/merqury
OUT=/work1/sponge/lht_njlan_merqury_20261019/NJ-lan
ASSEMBLY=/work1/sponge/lan_tier2_decontam_review/assembly/lan_host_review_filtered.fasta
R1=/public/home/sponge/sci_data/lan/00_input/3_R1.fq.gz
R2=/public/home/sponge/sci_data/lan/00_input/3_R2.fq.gz
cd "$OUT"
sha256sum "$ASSEMBLY" "$R1" "$R2" > input_sha256.txt
meryl --version > meryl_version.txt 2>&1
meryl k=19 count threads=16 memory=56 output NJlan_reads.meryl "$R1" "$R2"
merqury.sh NJlan_reads.meryl "$ASSEMBLY" NJ-lan
