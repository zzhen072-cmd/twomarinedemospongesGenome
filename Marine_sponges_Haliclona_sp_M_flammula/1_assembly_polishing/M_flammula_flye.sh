#!/bin/bash
#SBATCH -J lan_flye_sub118
#SBATCH -p long
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 30
#SBATCH --mem=108G
#SBATCH -t 7-00:00:00
#SBATCH -o /work1/sponge/lan_reassembly_compare/logs/flye_sub118-%j.out
#SBATCH -e /work1/sponge/lan_reassembly_compare/logs/flye_sub118-%j.err
set -euo pipefail

# The "11.8Gb" filename is a legacy label. This file contained the complete
# NJ-lan ONT dataset (27,612,556,815 bp), not an 11.8-Gb subset.
E=/work1/sponge/lht_reassembly_compare/envs/flye39/bin
export PATH=$E:/public/home/sponge/.conda/envs/gzx_qc/bin:$PATH
O=/work1/sponge/lan_reassembly_compare/03_assemblies/flye_sub118
mkdir -p "$O"
"$E/flye" \
  --nano-raw /work1/sponge/lan_reassembly_compare/02_read_filtering/Lanm.longest_11.8Gb.fq.gz \
  --genome-size 560m --threads 30 --out-dir "$O"
