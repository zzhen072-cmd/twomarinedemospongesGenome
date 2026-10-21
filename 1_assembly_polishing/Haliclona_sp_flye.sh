#!/bin/bash
#SBATCH -J lht_flye_sub118
#SBATCH -p long
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 30
#SBATCH --mem=108G
#SBATCH -t 7-00:00:00
#SBATCH -o /work1/sponge/lht_reassembly_compare/logs/flye_sub118-%j.out
#SBATCH -e /work1/sponge/lht_reassembly_compare/logs/flye_sub118-%j.err
set -eo pipefail

env=/work1/sponge/lht_reassembly_compare/envs/flye39/bin
export PATH=$env:/public/home/sponge/.conda/envs/gzx_qc/bin:$PATH
out=/work1/sponge/lht_reassembly_compare/03_assemblies/flye_sub118
mkdir -p "$out"
"$env/flye" \
  --nano-raw /work1/sponge/lht_reassembly_compare/02_read_filtering/LHT04.longest_11.8Gb.fq.gz \
  --genome-size 160m --threads 30 --out-dir "$out"
