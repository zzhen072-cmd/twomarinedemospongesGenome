#!/bin/bash
#SBATCH --job-name=lan_flye118_finalsgs
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=108G
#SBATCH --time=1-00:00:00
#SBATCH --output=/work1/sponge/lan_reassembly_compare/logs/flye118_finalsgs_%j.out
#SBATCH --error=/work1/sponge/lan_reassembly_compare/logs/flye118_finalsgs_%j.err
set -euo pipefail

export PATH=/work1/sponge/lht_reassembly_compare/envs/nextdenovo39/bin:/work1/sponge/lht_reassembly_compare/envs/polish/bin:$PATH
D=/work1/sponge/lan_reassembly_compare/08_flye_sub118_final_sgs_polish
mkdir -p "$D/work"
test -s /work1/sponge/lan_reassembly_compare/07_flye_sub118_purge_dups/lan_flye118_nextpolish_purged.fa
echo "[START] $(date)"
nextPolish "$D/run.cfg"
test -s "$D/work/genome.nextpolish.fasta"
ln -sfn "$D/work/genome.nextpolish.fasta" "$D/lan_flye118_final_polished.fasta"
echo "[DONE] $(date)"
