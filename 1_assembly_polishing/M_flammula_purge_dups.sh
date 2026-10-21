#!/bin/bash
#SBATCH --job-name=lan_flye118_purge
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=108G
#SBATCH --time=1-00:00:00
#SBATCH --output=/work1/sponge/lan_reassembly_compare/logs/flye118_purgedups_%j.out
#SBATCH --error=/work1/sponge/lan_reassembly_compare/logs/flye118_purgedups_%j.err
set -euo pipefail

export PATH=/public/home/sponge/.conda/envs/gzx_qc/bin:$PATH
D=/work1/sponge/lan_reassembly_compare/07_flye_sub118_purge_dups
ASM=/work1/sponge/lan_reassembly_compare/06_flye_sub118_nextpolish/lan_flye118_nextpolish.fasta
ONT=/public/home/sponge/sci_data/lan/00_input/Lanm.ont.fq.gz
mkdir -p "$D"
cd "$D"
test -s "$ASM"
echo "[STAGE] ONT coverage $(date)"
minimap2 -x map-ont -t 30 "$ASM" "$ONT" | gzip -c > ont.paf.gz
pbcstat ont.paf.gz
calcuts PB.stat > cutoffs
echo "[CUTOFFS] $(cat cutoffs)"
echo "[STAGE] assembly self-alignment $(date)"
split_fa "$ASM" > assembly.split.fa
minimap2 -x asm5 -DP -t 30 assembly.split.fa assembly.split.fa | gzip -c > self.paf.gz
echo "[STAGE] purge_dups $(date)"
purge_dups -2 -T cutoffs -c PB.base.cov self.paf.gz > dups.bed
get_seqs -e dups.bed "$ASM"
test -s purged.fa
ln -sfn "$D/purged.fa" "$D/lan_flye118_nextpolish_purged.fa"
echo "[DONE] $(date)"
