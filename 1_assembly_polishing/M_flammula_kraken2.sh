#!/bin/bash
#SBATCH -J lan_kraken2
#SBATCH -p normal
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 30
#SBATCH --mem=108G
#SBATCH -o /work1/sponge/lan_reassembly_compare/11_kraken2/kraken2-%j.out
#SBATCH -e /work1/sponge/lan_reassembly_compare/11_kraken2/kraken2-%j.err
set -eo pipefail

ASSEMBLY=/work1/sponge/lan_reassembly_compare/10_fcs_gx/lan_fcs_adaptor.fasta
OUT=/work1/sponge/lan_reassembly_compare/11_kraken2
DB=/public/home/sponge/sci_data/S07/software/standard16
source ~/.bashrc
conda activate kraken2_env
set -u
test -s "$ASSEMBLY"
test -s "$DB/hash.k2d"
test ! -e "$OUT/lan.kraken2.out"
echo "START $(date)"
kraken2 --db "$DB" --threads 30 --confidence 0.10 \
  --report "$OUT/lan.kraken2.report.tsv" \
  --output "$OUT/lan.kraken2.out" \
  --classified-out "$OUT/lan.kraken2.classified.fasta" \
  --unclassified-out "$OUT/lan.kraken2.unclassified.fasta" \
  "$ASSEMBLY" > "$OUT/lan.kraken2.stdout.log" 2> "$OUT/lan.kraken2.stderr.log"
test -s "$OUT/lan.kraken2.report.tsv"
test -s "$OUT/lan.kraken2.out"
awk -F '\t' 'BEGIN{OFS="\t"}{count[$1]++}END{for(s in count)print s,count[s]}' \
  "$OUT/lan.kraken2.out" | sort > "$OUT/lan.kraken2.status_counts.tsv"
echo "FINISH $(date)"
