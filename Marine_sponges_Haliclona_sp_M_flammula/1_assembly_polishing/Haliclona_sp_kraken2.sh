#!/bin/bash
#SBATCH -J lht.kraken2_c0
#SBATCH -p normal
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 30
#SBATCH --mem=108G
#SBATCH -o /work1/sponge/lht_reassembly_compare/11_kraken2/kraken2-%j.out
#SBATCH -e /work1/sponge/lht_reassembly_compare/11_kraken2/kraken2-%j.err
set -eo pipefail

ASSEMBLY=/work1/sponge/lht_reassembly_compare/10_fcs_gx/lht_fcs_adaptor.fasta
OUT=/work1/sponge/lht_reassembly_compare/11_kraken2
DB=/public/home/sponge/sci_data/S07/software/standard16
source ~/.bashrc
conda activate kraken2_env
set -u
test -s "$ASSEMBLY"
test -s "$DB/hash.k2d"
test ! -e "$OUT/lht.kraken2_c0.out"
echo "START $(date)"
kraken2 --db "$DB" --threads 30 --confidence 0 \
  --report "$OUT/lht.kraken2_c0.report.tsv" \
  --output "$OUT/lht.kraken2_c0.out" \
  --classified-out "$OUT/lht.kraken2_c0.classified.fasta" \
  --unclassified-out "$OUT/lht.kraken2_c0.unclassified.fasta" \
  "$ASSEMBLY" > "$OUT/lht.kraken2_c0.stdout.log" 2> "$OUT/lht.kraken2_c0.stderr.log"
test -s "$OUT/lht.kraken2_c0.report.tsv"
test -s "$OUT/lht.kraken2_c0.out"
awk -F '\t' 'BEGIN{OFS="\t"}{count[$1]++}END{for(s in count)print s,count[s]}' \
  "$OUT/lht.kraken2_c0.out" | sort > "$OUT/lht.kraken2_c0.status_counts.tsv"
echo "FINISH $(date)"
