#!/bin/bash
#SBATCH -J lht_tiara
#SBATCH -p normal
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -n 30
#SBATCH --mem=108G
#SBATCH -o /work1/sponge/lht_reassembly_compare/13_tiara/tiara-%j.out
#SBATCH -e /work1/sponge/lht_reassembly_compare/13_tiara/tiara-%j.err
set -euo pipefail

ASSEMBLY=/work1/sponge/lht_reassembly_compare/10_fcs_gx/lht_fcs_adaptor.fasta
OUT=/work1/sponge/lht_reassembly_compare/13_tiara
TIARA=/public/home/sponge/envs/tiara/bin/tiara
test -s "$ASSEMBLY"
test -x "$TIARA"
test ! -e "$OUT/lht.tiara.tsv"
cd "$OUT"
echo "START $(date)"
"$TIARA" -i "$ASSEMBLY" -o "$OUT/lht.tiara.tsv" -t 30 --to_fasta pro pla mit \
  > "$OUT/lht.tiara.stdout.log" 2> "$OUT/lht.tiara.stderr.log"
test -s "$OUT/lht.tiara.tsv"
echo "FINISH $(date)"
