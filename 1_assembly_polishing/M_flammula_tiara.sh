#!/bin/bash
#SBATCH -J lan_tiara
#SBATCH -p normal
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -n 30
#SBATCH --mem=108G
#SBATCH -o /work1/sponge/lan_reassembly_compare/13_tiara/tiara-%j.out
#SBATCH -e /work1/sponge/lan_reassembly_compare/13_tiara/tiara-%j.err
set -euo pipefail

ASSEMBLY=/work1/sponge/lan_reassembly_compare/10_fcs_gx/lan_fcs_adaptor.fasta
OUT=/work1/sponge/lan_reassembly_compare/13_tiara
TIARA=/public/home/sponge/envs/tiara/bin/tiara
test -s "$ASSEMBLY"
test -x "$TIARA"
test ! -e "$OUT/lan.tiara.tsv"
cd "$OUT"
echo "START $(date)"
"$TIARA" -i "$ASSEMBLY" -o "$OUT/lan.tiara.tsv" -t 30 --to_fasta pro pla mit \
  > "$OUT/lan.tiara.stdout.log" 2> "$OUT/lan.tiara.stderr.log"
test -s "$OUT/lan.tiara.tsv"
echo "FINISH $(date)"
