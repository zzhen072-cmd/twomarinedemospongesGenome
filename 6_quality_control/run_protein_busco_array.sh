#!/bin/bash
#SBATCH --job-name=finalann_pbusco
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=108G
#SBATCH --time=1-00:00:00
#SBATCH --array=0-1
#SBATCH --output=/work1/sponge/lht-njlan-final-annotation-20261007/logs/protein_busco_%A_%a.out
#SBATCH --error=/work1/sponge/lht-njlan-final-annotation-20261007/logs/protein_busco_%A_%a.err
set -euo pipefail

ROOT=/work1/sponge/lht-njlan-final-annotation-20261007
SAMPLES=(lht04 njlan)
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
RESULT="$ROOT/$SAMPLE/03_braker4/output/${SAMPLE}_final/results"
OUT="$ROOT/$SAMPLE/04_protein_busco"
ALL="$OUT/${SAMPLE}_final.all_isoforms.aa"
LONGEST="$OUT/${SAMPLE}_final.longest_isoform.aa"
BUSCO_RUN="${SAMPLE}_final_longest_metazoa"
LINEAGE=/public/home/sponge/metazoa_odb10

test -s "$RESULT/braker.aa.gz"
test -s "$LINEAGE/dataset.cfg"
test ! -e "$OUT/.complete"
test ! -e "$OUT/$BUSCO_RUN"
mkdir -p "$OUT"
gzip -dc "$RESULT/braker.aa.gz" > "$ALL"

awk '
BEGIN { RS=">"; ORS="" }
NR>1 {
  n=split($0, lines, "\n")
  id=lines[1]
  sub(/[[:space:]].*$/, "", id)
  gene=id
  sub(/\.t[0-9]+$/, "", gene)
  seq=""
  for (i=2; i<=n; i++) {
    gsub(/[[:space:]]/, "", lines[i])
    seq=seq lines[i]
  }
  if (!(gene in seen)) {
    seen[gene]=1
    order[++ng]=gene
  }
  if (length(seq) > bestlen[gene]) {
    bestlen[gene]=length(seq)
    bestid[gene]=id
    bestseq[gene]=seq
  }
}
END {
  for (j=1; j<=ng; j++) {
    gene=order[j]
    print ">" bestid[gene] "\n"
    for (p=1; p<=length(bestseq[gene]); p+=60)
      print substr(bestseq[gene], p, 60) "\n"
  }
}
' "$ALL" > "$LONGEST"

ALL_N=$(grep -c '^>' "$ALL")
LONG_N=$(grep -c '^>' "$LONGEST")
test "$ALL_N" -gt 0
test "$LONG_N" -gt 0
test "$LONG_N" -le "$ALL_N"
printf 'all_isoforms\t%s\nrepresentative_longest_isoforms\t%s\n' "$ALL_N" "$LONG_N" \
  > "$OUT/isoform_stats.tsv"

export PATH=/public/home/sponge/.conda/envs/gzx_qc/bin:/public/home/sponge/.conda/envs/busco_env/bin:$PATH
busco --version
busco -i "$LONGEST" -o "$BUSCO_RUN" -m proteins -l "$LINEAGE" \
  -c "$SLURM_CPUS_PER_TASK" --offline --out_path "$OUT"
test -n "$(find "$OUT/$BUSCO_RUN" -maxdepth 1 -name 'short_summary*.txt' -print -quit)"
touch "$OUT/.complete"
