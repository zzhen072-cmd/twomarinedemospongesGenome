#!/bin/bash
#SBATCH --job-name=lan_review_filter
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=lan_review_filter.%j.out
#SBATCH --error=lan_review_filter.%j.err

set -euo pipefail
export LC_ALL=C

module load apps/miniforge3/25.3.1-0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gzx_qc

ROOT=/work1/sponge/lan_tier2_decontam_review
INPUT=/work1/sponge/lan_reassembly_compare/14_manual_decontam/lan_manual_decontaminated.fasta
REMOVE_IDS="$ROOT/audit/remove_confirmed_contaminants.ids.txt"
OUT="$ROOT/assembly"
RETAINED="$OUT/lan_host_review_filtered.fasta"
REMOVED="$OUT/lan_removed_confirmed_contaminants.fasta"

test -s "$INPUT"
test -s "$REMOVE_IDS"
test -s "$ROOT/scripts/filter_fasta_exact.py"
if [[ -e "$ROOT/FILTER_COMPLETED.txt" ]]; then
    echo "Completed filter output already exists: $ROOT/FILTER_COMPLETED.txt" >&2
    exit 1
fi
mkdir -p "$OUT"

python "$ROOT/scripts/filter_fasta_exact.py" \
    --input-fasta "$INPUT" \
    --remove-ids "$REMOVE_IDS" \
    --retained-fasta "$RETAINED" \
    --removed-fasta "$REMOVED" \
    --audit-json "$OUT/filter_audit.json" \
    --expected-input-count 7366 \
    --expected-input-bp 302079020 \
    --expected-remove-count 249 \
    --expected-remove-bp 10189414

samtools faidx "$RETAINED"
samtools faidx "$REMOVED"
gzip -c "$RETAINED" > "$RETAINED.gz"
gzip -c "$REMOVED" > "$REMOVED.gz"
sha256sum "$RETAINED" "$RETAINED.gz" "$REMOVED" "$REMOVED.gz" \
    > "$OUT/fasta_sha256.tsv"
date -Is > "$ROOT/FILTER_COMPLETED.txt"
