#!/bin/bash
#SBATCH --job-name=lan_review_ont
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=108G
#SBATCH --time=24:00:00
#SBATCH --output=lan_review_ont.%j.out
#SBATCH --error=lan_review_ont.%j.err

set -euo pipefail
export LC_ALL=C

module load apps/miniforge3/25.3.1-0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gzx_qc

ROOT=/work1/sponge/lan_tier2_decontam_review
ASM="$ROOT/assembly/lan_host_review_filtered.fasta"
READS=/public/home/sponge/sci_data/lan/00_input/Lanm.ont.fq.gz
OUT="$ROOT/mapping"
BAM="$OUT/lan.review_filtered.ONT.sorted.bam"
PREFIX="$OUT/lan.review_filtered.ONT"

for TOOL in minimap2 samtools mosdepth python; do
    command -v "$TOOL" >/dev/null
done
test -s "$ROOT/FILTER_COMPLETED.txt"
test -s "$ASM"
test -s "$READS"
test -s "$ROOT/scripts/make_lan_coverage_table.py"
if [[ -e "$OUT/COMPLETED.txt" ]]; then
    echo "Completed mapping output already exists: $OUT" >&2
    exit 1
fi
mkdir -p "$OUT/tmp"

printf 'assembly\t%s\nreads\t%s\nbam\t%s\ncoverage_method\tmosdepth -n --by whole-contig BED\n' \
    "$ASM" "$READS" "$BAM" > "$OUT/input_manifest.tsv"
sha256sum "$ASM" "$READS" > "$OUT/input_sha256.tsv"
minimap2 --version > "$OUT/minimap2.version.txt"
samtools --version > "$OUT/samtools.version.txt"
mosdepth --version > "$OUT/mosdepth.version.txt" 2>&1

minimap2 -ax map-ont -t 22 -I 20G "$ASM" "$READS" \
    2> "$OUT/minimap2.log" \
    | samtools sort -@ 4 -m 2G -T "$OUT/tmp/lan.review.sort" -o "$BAM" -
samtools quickcheck -v "$BAM"
samtools index -@ 4 "$BAM"
samtools flagstat -@ 4 "$BAM" > "$PREFIX.flagstat.txt"

samtools idxstats "$BAM" \
    | awk 'BEGIN{OFS="\t"} $1!="*"{print $1,0,$2}' \
    > "$OUT/whole_contigs.bed"
test "$(wc -l < "$OUT/whole_contigs.bed")" -eq 7117

mosdepth -n -t 4 --by "$OUT/whole_contigs.bed" --thresholds 1,5,10,20 \
    "$PREFIX.by_contig" "$BAM"
test -s "$PREFIX.by_contig.regions.bed.gz"
test -s "$PREFIX.by_contig.thresholds.bed.gz"

python "$ROOT/scripts/make_lan_coverage_table.py" \
    --fasta "$ASM" \
    --regions "$PREFIX.by_contig.regions.bed.gz" \
    --out "$OUT/lan_review_filtered_classification_coverage.tsv" \
    --audit "$OUT/coverage_audit.json"

date -Is > "$OUT/COMPLETED.txt"
printf 'Finished. Coverage table: %s\n' "$OUT/lan_review_filtered_classification_coverage.tsv"
