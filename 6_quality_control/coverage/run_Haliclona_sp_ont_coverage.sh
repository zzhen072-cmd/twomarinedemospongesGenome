#!/bin/bash
#SBATCH --job-name=lht844_ontmap
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=108G
#SBATCH --time=24:00:00
#SBATCH --output=/work1/sponge/lht-njlan-final-coverage-20261009/logs/lht844_ontmap.%j.out
#SBATCH --error=/work1/sponge/lht-njlan-final-coverage-20261009/logs/lht844_ontmap.%j.err

set -euo pipefail
export LC_ALL=C

module load apps/miniforge3/25.3.1-0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gzx_qc

ROOT=/work1/sponge/lht-njlan-final-coverage-20261009
SAMPLE=LHT-04
ASM=/work1/sponge/lht_reassembly_compare/15_read_validated_curated/LHT_curated_primary.fasta
READS=/public/home/sponge/sci_data/lht/00_input/LHT04.ont.fq.gz
FINAL_OUT="$ROOT/$SAMPLE"
OUT="$ROOT/.${SAMPLE}.work.${SLURM_JOB_ID}"
BAM="$OUT/LHT04.final844.ONT.sorted.bam"
BED="$OUT/whole_contigs.bed"
DEFAULT_PREFIX="$OUT/LHT04.final844.ONT.default"
PRIMARY_PREFIX="$OUT/LHT04.final844.ONT.primary"
PRIMARY_Q20_PREFIX="$OUT/LHT04.final844.ONT.primary_q20"
EXPECTED_FASTA_SHA=026d8b2e3e9888de0763d42ce226fcf7c63062a145e6ad8478e2e8c1725b2499
EXPECTED_CONTIGS=844
EXPECTED_SPAN=115734075
EXPECTED_READS=7319340
EXPECTED_READ_BASES=16987937488
EXPECTED_READS_SHA=d75566c822e74e5807eb2b218de4bb8f585ef68a5600d6062bd7ecb4f9929a4b

for tool in minimap2 samtools mosdepth seqkit sha256sum python; do
    command -v "$tool" >/dev/null
done
test -s "$ASM"
test -s "$READS"
gzip -t "$READS"
if [[ -e "$FINAL_OUT" ]]; then
    echo "Refusing to overwrite existing output: $FINAL_OUT" >&2
    exit 1
fi
mkdir "$OUT"
mkdir "$OUT/tmp"

actual_fasta_sha=$(sha256sum "$ASM" | awk '{print $1}')
[[ "$actual_fasta_sha" == "$EXPECTED_FASTA_SHA" ]]
printf '%s  %s\n' "$actual_fasta_sha" "$ASM" > "$OUT/assembly.sha256"

# The raw-read checksum is deliberately recomputed for the final provenance record.
actual_reads_sha=$(sha256sum "$READS" | awk '{print $1}')
[[ "$actual_reads_sha" == "$EXPECTED_READS_SHA" ]]
printf '%s  %s\n' "$actual_reads_sha" "$READS" > "$OUT/reads.sha256"
seqkit stats -T -j 8 "$READS" > "$OUT/reads.seqkit.tsv"
awk -F '\t' -v n="$EXPECTED_READS" -v b="$EXPECTED_READ_BASES" \
    'NR==2 {if ($4!=n || $5!=b) exit 1; ok=1} END{if(!ok) exit 1}' \
    "$OUT/reads.seqkit.tsv"

minimap2 --version > "$OUT/minimap2.version.txt"
samtools --version > "$OUT/samtools.version.txt"
mosdepth --version > "$OUT/mosdepth.version.txt" 2>&1
seqkit version > "$OUT/seqkit.version.txt" 2>&1
python --version > "$OUT/python.version.txt" 2>&1

printf 'sample_id\t%s\nassembly\t%s\nreads\t%s\nbam\t%s\nminimap2\t-ax map-ont -t 22 -I 20G\ncoverage_default\tmosdepth -n -F 1796 -Q 0 --by whole-contig BED --thresholds 1,5,10,20\ncoverage_primary\tmosdepth -n -F 3844 -Q 0 --by whole-contig BED --thresholds 1,5,10,20\ncoverage_primary_q20\tmosdepth -n -F 3844 -Q 20 --by whole-contig BED --thresholds 1,5,10,20\n' \
    "$SAMPLE" "$ASM" "$READS" "$FINAL_OUT/LHT04.final844.ONT.sorted.bam" > "$OUT/input_manifest.tsv"

minimap2 -ax map-ont -t 22 -I 20G "$ASM" "$READS" \
    2> "$OUT/minimap2.log" \
    | samtools sort -@ 4 -m 2G -T "$OUT/tmp/LHT04.final844.sort" -o "$BAM" -

samtools quickcheck -v "$BAM"
samtools index -@ 4 "$BAM"
samtools flagstat -@ 4 "$BAM" > "$OUT/LHT04.final844.ONT.flagstat.txt"
primary_records=$(samtools view -@ 4 -c -F 2816 "$BAM")
primary_mapped=$(samtools view -@ 4 -c -F 2820 "$BAM")
primary_mapped_q20=$(samtools view -@ 4 -c -F 2820 -q 20 "$BAM")
printf 'metric\tcount\nprimary_records_F2816\t%s\nprimary_mapped_F2820\t%s\nprimary_mapped_q20_F2820\t%s\n' \
    "$primary_records" "$primary_mapped" "$primary_mapped_q20" \
    > "$OUT/mapping_counts.tsv"
[[ "$primary_records" -eq "$EXPECTED_READS" ]]
samtools idxstats "$BAM" > "$OUT/LHT04.final844.ONT.idxstats.tsv"
awk 'BEGIN{OFS="\t"} $1!="*"{print $1,0,$2}' \
    "$OUT/LHT04.final844.ONT.idxstats.tsv" > "$BED"

[[ "$(wc -l < "$BED")" -eq "$EXPECTED_CONTIGS" ]]
awk -F '\t' -v expected="$EXPECTED_SPAN" \
    '{s+=$3-$2} END{if(s!=expected) exit 1}' "$BED"

# Compatibility output reproduces the earlier NJ-lan mosdepth default, which
# includes supplementary alignments but excludes unmapped, secondary,
# QC-failed and duplicate records (mask 1796).
mosdepth -n -t 4 -F 1796 --by "$BED" --thresholds 1,5,10,20 \
    "$DEFAULT_PREFIX" "$BAM"

# Article-primary output additionally excludes supplementary alignments
# (mask 3844) and is the primary coverage definition used in the manuscript.
mosdepth -n -t 4 -F 3844 --by "$BED" --thresholds 1,5,10,20 \
    "$PRIMARY_PREFIX" "$BAM"

# A high-confidence primary-only diagnostic additionally requires MAPQ >=20.
mosdepth -n -t 4 -F 3844 -Q 20 --by "$BED" --thresholds 1,5,10,20 \
    "$PRIMARY_Q20_PREFIX" "$BAM"

for path in \
    "$DEFAULT_PREFIX.regions.bed.gz" \
    "$DEFAULT_PREFIX.thresholds.bed.gz" \
    "$PRIMARY_PREFIX.regions.bed.gz" \
    "$PRIMARY_PREFIX.thresholds.bed.gz" \
    "$PRIMARY_Q20_PREFIX.regions.bed.gz" \
    "$PRIMARY_Q20_PREFIX.thresholds.bed.gz"; do
    test -s "$path"
done

printf 'bam_path\t%s\nbam_size_bytes\t%s\nbam_mtime_epoch\t%s\n' \
    "$FINAL_OUT/LHT04.final844.ONT.sorted.bam" \
    "$(stat -c %s "$BAM")" "$(stat -c %Y "$BAM")" > "$OUT/bam.stat.tsv"
bam_sha=$(sha256sum "$BAM" | awk '{print $1}')
bai_sha=$(sha256sum "$BAM.bai" | awk '{print $1}')
printf '%s  %s\n%s  %s\n' \
    "$bam_sha" "$FINAL_OUT/LHT04.final844.ONT.sorted.bam" \
    "$bai_sha" "$FINAL_OUT/LHT04.final844.ONT.sorted.bam.bai" \
    > "$OUT/bam.sha256"
date -Is > "$OUT/COMPLETED.txt"
mv -T "$OUT" "$FINAL_OUT"
echo "LHT-04 final 844-contig ONT mapping and three coverage definitions completed."
