#!/bin/bash
#SBATCH --job-name=njlan7117_cov
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --time=12:00:00
#SBATCH --output=/work1/sponge/lht-njlan-final-coverage-20261009/logs/njlan7117_cov.%j.out
#SBATCH --error=/work1/sponge/lht-njlan-final-coverage-20261009/logs/njlan7117_cov.%j.err

set -euo pipefail
export LC_ALL=C

module load apps/miniforge3/25.3.1-0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gzx_qc

ROOT=/work1/sponge/lht-njlan-final-coverage-20261009
SOURCE=/work1/sponge/lan_tier2_decontam_review/mapping
SAMPLE=NJ-lan
ASM=/work1/sponge/lan_tier2_decontam_review/assembly/lan_host_review_filtered.fasta
READS=/public/home/sponge/sci_data/lan/00_input/Lanm.ont.fq.gz
BAM="$SOURCE/lan.review_filtered.ONT.sorted.bam"
FINAL_OUT="$ROOT/$SAMPLE"
OUT="$ROOT/.${SAMPLE}.work.${SLURM_JOB_ID}"
BED="$OUT/whole_contigs.bed"
PRIMARY_PREFIX="$OUT/NJ-lan.final7117.ONT.primary"
PRIMARY_Q20_PREFIX="$OUT/NJ-lan.final7117.ONT.primary_q20"
EXPECTED_FASTA_SHA=cea9b34517164f02103111a3a1dfa7cfd65ec95a61c224da1901036a4cca120a
EXPECTED_READS_SHA=3453e91714f7cc1b9bf9d18fbaac6dc82eed51e4f300644d7464bf08a0d00587
EXPECTED_BAM_SHA=d607aeac091f34ef4eca9e7cc9c091600de8bb857e99e1923713d2074b4cb949
EXPECTED_DEFAULT_REGIONS_SHA=b289b40748a98b2c1ed0e958f0b014e01017022fcda4f628e2444e32cf7b85df
EXPECTED_DEFAULT_THRESHOLDS_SHA=c483d2b698dd8cc3aa7b6f92233bf175ca3544d25b600f6491e196cf18812669
EXPECTED_DEFAULT_SUMMARY_SHA=7055b198de2ce23fbc2b955de0de04e61e80e4ecb865a83f41e41484b01839a3
EXPECTED_CONTIGS=7117
EXPECTED_SPAN=291889606

for tool in samtools mosdepth seqkit sha256sum python; do
    command -v "$tool" >/dev/null
done
test -s "$ASM"
test -s "$READS"
test -s "$BAM"
test -s "$BAM.bai"
test -s "$SOURCE/COMPLETED.txt"
test -s "$SOURCE/input_manifest.tsv"
test -s "$SOURCE/input_sha256.tsv"
gzip -t "$READS"
grep -Fqx $'assembly\t/work1/sponge/lan_tier2_decontam_review/assembly/lan_host_review_filtered.fasta' "$SOURCE/input_manifest.tsv"
grep -Fqx $'reads\t/public/home/sponge/sci_data/lan/00_input/Lanm.ont.fq.gz' "$SOURCE/input_manifest.tsv"
grep -Fqx "$EXPECTED_FASTA_SHA  $ASM" "$SOURCE/input_sha256.tsv"
grep -Fqx "$EXPECTED_READS_SHA  $READS" "$SOURCE/input_sha256.tsv"
grep -Fq 'minimap2 -ax map-ont -t 22 -I 20G' "$SOURCE/minimap2.log"
if [[ -e "$FINAL_OUT" ]]; then
    echo "Refusing to overwrite existing output: $FINAL_OUT" >&2
    exit 1
fi
mkdir "$OUT"

actual_fasta_sha=$(sha256sum "$ASM" | awk '{print $1}')
[[ "$actual_fasta_sha" == "$EXPECTED_FASTA_SHA" ]]
printf '%s  %s\n' "$actual_fasta_sha" "$ASM" > "$OUT/assembly.sha256"

# Revalidate that the extant 7117-contig BAM is tied to the same raw reads.
actual_reads_sha=$(sha256sum "$READS" | awk '{print $1}')
[[ "$actual_reads_sha" == "$EXPECTED_READS_SHA" ]]
printf '%s  %s\n' "$actual_reads_sha" "$READS" > "$OUT/reads.sha256"
actual_bam_sha=$(sha256sum "$BAM" | awk '{print $1}')
[[ "$actual_bam_sha" == "$EXPECTED_BAM_SHA" ]]
printf '%s  %s\n' "$actual_bam_sha" "$BAM" > "$OUT/bam.sha256"
[[ "$(sha256sum "$SOURCE/lan.review_filtered.ONT.by_contig.regions.bed.gz" | awk '{print $1}')" == "$EXPECTED_DEFAULT_REGIONS_SHA" ]]
[[ "$(sha256sum "$SOURCE/lan.review_filtered.ONT.by_contig.thresholds.bed.gz" | awk '{print $1}')" == "$EXPECTED_DEFAULT_THRESHOLDS_SHA" ]]
[[ "$(sha256sum "$SOURCE/lan.review_filtered.ONT.by_contig.mosdepth.summary.txt" | awk '{print $1}')" == "$EXPECTED_DEFAULT_SUMMARY_SHA" ]]
seqkit stats -T -j 8 "$READS" > "$OUT/reads.seqkit.tsv"

samtools quickcheck -v "$BAM"
samtools flagstat -@ 4 "$BAM" > "$OUT/NJ-lan.final7117.ONT.flagstat.txt"
primary_records=$(samtools view -@ 4 -c -F 2816 "$BAM")
primary_mapped=$(samtools view -@ 4 -c -F 2820 "$BAM")
primary_mapped_q20=$(samtools view -@ 4 -c -F 2820 -q 20 "$BAM")
printf 'metric\tcount\nprimary_records_F2816\t%s\nprimary_mapped_F2820\t%s\nprimary_mapped_q20_F2820\t%s\n' \
    "$primary_records" "$primary_mapped" "$primary_mapped_q20" \
    > "$OUT/mapping_counts.tsv"
samtools idxstats "$BAM" > "$OUT/NJ-lan.final7117.ONT.idxstats.tsv"
awk 'BEGIN{OFS="\t"} $1!="*"{print $1,0,$2}' \
    "$OUT/NJ-lan.final7117.ONT.idxstats.tsv" > "$BED"
[[ "$(wc -l < "$BED")" -eq "$EXPECTED_CONTIGS" ]]
awk -F '\t' -v expected="$EXPECTED_SPAN" \
    '{s+=$3-$2} END{if(s!=expected) exit 1}' "$BED"

cp "$SOURCE/lan.review_filtered.ONT.by_contig.regions.bed.gz" \
    "$OUT/NJ-lan.final7117.ONT.default.regions.bed.gz"
cp "$SOURCE/lan.review_filtered.ONT.by_contig.thresholds.bed.gz" \
    "$OUT/NJ-lan.final7117.ONT.default.thresholds.bed.gz"
cp "$SOURCE/lan.review_filtered.ONT.by_contig.mosdepth.summary.txt" \
    "$OUT/NJ-lan.final7117.ONT.default.mosdepth.summary.txt"
cp "$SOURCE/minimap2.version.txt" "$OUT/minimap2.version.txt"
cp "$SOURCE/samtools.version.txt" "$OUT/samtools.version.txt"
cp "$SOURCE/mosdepth.version.txt" "$OUT/mosdepth.version.txt"
cp "$SOURCE/minimap2.log" "$OUT/minimap2.log"
seqkit version > "$OUT/seqkit.version.txt" 2>&1
python --version > "$OUT/python.version.txt" 2>&1

printf 'sample_id\t%s\nassembly\t%s\nreads\t%s\nbam\t%s\nsource_mapping_job\t22273558\nminimap2\t-ax map-ont -t 22 -I 20G\ncoverage_default\tmosdepth -n -F 1796 -Q 0 --by whole-contig BED --thresholds 1,5,10,20 (copied from completed source mapping)\ncoverage_primary\tmosdepth -n -F 3844 -Q 0 --by whole-contig BED --thresholds 1,5,10,20\ncoverage_primary_q20\tmosdepth -n -F 3844 -Q 20 --by whole-contig BED --thresholds 1,5,10,20\n' \
    "$SAMPLE" "$ASM" "$READS" "$BAM" > "$OUT/input_manifest.tsv"

mosdepth -n -t 4 -F 3844 --by "$BED" --thresholds 1,5,10,20 \
    "$PRIMARY_PREFIX" "$BAM"

mosdepth -n -t 4 -F 3844 -Q 20 --by "$BED" --thresholds 1,5,10,20 \
    "$PRIMARY_Q20_PREFIX" "$BAM"

for path in \
    "$OUT/NJ-lan.final7117.ONT.default.regions.bed.gz" \
    "$OUT/NJ-lan.final7117.ONT.default.thresholds.bed.gz" \
    "$PRIMARY_PREFIX.regions.bed.gz" \
    "$PRIMARY_PREFIX.thresholds.bed.gz" \
    "$PRIMARY_Q20_PREFIX.regions.bed.gz" \
    "$PRIMARY_Q20_PREFIX.thresholds.bed.gz"; do
    test -s "$path"
done

stat -c 'bam_path\t%n\nbam_size_bytes\t%s\nbam_mtime_epoch\t%Y' "$BAM" \
    > "$OUT/bam.stat.tsv"
date -Is > "$OUT/COMPLETED.txt"
mv -T "$OUT" "$FINAL_OUT"
echo "NJ-lan final 7117-contig primary-only coverage audit completed."
