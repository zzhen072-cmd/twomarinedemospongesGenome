#!/bin/bash
#SBATCH -J lan_blobtoolkit
#SBATCH -p normal
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 30
#SBATCH --mem=108G
#SBATCH -o /work1/sponge/lan_reassembly_compare/12_blobtoolkit/blobtoolkit-%j.out
#SBATCH -e /work1/sponge/lan_reassembly_compare/12_blobtoolkit/blobtoolkit-%j.err
set -eo pipefail

ASSEMBLY=/work1/sponge/lan_reassembly_compare/10_fcs_gx/lan_fcs_adaptor.fasta
ONT=/public/home/sponge/sci_data/lan/00_input/Lanm.ont.fq.gz
KRAKEN=/work1/sponge/lan_reassembly_compare/11_kraken2/lan.kraken2.out
OUT=/work1/sponge/lan_reassembly_compare/12_blobtoolkit
BAM=$OUT/lan.ont_to_fcs_adaptor.bam
FIELDS=$OUT/lan.kraken2_fields.tsv
BLOBDIR=$OUT/lan_fcs_adaptor.blobdir
BLOB_ENV=/public/home/sponge/conda_envs/blobtoolkit_env
source ~/.bashrc
conda activate gzx_qc
set -u
test -s "$ASSEMBLY"
test -s "$ONT"
test -s "$KRAKEN"
test -x "$BLOB_ENV/bin/blobtools"
test ! -e "$BLOBDIR"
echo "START $(date)"
TMP_PREFIX=${TMPDIR:-/tmp}/lan_blob_${SLURM_JOB_ID}
minimap2 -ax map-ont -t 26 "$ASSEMBLY" "$ONT" | \
  samtools sort -@ 4 -m 4G -T "$TMP_PREFIX" -o "$BAM" -
samtools index -@ 4 -c "$BAM"
samtools quickcheck -v "$BAM"
{ printf 'scaffold_identifier\tkraken2_status\tkraken2_taxid\tkraken2_query_length\n'; awk -F '\t' 'BEGIN{OFS="\t"}{print $2,$1,$3,$4}' "$KRAKEN"; } > "$FIELDS"
conda run -p "$BLOB_ENV" blobtools add --create --fasta "$ASSEMBLY" --cov "$BAM" \
  --text "$FIELDS" --text-delimiter "\t" --text-header --text-no-array \
  --threads 30 "$BLOBDIR"
test -d "$BLOBDIR"
find "$BLOBDIR" -maxdepth 1 -type f -print | sort > "$OUT/lan.blobdir_files.txt"
echo "FINISH $(date)"
