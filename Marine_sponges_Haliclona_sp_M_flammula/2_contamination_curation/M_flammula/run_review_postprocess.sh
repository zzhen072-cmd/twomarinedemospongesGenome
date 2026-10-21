#!/bin/bash
#SBATCH --job-name=lan_review_plot
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=lan_review_plot.%j.out
#SBATCH --error=lan_review_plot.%j.err

set -euo pipefail
export LC_ALL=C

module load apps/miniforge3/25.3.1-0
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gzx_qc

ROOT=/work1/sponge/lan_tier2_decontam_review
PLOT="$ROOT/plot"
ASM="$ROOT/assembly/lan_host_review_filtered.fasta"
RAW_BLAST=/work1/sponge/lan_blobplot_final_coverage/source_data/lan_cleaned.nt.blastn.tsv
TAXONOMY=/work1/sponge/lan_blobplot_final_coverage/taxonomy/taxonomy_records.json
COVERAGE="$ROOT/mapping/lan_review_filtered_classification_coverage.tsv"
RETAINED_IDS="$ROOT/audit/retained_derivative_assembly.ids.txt"
FILTERED_BLAST="$ROOT/source_data/lan_review_filtered.nt.blastn.tsv"
PREFIX="$PLOT/lan_review_filtered_blobplot"

test -s "$ROOT/mapping/COMPLETED.txt"
test -s "$ASM"
test -s "$RAW_BLAST"
test -s "$TAXONOMY"
test -s "$COVERAGE"
test -s "$RETAINED_IDS"
if [[ -e "$ROOT/PLOT_COMPLETED.txt" ]]; then
    echo "Completed plot output already exists: $ROOT/PLOT_COMPLETED.txt" >&2
    exit 1
fi
mkdir -p "$PLOT" "$ROOT/source_data" "$ROOT/mplconfig"
export MPLCONFIGDIR="$ROOT/mplconfig"

python "$ROOT/scripts/filter_blast_queries.py" \
    --input "$RAW_BLAST" \
    --retained-ids "$RETAINED_IDS" \
    --output "$FILTERED_BLAST" \
    --audit "$ROOT/source_data/blast_filter_audit.json"

python "$ROOT/scripts/build_lan_blobplot_data.py" \
    --blast "$FILTERED_BLAST" \
    --fasta "$ASM" \
    --coverage-table "$COVERAGE" \
    --taxonomy-json "$TAXONOMY" \
    --outdir "$PLOT" \
    --blast-complete

python "$ROOT/scripts/plot_lan_blob.py" \
    --table "$PLOT/lan_final_plot_data.tsv" \
    --out-prefix "$PREFIX" \
    --title "LAN evidence-reviewed host assembly: GC–coverage plot" \
    --subtitle "Confirmed prokaryotic contaminants removed as whole contigs; ONT reads remapped to the filtered derivative assembly" \
    --coverage-label "ONT coverage on filtered host assembly (×)" \
    --note "The original assembly remains unchanged; excluded contigs are archived in a separate FASTA."

tar -czf "$ROOT/lan_review_filtered_results.tar.gz" \
    -C "$ROOT" \
    audit/lan_final_decontam_review.tsv \
    audit/decontam_review_audit.json \
    audit/host_gene_guard.tsv \
    audit/host_gene_guard_audit.json \
    audit/lan_final_decontam_summary.tsv \
    audit/remove_confirmed_contaminants.ids.txt \
    audit/removed_confirmed_contaminants.tsv \
    audit/tier2_removed_confirmed_contaminants.ids.txt \
    audit/tier2_removed_confirmed_contaminants.tsv \
    audit/tier2_retained_host_supported.ids.txt \
    audit/tier2_retained_host_supported.tsv \
    audit/tier2_retained_ambiguous_review.ids.txt \
    audit/tier2_retained_ambiguous_review.tsv \
    assembly/filter_audit.json \
    assembly/fasta_sha256.tsv \
    mapping/coverage_audit.json \
    mapping/input_manifest.tsv \
    mapping/input_sha256.tsv \
    mapping/lan.review_filtered.ONT.flagstat.txt \
    mapping/minimap2.log \
    mapping/minimap2.version.txt \
    mapping/mosdepth.version.txt \
    mapping/samtools.version.txt \
    mapping/lan_review_filtered_classification_coverage.tsv \
    source_data/blast_filter_audit.json \
    source_data/lan_review_filtered.nt.blastn.tsv \
    plot/analysis_audit.json \
    plot/duplicate_blast_rows.tsv \
    plot/lan_blast_hits_annotated.tsv \
    plot/lan_final_phylum_summary.tsv \
    plot/lan_final_plot_data.tsv \
    plot/lan_review_filtered_blobplot.png \
    plot/lan_review_filtered_blobplot.svg \
    plot/lan_review_filtered_blobplot.pdf \
    plot/lan_review_filtered_blobplot.alignment.json \
    plot/lan_review_filtered_blobplot.alignment-overlay.svg \
    plot/lan_review_filtered_blobplot.category_summary.tsv \
    plot/lan_review_filtered_blobplot.histogram_bins.tsv \
    plot/lan_review_filtered_blobplot.plotted_points.tsv \
    scripts/audit_panel_alignment.py \
    scripts/build_lan_blobplot_data.py \
    scripts/filter_blast_queries.py \
    scripts/filter_fasta_exact.py \
    scripts/make_lan_coverage_table.py \
    scripts/plot_lan_blob.py

sha256sum "$ROOT/lan_review_filtered_results.tar.gz" \
    > "$ROOT/lan_review_filtered_results.tar.gz.sha256"
date -Is > "$ROOT/PLOT_COMPLETED.txt"
