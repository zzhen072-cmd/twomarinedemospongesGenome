#!/bin/bash
#SBATCH --job-name=lhtlan_fig_repeat
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=/work1/sponge/lht-njlan-final-figures-20261009/logs/repeat_landscape_%j.out
#SBATCH --error=/work1/sponge/lht-njlan-final-figures-20261009/logs/repeat_landscape_%j.err
set -euo pipefail

export LC_ALL=C
export PATH=/work1/sponge/conda_envs/repeatmodeler_ltr/bin:$PATH
RM_ROOT=/work1/sponge/conda_envs/repeatmodeler_ltr/share/RepeatMasker
SOURCE_ROOT=/work1/sponge/lht-njlan-final-annotation-20261007
DEST=/work1/sponge/lht-njlan-final-figures-20261009/repeat_landscape

process_sample() {
    local sample=$1
    local stem=$2
    local genome_bp=$3
    local library=$4
    local source="$SOURCE_ROOT/$sample/02_repeatmasker"
    local work="$DEST/$sample"

    mkdir -p "$work"
    cd "$work"
    test ! -e DONE
    sha256sum "$source/$stem.cat.gz" "$source/$stem.out" "$source/$stem.tbl" "$library" > input_sha256.before.txt
    cp -p "$source/$stem.cat.gz" .
    ProcessRepeats -a -lib "$library" "$stem.cat.gz" > process_repeats.log 2>&1
    test -s "$stem.align"
    perl "$RM_ROOT/util/calcDivergenceFromAlign.pl" -s "$sample.divsum" "$stem.align" > calc_divergence.log 2>&1
    test -s "$sample.divsum"
    sha256sum -c input_sha256.before.txt > input_unchanged.txt
    sha256sum "$sample.divsum" "$stem.align" > output_sha256.txt
    printf '%s\t%s\n' "$sample" "$genome_bp" > genome_span.tsv
    date -Is > DONE
}

process_sample lht04 lht04_final.shortids.fasta 115734075 /work1/sponge/genome_annotation/lht_flye118_manual/repeatmodeler/lht_db-families.fa &
LHT_PID=$!
process_sample njlan njlan_final.shortids.fasta 291889606 /work1/sponge/genome_annotation/lan_flye118_manual/repeatmodeler/lan_db-families.fa &
NJ_PID=$!
wait "$LHT_PID"
wait "$NJ_PID"
date -Is > "$DEST/DONE"
