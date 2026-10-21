#!/bin/bash
#SBATCH --job-name=finalann_braker
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=7-00:00:00
#SBATCH --array=0-1
#SBATCH --output=/work1/sponge/lht-njlan-final-annotation-20261007/logs/braker_master_%A_%a.out
#SBATCH --error=/work1/sponge/lht-njlan-final-annotation-20261007/logs/braker_master_%A_%a.err
set -eo pipefail

ROOT=/work1/sponge/lht-njlan-final-annotation-20261007
SAMPLES=(lht04 njlan)
OLD_AUG=(
  /work1/sponge/genome_annotation/lht_flye118_manual/braker4/augustus_config
  /work1/sponge/genome_annotation/lan_flye118_manual/braker4/augustus_config
)
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
RUN="$ROOT/$SAMPLE/03_braker4"
SOURCE_AUG=${OLD_AUG[$SLURM_ARRAY_TASK_ID]}
MASKED="$ROOT/$SAMPLE/02_repeatmasker/${SAMPLE}_final.shortids.fasta.masked"
RESULT="$RUN/output/${SAMPLE}_final/results/braker.aa.gz"

test -s "$MASKED"
test -s "$RUN/samples.csv"
test -s "$RUN/config.ini"
test -s /work1/sponge/genome_annotation/braker4_containers/braker3_v3.1.1.sif
test -s /work1/sponge/genome_annotation/braker4_containers/agat_1.4.1.sif
test -s /work1/sponge/genome_annotation/protein_evidence/Porifera_2species.faa
test ! -e "$RESULT"
if [ ! -d "$RUN/augustus_config" ]; then
  cp -a "$SOURCE_AUG" "$RUN/augustus_config"
fi

source "$HOME/.bashrc"
conda activate /work1/sponge/conda_envs/braker4
module load singularity/3.8.0
cd "$RUN"
export BRAKER4_CONFIG="$RUN/config.ini"
snakemake -s /work1/sponge/software/BRAKER4/Snakefile -n --cores 30 \
  --use-singularity \
  --singularity-prefix /work1/sponge/genome_annotation/braker4_singularity \
  > "$RUN/dryrun.log"
snakemake -s /work1/sponge/software/BRAKER4/Snakefile \
  --executor slurm --jobs 15 --cores 30 \
  --default-resources slurm_partition=normal \
  --use-singularity \
  --singularity-prefix /work1/sponge/genome_annotation/braker4_singularity \
  --singularity-args '-B /public,/work1' \
  --latency-wait 120 --restart-times 2 --rerun-incomplete --keep-going --printshellcmds
test -s "$RESULT"
