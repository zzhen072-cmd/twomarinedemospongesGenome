#!/bin/bash
# Successful final InterProScan launcher. This is the Python 3.11 retry that
# replaced two failed launch attempts; those failed launchers are not included.
#SBATCH --job-name=finalann_ipr_py311
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=96G
#SBATCH --time=2-00:00:00
#SBATCH --array=0-1
#SBATCH --export=NONE
#SBATCH --output=/work1/sponge/lht-njlan-final-annotation-20261007/logs/interpro_retry_py311_%A_%a.out
#SBATCH --error=/work1/sponge/lht-njlan-final-annotation-20261007/logs/interpro_retry_py311_%A_%a.err
set -euo pipefail

ROOT=/work1/sponge/lht-njlan-final-annotation-20261007
SAMPLES=(lht04 njlan)
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
INPUT="$ROOT/$SAMPLE/04_protein_busco/${SAMPLE}_final.longest_isoform.aa"
OUT="$ROOT/$SAMPLE/07_interpro"
RESULT="$OUT/${SAMPLE}_final.longest.interproscan.tsv"
STAGED="$OUT/${SAMPLE}_final.longest.interproscan.part.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.tsv"
INTERPRO_ROOT=/work1/sponge/software/interproscan-5.77-108.0
INTERPRO="$INTERPRO_ROOT/interproscan.sh"
JAVA_ENV=/work1/sponge/conda_envs/interproscan_java11
PYTHON_ENV=/work1/sponge/conda_envs/eggnog_mapper_2.1.13
TMP="$OUT/tmp_retry_py311_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
ENV_AUDIT="$OUT/${SAMPLE}.interpro.retry_py311.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.environment.txt"
HASH_AUDIT="$OUT/${SAMPLE}.interpro.retry_py311.${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.sha256"

test -s "$INPUT"
test -x "$INTERPRO"
test -x "$JAVA_ENV/bin/java"
test -x "$PYTHON_ENV/bin/python3"
test ! -e "$RESULT"
test ! -e "$RESULT.gz"
test ! -e "$STAGED"
test ! -e "$STAGED.gz"
mkdir -p "$OUT" "$TMP"

unset PYTHONHOME PYTHONPATH CONDA_PREFIX CONDA_DEFAULT_ENV
export JAVA_HOME="$JAVA_ENV"
export PATH="$JAVA_ENV/bin:$PYTHON_ENV/bin:/usr/local/bin:/usr/bin:/bin"
hash -r

test "$(command -v python3)" = "$PYTHON_ENV/bin/python3"
python3 -c 'import sys; assert sys.version_info >= (3, 8)'
(
  cd "$INTERPRO_ROOT"
  PYTHONPATH=bin/panther python3 -c 'from Bio.Align import _aligners; print(_aligners.__file__)'
)
(
  date -Is
  hostname
  printf 'python3\t'
  python3 --version
  printf 'java\t'
  java -version 2>&1 | head -1
  "$INTERPRO" --version 2>&1 || true
) > "$ENV_AUDIT"

cd "$INTERPRO_ROOT"
"$INTERPRO" -i "$INPUT" -f TSV -o "$STAGED" \
  -cpu "$SLURM_CPUS_PER_TASK" -goterms -pa -dp \
  -exclappl MobiDBLite,ProSiteProfiles,ProSitePatterns -T "$TMP"

test -s "$STAGED"
gzip -n "$STAGED"
test -s "$STAGED.gz"
ln "$STAGED.gz" "$RESULT.gz"
test -s "$RESULT.gz"
sha256sum "$INPUT" "$STAGED.gz" "$RESULT.gz" > "$HASH_AUDIT"
printf 'PASS\t%s\t%s\n' "$SAMPLE" "$RESULT.gz"
