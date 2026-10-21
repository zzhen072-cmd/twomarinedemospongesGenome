#!/usr/bin/env bash
#SBATCH --job-name=lht_njlan_quast
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --output=/work1/sponge/lht_njlan_final_quast_20261014/slurm-%j.out
#SBATCH --error=/work1/sponge/lht_njlan_final_quast_20261014/slurm-%j.err

set -euo pipefail

run_root=/work1/sponge/lht_njlan_final_quast_20261014
quast_env=/work1/sponge/lht_njlan_quast_5.3.0_env
quast_bin="${quast_env}/bin/quast.py"
lht_fasta=/work1/sponge/lht_reassembly_compare/15_read_validated_curated/LHT_curated_primary.fasta
nj_fasta=/work1/sponge/lan_tier2_decontam_review/assembly/lan_host_review_filtered.fasta

mkdir -p "${run_root}/results/LHT-04" "${run_root}/results/NJ-lan"

"${quast_bin}" --version > "${run_root}/quast_version.txt" 2>&1
sha256sum "${lht_fasta}" "${nj_fasta}" > "${run_root}/input_fasta.sha256"

printf 'Sample_ID\tInput_FASTA\tInput_SHA256\tQUAST_Version\tCommand\n' > "${run_root}/run_manifest.tsv"
printf 'LHT-04\t%s\t%s\tQUAST v5.3.0\t%s\n' \
  "${lht_fasta}" \
  '026d8b2e3e9888de0763d42ce226fcf7c63062a145e6ad8478e2e8c1725b2499' \
  'quast.py --min-contig 0 --threads 16 --report-all-metrics -o results/LHT-04 <locked-final-FASTA>' \
  >> "${run_root}/run_manifest.tsv"
printf 'NJ-lan\t%s\t%s\tQUAST v5.3.0\t%s\n' \
  "${nj_fasta}" \
  'cea9b34517164f02103111a3a1dfa7cfd65ec95a61c224da1901036a4cca120a' \
  'quast.py --min-contig 0 --threads 16 --report-all-metrics -o results/NJ-lan <locked-final-FASTA>' \
  >> "${run_root}/run_manifest.tsv"

"${quast_bin}" \
  --min-contig 0 \
  --threads 16 \
  --report-all-metrics \
  -o "${run_root}/results/LHT-04" \
  "${lht_fasta}"

"${quast_bin}" \
  --min-contig 0 \
  --threads 16 \
  --report-all-metrics \
  -o "${run_root}/results/NJ-lan" \
  "${nj_fasta}"

date -u +'%Y-%m-%dT%H:%M:%SZ' > "${run_root}/completed_utc.txt"
sha256sum "${run_root}/results/LHT-04/report.tsv" "${run_root}/results/NJ-lan/report.tsv" > "${run_root}/report_tsv.sha256"
