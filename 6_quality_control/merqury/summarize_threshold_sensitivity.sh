#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 LHT.spectra-cn.hist NJ-lan.spectra-cn.hist" >&2
  exit 2
fi

summarize() {
  local label=$1
  local hist=$2
  awk -F '\t' -v label="$label" '
    NR > 1 {
      for (threshold = 1; threshold <= 10; threshold++) {
        if ($2 > threshold) {
          total[threshold] += $3
          if ($1 != "read-only") represented[threshold] += $3
        }
      }
    }
    END {
      for (threshold = 1; threshold <= 10; threshold++) {
        printf "%s_gt%d\t%.0f\t%.0f\t%.4f\n", label, threshold,
               represented[threshold], total[threshold],
               100 * represented[threshold] / total[threshold]
      }
    }
  ' "$hist"
}

printf 'sample_threshold\trepresented_distinct_19mers\treliable_distinct_19mers\tcompleteness_percent\n'
summarize LHT "$1"
summarize NJlan "$2"
