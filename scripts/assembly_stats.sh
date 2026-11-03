#!/usr/bin/env bash
# =============================================================================
#  assembly_stats.sh
#  Statistics for one or more assemblies (FASTA or FASTA.gz) in a single table,
#  for easy comparison.
#
#  Usage:
#    bash assembly_stats.sh [-m MINLEN] [-t] assembly1.fasta [assembly2.fasta.gz ...]
#
#    -m MINLEN  Count only contigs >= MINLEN bp (default: all)
#    -t         Output as TSV (for Excel/R) instead of an aligned table
#
#  Example:
#    bash assembly_stats.sh -m 1000 \
#      soil_assembly/05_metaspades/contigs.fasta \
#      megahit_out/final.contigs.fa \
#      | tee comparison.txt
# =============================================================================
set -euo pipefail

MINLEN=0
TSV=0

usage() { sed -n '6,12p' "$0" | sed 's/^# \{0,1\}//'; }

while getopts "m:th" opt; do
  case "$opt" in
    m) MINLEN="$OPTARG" ;;
    t) TSV=1 ;;
    h) usage; exit 0 ;;
    *) usage; exit 1 ;;
  esac
done
shift $(( OPTIND - 1 ))
[[ $# -ge 1 ]] || { usage; exit 1; }

stats() {
  printf "Assembly\tContigs\tTotal_Mbp\tLargest_bp\tN50_bp\tL50\tN90_bp\t>=1kb\t>=5kb\t>=10kb\tMbp_in_>=1kb\tGC_%%\n"

  for f in "$@"; do
    if [[ ! -s "$f" ]]; then
      echo "Skipping: file not found or empty: $f" >&2
      continue
    fi
    # Label: directory/file, because many assemblies are simply named contigs.fasta
    local label reader="cat"
    label="$(basename "$(dirname "$(realpath "$f")")")/$(basename "$f")"
    [[ "$f" == *.gz ]] && reader="zcat"

    "$reader" "$f" \
      | awk -v min="$MINLEN" '
          /^>/ { if (seen && len >= min) print len, gc; len = 0; gc = 0; seen = 1; next }
          { len += length($0); gc += gsub(/[GCgcSs]/, "") }
          END  { if (seen && len >= min) print len, gc }' \
      | sort -k1,1nr \
      | awk -v label="$label" '
          { l[NR] = $1; tot += $1; gc += $2
            if ($1 >= 1000)  { k1++; tot1k += $1 }
            if ($1 >= 5000)  k5++
            if ($1 >= 10000) k10++ }
          END {
            if (NR == 0) { printf "%s\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tNA\n", label; exit }
            run = 0
            for (i = 1; i <= NR; i++) {
              run += l[i]
              if (!n50 && run >= tot * 0.5) { n50 = l[i]; l50 = i }
              if (!n90 && run >= tot * 0.9) { n90 = l[i]; break }
            }
            printf "%s\t%d\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.1f\n",
                   label, NR, tot / 1e6, l[1], n50, l50, n90,
                   k1, k5, k10, tot1k / 1e6, 100 * gc / tot
          }'
  done
}

if (( TSV )); then
  stats "$@"
else
  (( MINLEN > 0 )) && echo "Only contigs >= ${MINLEN} bp"
  if command -v column >/dev/null 2>&1; then
    stats "$@" | column -t -s $'\t'
  else
    stats "$@"
  fi
fi
