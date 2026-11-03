#!/usr/bin/env bash
# =============================================================================
#  run_all.sh
#  Runs soil_assembly.sh for every Illumina sample, one after the other.
#  Every 3 consecutive Illumina samples are paired with 1 PacBio sample:
#      Met_7  -> Met_19, Met_20, Met_21
#      Met_8  -> Met_22, Met_23, Met_24
#      ...
#      Met_14 -> Met_40, Met_41, Met_42
#
#  Usage:
#    bash run_all.sh -n     # dry run: shows the pairs and checks the input files
#    nohup bash run_all.sh > run_all.out 2>&1 &
#
#  If it stops (e.g. reboot), run it again: completed samples are skipped and
#  the unfinished one resumes from its last checkpoint.
# =============================================================================
set -uo pipefail

# ---------------------------- Settings --------------------------------------
PIPELINE="$(dirname "$(realpath "$0")")/soil_assembly.sh"

ILLUMINA_DIR="/path/to/illumina"     # contains subdirectories Met_19/, Met_20/, ...
PACBIO_DIR="/path/to/pacbio"         # contains Met_7.hifi_reads.fastq, ...
OUT_BASE="/path/to/output/batch"

FIRST_PACBIO=7        # Met_7
LAST_PACBIO=14        # Met_14
FIRST_ILLUMINA=19     # Met_19
PER_GROUP=3           # Illumina samples per PacBio sample

# Samples to skip, e.g. SKIP_SAMPLES=(Met_19 Met_25)
SKIP_SAMPLES=()

THREADS=32
MEM_GB=220
NORM_TARGET=50

# 1 = after a successful assembly, delete the large intermediate files
#     (tadpole/bbnorm reads, SPAdes K* and tmp directories).
#     Contigs, scaffolds, the assembly graph and all logs are kept.
CLEANUP=1

# 1 = before starting, wait for any metaSPAdes / pipeline run already in progress
WAIT_FOR_RUNNING=1

# File names (edit to match your data)
r1_path()     { echo "$ILLUMINA_DIR/$1/${1}_R1_repaired.fq.gz"; }
r2_path()     { echo "$ILLUMINA_DIR/$1/${1}_R2_repaired.fq.gz"; }
pacbio_path() { echo "$PACBIO_DIR/${1}.hifi_reads.fastq"; }

# ---------------------------- Helpers ---------------------------------------
DRY_RUN=0
[[ "${1:-}" == "-n" ]] && DRY_RUN=1

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

in_skip_list() {
  local s
  for s in ${SKIP_SAMPLES[@]+"${SKIP_SAMPLES[@]}"}; do [[ "$s" == "$1" ]] && return 0; done
  return 1
}

cleanup_sample() {
  local out="$1"
  rm -f  "$out"/02_tadpole/*.fq.gz "$out"/03_bbnorm/*.fq.gz
  rm -rf "$out"/05_metaspades/K* "$out"/05_metaspades/tmp "$out"/tmp
}

# ---------------------------- Sample pairs ----------------------------------
PAIRS=()
for (( pb = FIRST_PACBIO; pb <= LAST_PACBIO; pb++ )); do
  for (( j = 0; j < PER_GROUP; j++ )); do
    il=$(( FIRST_ILLUMINA + (pb - FIRST_PACBIO) * PER_GROUP + j ))
    PAIRS+=("Met_${il}:Met_${pb}")
  done
done

log "Illumina -> PacBio pairs (${#PAIRS[@]} in total):"
missing=0
for pair in "${PAIRS[@]}"; do
  il="${pair%%:*}"; pb="${pair##*:}"
  status="OK"
  for f in "$(r1_path "$il")" "$(r2_path "$il")" "$(pacbio_path "$pb")"; do
    if [[ ! -s "$f" ]]; then status="MISSING: $f"; missing=1; fi
  done
  in_skip_list "$il" && status="skipped (SKIP_SAMPLES)"
  printf "    %-7s + %-7s  %s\n" "$il" "$pb" "$status"
done

if (( missing )); then
  log "Some files are missing. Are all disks mounted? Fix and run again."
  exit 1
fi
(( DRY_RUN )) && { log "Dry run: all files found. Nothing was executed."; exit 0; }

[[ -f "$PIPELINE" ]] || { log "$PIPELINE not found"; exit 1; }
command -v metaspades.py >/dev/null 2>&1 \
  || { log "metaspades.py not found. Is the conda environment activated?"; exit 1; }

# ---------------------------- Wait for a run already in progress -----------
if (( WAIT_FOR_RUNNING )); then
  while pgrep -u "$USER" -f "spades-core" >/dev/null \
     || pgrep -u "$USER" -f "soil_assembly\.sh" >/dev/null; do
    log "Another assembly is running. Waiting 10 minutes..."
    sleep 600
  done
fi

# ---------------------------- Run -------------------------------------------
mkdir -p "$OUT_BASE"
SUMMARY="$OUT_BASE/batch_summary.tsv"
[[ -f "$SUMMARY" ]] || printf "illumina\tpacbio\tstatus\tstart\tend\n" > "$SUMMARY"

for pair in "${PAIRS[@]}"; do
  il="${pair%%:*}"; pb="${pair##*:}"
  out="$OUT_BASE/$il"

  if in_skip_list "$il"; then
    log "$il: skipped (SKIP_SAMPLES)"
    continue
  fi
  if [[ -f "$out/.checkpoints/metaspades.done" ]]; then
    log "$il: already completed, skipping"
    continue
  fi

  start=$(date '+%Y-%m-%d %H:%M')
  log "$il + $pb: starting (log: $out.run.out)"
  mkdir -p "$out"

  if bash "$PIPELINE" -S \
       -1 "$(r1_path "$il")" -2 "$(r2_path "$il")" \
       -l "$(pacbio_path "$pb")" \
       -t "$THREADS" -m "$MEM_GB" -n "$NORM_TARGET" \
       -o "$out" > "$out.run.out" 2>&1; then
    status="OK"
    log "$il: completed"
    if (( CLEANUP )); then
      cleanup_sample "$out"
      log "$il: intermediate files deleted"
    fi
  else
    status="FAILED"
    log "$il: FAILED. See $out.run.out and $out/logs/. Continuing with the next sample."
  fi

  printf "%s\t%s\t%s\t%s\t%s\n" "$il" "$pb" "$status" "$start" "$(date '+%Y-%m-%d %H:%M')" >> "$SUMMARY"
done

log "Batch finished. Summary: $SUMMARY"
column -t -s $'\t' "$SUMMARY" 2>/dev/null || cat "$SUMMARY"
