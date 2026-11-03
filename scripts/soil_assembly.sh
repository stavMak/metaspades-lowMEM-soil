#!/usr/bin/env bash
# =============================================================================
#  soil_assembly.sh
#  Metagenomic assembly of soil samples (short reads + PacBio long reads)
#
#  Steps:
#    1. fastp       -> trimming, adapter removal, QC report (skipped with -S)
#    2. tadpole.sh  -> error correction (BBTools)
#    3. bbnorm.sh   -> digital normalization (BBTools)
#    4. metaMDBG    -> HiFi-only assembly (ONLY with -L hifi)
#    5. metaSPAdes  -> hybrid assembly (--only-assembler) with the long reads
#    6. Assembly statistics
#
#  Every completed step is checkpointed. If the script stops, run it again with
#  the same parameters and it resumes where it left off.
# =============================================================================
set -euo pipefail

# ---------------------------- Defaults --------------------------------------
R1=""
R2=""
LONG=""
LONG_TYPE="clr"          # clr | hifi | nanopore
OUTDIR="soil_assembly"
THREADS=32
MEM_GB=220
TMP_DIR=""
KMERS="21,33,55,77"
NORM_TARGET=50
NORM_MIN=2
SKIP_FASTP=0

usage() {
cat <<USAGE
Usage:
  bash $(basename "$0") -1 R1.fastq.gz -2 R2.fastq.gz -l pacbio.fastq.gz [options]

Required:
  -1 FILE   Forward short reads (R1)
  -2 FILE   Reverse short reads (R2)

Optional:
  -l FILE   Long reads (FASTQ or FASTQ.gz)
  -L TYPE   Long-read type: clr | hifi | nanopore     (default: ${LONG_TYPE})
              clr      -> metaSPAdes hybrid (--pacbio)
              hifi     -> metaMDBG HiFi-only assembly + metaSPAdes hybrid (--pacbio)
              nanopore -> metaSPAdes hybrid (--nanopore)
  -o DIR    Output directory                      (default: ${OUTDIR})
  -t INT    Threads                               (default: ${THREADS})
  -m INT    Memory limit for metaSPAdes, in GB    (default: ${MEM_GB})
  -T DIR    Temporary directory on a fast disk    (default: OUTDIR/tmp)
  -k LIST   k-mer sizes for metaSPAdes            (default: ${KMERS})
  -n INT    bbnorm target depth                   (default: ${NORM_TARGET})
  -S        Skip fastp: reads are already trimmed
  -h        Show this help
USAGE
}

while getopts "1:2:l:L:o:t:m:T:k:n:Sh" opt; do
  case "$opt" in
    1) R1="$OPTARG" ;;
    2) R2="$OPTARG" ;;
    l) LONG="$OPTARG" ;;
    L) LONG_TYPE="${OPTARG,,}" ;;
    o) OUTDIR="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    m) MEM_GB="$OPTARG" ;;
    T) TMP_DIR="$OPTARG" ;;
    k) KMERS="$OPTARG" ;;
    n) NORM_TARGET="$OPTARG" ;;
    S) SKIP_FASTP=1 ;;
    h) usage; exit 0 ;;
    *) usage; exit 1 ;;
  esac
done

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
die() { echo "[ERROR] $*" >&2; exit 1; }

# ---------------------------- Checks ----------------------------------------
[[ -n "$R1" && -n "$R2" ]] || { usage; die "Both -1 and -2 are required."; }

case "$LONG_TYPE" in
  clr|hifi) SPADES_LONG_OPT="--pacbio" ;;
  nanopore) SPADES_LONG_OPT="--nanopore" ;;
  *) die "Unknown long-read type '$LONG_TYPE'. Choices: clr, hifi, nanopore" ;;
esac

for f in "$R1" "$R2" ${LONG:+"$LONG"}; do
  [[ -s "$f" ]] || die "File not found or empty: $f"
done

R1=$(realpath "$R1")
R2=$(realpath "$R2")
[[ -n "$LONG" ]] && LONG=$(realpath "$LONG")

RUN_METAMDBG=0
[[ -n "$LONG" && "$LONG_TYPE" == "hifi" ]] && RUN_METAMDBG=1

TOOLS=(tadpole.sh bbnorm.sh metaspades.py)
(( SKIP_FASTP ))   || TOOLS+=(fastp)
(( RUN_METAMDBG )) && TOOLS+=(metaMDBG)
for tool in "${TOOLS[@]}"; do
  command -v "$tool" >/dev/null 2>&1 \
    || die "'$tool' not found. Is the conda environment activated?"
done

# ---------------------------- Directories & logs ----------------------------
mkdir -p "$OUTDIR"
OUTDIR=$(realpath "$OUTDIR")
TMP_DIR=${TMP_DIR:-"$OUTDIR/tmp"}
mkdir -p "$TMP_DIR"

LOGDIR="$OUTDIR/logs"
CKDIR="$OUTDIR/.checkpoints"
QC="$OUTDIR/01_fastp"
EC="$OUTDIR/02_tadpole"
NORM="$OUTDIR/03_bbnorm"
MDBG="$OUTDIR/04_metamdbg"
ASM="$OUTDIR/05_metaspades"
mkdir -p "$LOGDIR" "$CKDIR" "$QC" "$EC" "$NORM"

MAIN_LOG="$LOGDIR/pipeline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$MAIN_LOG") 2>&1

trap 'log "The pipeline stopped with an error. See the logs in $LOGDIR"' ERR

# Memory for BBTools (Java heap): 70% of the limit. In practice Java uses
# 15-25 GB more than its heap, so this leaves comfortable headroom for the system.
BB_MEM=$(( MEM_GB * 70 / 100 ))
# fastp does not scale beyond 16 threads
FASTP_THREADS=$(( THREADS < 16 ? THREADS : 16 ))

# Record run time / peak memory if GNU time is available
USE_TIME=0
if [[ -x /usr/bin/time ]] && /usr/bin/time -v true >/dev/null 2>&1; then
  USE_TIME=1
fi

# run_step NAME command...
# Runs one step, writes its output to logs/NAME.log and marks it as completed.
run_step() {
  local name="$1"; shift
  if [[ -f "$CKDIR/$name.done" ]]; then
    log "Skipping '$name' (already completed)"
    return 0
  fi
  log "Starting '$name' (log: $LOGDIR/$name.log)"
  local start=$SECONDS
  if (( USE_TIME )); then
    /usr/bin/time -v -o "$LOGDIR/$name.time.txt" "$@" > "$LOGDIR/$name.log" 2>&1
  else
    "$@" > "$LOGDIR/$name.log" 2>&1
  fi
  local mins=$(( (SECONDS - start) / 60 ))
  local peak=""
  if (( USE_TIME )); then
    peak=$(awk -F': ' '/Maximum resident set size/ {printf "%.1f GB", $2/1024/1024}' \
           "$LOGDIR/$name.time.txt")
  fi
  touch "$CKDIR/$name.done"
  log "Finished '$name' in ${mins} min${peak:+, peak RAM: $peak}"
}

# assembly_stats LABEL FILE(.fasta or .fasta.gz)
assembly_stats() {
  local label="$1" fasta="$2"
  if [[ ! -s "$fasta" ]]; then
    log "$fasta not found, skipping statistics for $label"
    return 0
  fi
  local reader="cat"
  [[ "$fasta" == *.gz ]] && reader="zcat"
  echo "  $label  ($fasta)"
  "$reader" "$fasta" \
    | awk '/^>/ { if (len) print len; len = 0; next } { len += length($0) } END { if (len) print len }' \
    | sort -rn \
    | awk '{ l[NR] = $1; tot += $1; if ($1 >= 1000) k1++; if ($1 >= 5000) k5++; if ($1 >= 10000) k10++ }
           END {
             half = tot / 2; run = 0
             for (i = 1; i <= NR; i++) { run += l[i]; if (run >= half) { n50 = l[i]; break } }
             printf "    Contigs:            %d\n",       NR
             printf "    Total length:       %.2f Mbp\n", tot / 1e6
             printf "    Largest contig:     %d bp\n",    l[1]
             printf "    N50:                %d bp\n",    n50
             printf "    Contigs >= 1 kb:    %d\n",       k1
             printf "    Contigs >= 5 kb:    %d\n",       k5
             printf "    Contigs >= 10 kb:   %d\n",       k10
           }'
  echo
}

# ---------------------------- Settings summary ------------------------------
log "========== Soil metagenome assembly =========="
log "R1:          $R1"
log "R2:          $R2"
log "Long reads:  ${LONG:-(none)}${LONG:+ ($LONG_TYPE)}"
log "metaMDBG:    $( (( RUN_METAMDBG )) && echo yes || echo no )"
log "Output:      $OUTDIR"
log "Tmp:         $TMP_DIR"
log "Threads:     $THREADS"
log "RAM:         ${MEM_GB} GB (metaSPAdes), ${BB_MEM} GB (BBTools)"
log "k-mers:      $KMERS"
log "bbnorm:      target=$NORM_TARGET min=$NORM_MIN"
log "fastp:       $( (( SKIP_FASTP )) && echo 'no (-S, reads already trimmed)' || echo yes )"
log "Free disk space:"
df -h "$OUTDIR" "$TMP_DIR" | sed 's/^/    /'

# ---------------------------- 1. fastp --------------------------------------
if (( SKIP_FASTP )); then
  # Reads are already trimmed: just link to the input files, no copying
  ext="fq"; [[ "$R1" == *.gz ]] && ext="fq.gz"
  QC_R1="$QC/R1.input.$ext"
  QC_R2="$QC/R2.input.$ext"
  ln -sf "$R1" "$QC_R1"
  ln -sf "$R2" "$QC_R2"
  log "Skipping 'fastp' (-S): using the reads as they are"
else
QC_R1="$QC/R1.trim.fq.gz"
QC_R2="$QC/R2.trim.fq.gz"
run_step fastp fastp \
  -i "$R1" -I "$R2" \
  -o "$QC_R1" -O "$QC_R2" \
  --detect_adapter_for_pe \
  --cut_right --cut_right_mean_quality 20 \
  --length_required 50 \
  --thread "$FASTP_THREADS" \
  --html "$QC/fastp.html" --json "$QC/fastp.json"
fi

# ---------------------------- 2. tadpole (error correction) -----------------
# prefilter=2: keeps k-mers with depth <= 2 out of the main table, greatly reducing memory
run_step tadpole tadpole.sh "-Xmx${BB_MEM}g" \
  in="$QC_R1" in2="$QC_R2" \
  out="$EC/R1.ec.fq.gz" out2="$EC/R2.ec.fq.gz" \
  mode=correct prefilter=2 \
  threads="$THREADS" overwrite=t

# ---------------------------- 3. bbnorm (normalization) ---------------------
run_step bbnorm bbnorm.sh "-Xmx${BB_MEM}g" \
  in="$EC/R1.ec.fq.gz" in2="$EC/R2.ec.fq.gz" \
  out="$NORM/R1.norm.fq.gz" out2="$NORM/R2.norm.fq.gz" \
  target="$NORM_TARGET" min="$NORM_MIN" \
  prefilter=t \
  hist="$NORM/khist_before.txt" histout="$NORM/khist_after.txt" \
  threads="$THREADS" overwrite=t

# ---------------------------- 4. metaMDBG (HiFi only) -----------------------
# Runs before metaSPAdes: it is fast and needs little memory, so you get a
# HiFi assembly even if metaSPAdes fails.
if (( RUN_METAMDBG )); then
  mkdir -p "$MDBG"
  run_step metamdbg metaMDBG asm \
    --out-dir "$MDBG" \
    --in-hifi "$LONG" \
    --threads "$THREADS"
fi

# ---------------------------- 5. metaSPAdes ---------------------------------
if [[ ! -f "$CKDIR/metaspades.done" && -f "$ASM/params.txt" ]]; then
  # An unfinished run exists: resume from the last checkpoint.
  # --restart-from allows changing -t and -m (e.g. after running out of memory).
  log "Found an unfinished metaSPAdes run, resuming from the last checkpoint"
  run_step metaspades metaspades.py \
    --restart-from last -t "$THREADS" -m "$MEM_GB" -o "$ASM"
else
  SPADES_ARGS=(
    -1 "$NORM/R1.norm.fq.gz" -2 "$NORM/R2.norm.fq.gz"
    --only-assembler
    -k "$KMERS"
    -t "$THREADS" -m "$MEM_GB"
    --tmp-dir "$TMP_DIR"
    -o "$ASM"
  )
  if [[ -n "$LONG" ]]; then
    SPADES_ARGS+=( "$SPADES_LONG_OPT" "$LONG" )
  fi
  run_step metaspades metaspades.py "${SPADES_ARGS[@]}"
fi

# ---------------------------- 6. Statistics ---------------------------------
log "Assembly statistics:"
{
  assembly_stats "metaSPAdes (hybrid)" "$ASM/contigs.fasta"
  (( RUN_METAMDBG )) && assembly_stats "metaMDBG (HiFi-only)" "$MDBG/contigs.fasta.gz"
  true
} | tee "$OUTDIR/assembly_stats.txt"

log "Done!"
log "metaSPAdes contigs:   $ASM/contigs.fasta"
log "metaSPAdes scaffolds: $ASM/scaffolds.fasta"
(( RUN_METAMDBG )) && log "metaMDBG contigs:     $MDBG/contigs.fasta.gz"
log "Statistics:           $OUTDIR/assembly_stats.txt"
