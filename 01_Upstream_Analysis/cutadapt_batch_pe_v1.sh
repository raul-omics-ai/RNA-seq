#!/usr/bin/env bash
# cutadapt_batch_pe.sh
# Batch adapter trimming for paired-end FASTQ files using cutadapt
# Requirements: cutadapt (https://cutadapt.readthedocs.io), bash >= 4.3 (for wait -n)
# Created for: Raúl

set -euo pipefail
IFS=$'\n\t'

print_help() {
  cat <<'EOF'
Usage: cutadapt_batch_pe.sh -i INPUT -o OUTDIR [options]

Required:
  -i INPUT       Input: directory containing FASTQ pairs, a glob (e.g. "data/*_R1.fastq.gz")
                 or a text file listing one sample per line. Lines may be:
                   /full/path/sample_R1.fastq.gz    (script will infer R2 by replacing _R1 with _R2 or _1 with _2)
                   /path/R1.fastq.gz /path/R2.fastq.gz   (two columns: explicit pair)
  -o OUTDIR      Output directory (will be created if missing)

Options:
  -a ADAPTER_FWD Adapter for forward reads (passed to -a)
  -A ADAPTER_REV Adapter for reverse reads (passed to -A)
  -m MINLEN      Minimum length to keep (cutadapt -m). Default: 20
  -q QUALITY     Quality cutoff for trimming (cutadapt -q). Default: 20
  -t THREADS     Threads per cutadapt job (-j). Default: 1
  -j JOBS        Number of parallel samples to run simultaneously. Default: 1
  -e EXTRA       Additional extra arguments to append to cutadapt command (quoted)
  -n             Dry-run: show planned commands but do not execute
  -x             Stop-on-error: if any sample fails, kill remaining jobs and exit
  -h             Show this help and exit

Example:
  ./cutadapt_batch_pe.sh -i ./fastq/ -o ./trimmed/ -a AGATCGGAAGAGC -A AGATCGGAAGAGC -m 25 -q 20 -t 4 -j 4

EOF
}

# Default parameters
ADAPTER_FWD=""
ADAPTER_REV=""
MINLEN=20
QUALITY=20
THREADS=1
JOBS=1
EXTRA_ARGS=""
DRYRUN=0
STOP_ON_ERROR=0

while getopts ":i:o:a:A:m:q:t:j:e:nxh" opt; do
  case $opt in
    i) INPUT="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    a) ADAPTER_FWD="$OPTARG" ;;
    A) ADAPTER_REV="$OPTARG" ;;
    m) MINLEN="$OPTARG" ;;
    q) QUALITY="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    j) JOBS="$OPTARG" ;;
    e) EXTRA_ARGS="$OPTARG" ;;
    n) DRYRUN=1 ;;
    x) STOP_ON_ERROR=1 ;;
    h) print_help; exit 0 ;;
    \?) echo "Unknown option -$OPTARG" >&2; print_help; exit 1 ;;
  esac
done

if [ -z "${INPUT:-}" ] || [ -z "${OUTDIR:-}" ]; then
  echo "ERROR: -i and -o are required." >&2
  print_help
  exit 1
fi

mkdir -p "$OUTDIR"
LOGDIR="$OUTDIR/logs"
mkdir -p "$LOGDIR"

# Build list of pairs (arrays in1_list and in2_list)
declare -a in1_list in2_list sample_names

# Helper: try to infer R2 from R1 by replacing common patterns
infer_r2() {
  local r1="$1"
  if [[ "$r1" =~ _R1 ]]; then
    echo "${r1/_R1/_R2}"
  elif [[ "$r1" =~ _1\.fastq ]]; then
    echo "${r1/_1./_2.}"
  elif [[ "$r1" =~ \.1\.(fastq|fq) ]]; then
    echo "${r1/.1./.2.}"
  else
    # fallback: try replacing _R1 with _R2 anywhere
    echo "${r1/_R1/_R2}"
  fi
}

# If INPUT is a directory
if [ -d "$INPUT" ]; then
  mapfile -t candidates < <(find "$INPUT" -maxdepth 1 -type f \( -iname "*_R1*.fastq*" -o -iname "*_1.fastq*" \) | sort)
  for r1 in "${candidates[@]}"; do
    r2=$(infer_r2 "$r1")
    if [ -f "$r2" ]; then
      in1_list+=("$r1")
      in2_list+=("$r2")
      # sample name: remove path and _R1 / _1 suffix
      base=$(basename "$r1")
      sname=$(echo "$base" | sed -E 's/(_R?1|\.1)(\.|_).*//; s/(_R?1|\.1)//')
      # fallback: remove known suffixes
      sname=${sname%%.*}
      sample_names+=("$sname")
    fi
  done
elif [[ "$INPUT" == *"*"* || "$INPUT" == *"?"* ]]; then
  # glob pattern
  eval "mapfile -t matched < <(ls -1 $INPUT 2>/dev/null || true)"
  for r1 in "${matched[@]}"; do
    r2=$(infer_r2 "$r1")
    if [ -f "$r2" ]; then
      in1_list+=("$r1")
      in2_list+=("$r2")
      base=$(basename "$r1")
      sname=$(echo "$base" | sed -E 's/(_R?1|\.1)(\.|_).*//; s/(_R?1|\.1)//')
      sname=${sname%%.*}
      sample_names+=("$sname")
    fi
  done
elif [ -f "$INPUT" ]; then
  # treat as list file: each line either: r1 [r2]
  while read -r line || [ -n "$line" ]; do
    # skip empty / comments
    [[ -z "$line" || "$line" =~ ^# ]] && continue
    # split
    read -r col1 col2 <<<"$line"
    if [ -n "$col2" ]; then
      r1="$col1"; r2="$col2"
    else
      r1="$col1"; r2=$(infer_r2 "$r1")
    fi
    if [ -f "$r1" ] && [ -f "$r2" ]; then
      in1_list+=("$r1")
      in2_list+=("$r2")
      base=$(basename "$r1")
      sname=$(echo "$base" | sed -E 's/(_R?1|\.1)(\.|_).*//; s/(_R?1|\.1)//')
      sname=${sname%%.*}
      sample_names+=("$sname")
    else
      echo "Warning: pair not found or missing: $r1 $r2" >&2
    fi
  done < "$INPUT"
else
  echo "ERROR: -i value not recognized (not a dir, glob, or file): $INPUT" >&2
  exit 1
fi

TOTAL=${#in1_list[@]}
if [ "$TOTAL" -eq 0 ]; then
  echo "No paired files found. Exiting." >&2
  exit 1
fi

echo "Found $TOTAL paired samples. Output -> $OUTDIR. Using $JOBS parallel jobs, $THREADS threads per cutadapt."

# Build and optionally run commands
cmds=()
for i in "${!in1_list[@]}"; do
  in1="${in1_list[$i]}"; in2="${in2_list[$i]}"; sname="${sample_names[$i]}"
  out1="$OUTDIR/${sname}_R1.trim.fastq.gz"
  out2="$OUTDIR/${sname}_R2.trim.fastq.gz"

  # build cutadapt command
  cmd=(cutadapt)
  if [ -n "$ADAPTER_FWD" ]; then cmd+=( -a "$ADAPTER_FWD" ); fi
  if [ -n "$ADAPTER_REV" ]; then cmd+=( -A "$ADAPTER_REV" ); fi
  cmd+=( -m "$MINLEN" -q "$QUALITY" -j "$THREADS" -o "$out1" -p "$out2" )
  # add any extra args as words
  if [ -n "$EXTRA_ARGS" ]; then
    # shell-split EXTRA_ARGS conservatively
    read -r -a extra_arr <<< "$EXTRA_ARGS"
    for x in "${extra_arr[@]}"; do cmd+=( "$x" ); done
  fi
  cmd+=( "$in1" "$in2" )

  # join into a single, properly escaped string for display/execution
  cmdstr=$(printf '%q ' "${cmd[@]}")
  cmds+=("$cmdstr")
done

if [ "$DRYRUN" -eq 1 ]; then
  echo "Dry-run mode. Commands to be executed:"
  for c in "${cmds[@]}"; do echo "$c"; done
  exit 0
fi

# Run jobs with simple parallelization and progress reporting using wait -n
completed=0
failed=0

# Trap to ensure waiting and cleanup on SIGINT
trap 'echo; echo "Interrupted by user. Waiting for running jobs to finish/cancel..."; exit 1' SIGINT SIGTERM

for idx in "${!cmds[@]}"; do
  cmd=${cmds[$idx]}
  sname=${sample_names[$idx]}
  logfile="$LOGDIR/${sname}.cutadapt.log"

  echo "\n[START] ($((idx+1))/$TOTAL) $sname"
  # run command in background, redirect stdout+stderr to logfile
  bash -lc "$cmd" >"$logfile" 2>&1 &

  # enforce max parallel jobs
  while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do
    # wait for any job to end
    if wait -n; then
      completed=$((completed+1))
      pct=$(( 100 * completed / TOTAL ))
      echo "[PROGRESS] Completed: $completed/$TOTAL ($pct%)"
    else
      failed=$((failed+1))
      completed=$((completed+1))
      pct=$(( 100 * completed / TOTAL ))
      echo "[PROGRESS] Completed: $completed/$TOTAL ($pct%) -- some job failed"
      if [ "$STOP_ON_ERROR" -eq 1 ]; then
        echo "Stopping due to error (stop-on-error enabled). Killing remaining jobs..."
        pids=$(jobs -pr)
        if [ -n "$pids" ]; then
          kill $pids 2>/dev/null || true
        fi
        wait
        echo "Exiting due to failure. See logs in $LOGDIR"
        exit 2
      fi
    fi
  done

done

# wait for remaining jobs
while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
  if wait -n; then
    completed=$((completed+1))
    pct=$(( 100 * completed / TOTAL ))
    echo "[PROGRESS] Completed: $completed/$TOTAL ($pct%)"
  else
    failed=$((failed+1))
    completed=$((completed+1))
    pct=$(( 100 * completed / TOTAL ))
    echo "[PROGRESS] Completed: $completed/$TOTAL ($pct%) -- some job failed"
    if [ "$STOP_ON_ERROR" -eq 1 ]; then
      echo "Stopping due to error (stop-on-error enabled). Killing remaining jobs..."
      pids=$(jobs -pr)
      if [ -n "$pids" ]; then
        kill $pids 2>/dev/null || true
      fi
      wait
      echo "Exiting due to failure. See logs in $LOGDIR"
      exit 2
    fi
  fi
done

# All done
echo "\nAll jobs finished. Total: $TOTAL. Succeeded: $((TOTAL - failed)). Failed: $failed."
echo "Log files: $LOGDIR"

if [ "$failed" -gt 0 ]; then
  echo "Some samples failed. Inspect the corresponding logs in $LOGDIR"
  exit 2
fi

exit 0

