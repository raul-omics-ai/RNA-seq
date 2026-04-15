#!/bin/bash

set -euo pipefail

#############################################
# Script: rename_fastq.sh
# Description:
#   Renames FASTQ files of the form:
#   sample_1.fastq.gz / sample_2.fastq.gz
#   to:
#   sample_R1_001.fastq.gz / sample_R2_001.fastq.gz
#############################################

# Default values
INPUT_DIR=""
DRY_RUN=false

usage() {
    cat <<EOF
Usage: $(basename "$0") -i <directory> [options]

Options:
  -i DIR    Directory containing FASTQ files (required)
  -n        Dry-run: show changes without renaming files
  -h        Show this help message and exit

Examples:
  $(basename "$0") -i ./fastq_files
  $(basename "$0") -i /data/reads -n
EOF
}

# -----------------------
# Argument parsing
# -----------------------
while getopts ":i:nh" opt; do
    case $opt in
        i)
            INPUT_DIR="$OPTARG"
            ;;
        n)
            DRY_RUN=true
            ;;
        h)
            usage
            exit 0
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            exit 1
            ;;
    esac
done

# -----------------------
# Checks
# -----------------------
if [[ -z "$INPUT_DIR" ]]; then
    echo "❌ Error: You must specify a directory with -i"
    usage
    exit 1
fi

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "❌ Error: Directory does not exist: $INPUT_DIR"
    exit 1
fi

# -----------------------
# Renaming
# -----------------------
shopt -s nullglob

for f in "$INPUT_DIR"/*.fastq.gz; do
    base=$(basename "$f")

    if [[ "$base" =~ _1\.fastq\.gz$ ]]; then
        sample="${base/_1.fastq.gz/}"
        new_name="${sample}_R1_001.fastq.gz"
    elif [[ "$base" =~ _2\.fastq\.gz$ ]]; then
        sample="${base/_2.fastq.gz/}"
        new_name="${sample}_R2_001.fastq.gz"
    else
        echo "⚠️  Skipping unrecognized file: $base"
        continue
    fi

    if [[ "$DRY_RUN" == true ]]; then
        echo "[DRY-RUN] mv \"$f\" \"$INPUT_DIR/$new_name\""
    else
        mv "$f" "$INPUT_DIR/$new_name"
        echo "$base → $new_name"
    fi
done

echo " Process completed."
