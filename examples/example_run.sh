#!/usr/bin/env bash

# Advanced example run script for the VCF Trio Analyzer.
#
# Features:
#   - Argument parsing (VCF path, output dir, databases dir)
#   - Checks for required tools: python, awk, bedtools
#   - Checks for required BED files in the databases directory
#   - Creates the output directory if needed
#
# Usage:
#   ./example_run.sh [-v VCF_PATH] [-o OUT_DIR] [-d DATABASES_DIR]
#
# Defaults:
#   VCF_PATH      = ../DellyVariation.vcf
#   OUT_DIR       = ../out
#   DATABASES_DIR = ../databases

set -euo pipefail

# ----------- defaults ----------- #
DEFAULT_VCF="../DellyVariation.vcf"
DEFAULT_OUT="../out"
DEFAULT_DB_DIR="../databases"

VCF_PATH="${DEFAULT_VCF}"
OUT_DIR="${DEFAULT_OUT}"
DB_DIR="${DEFAULT_DB_DIR}"

# Required BED files (filenames only, expected inside DB_DIR)
REQUIRED_BEDS=(
  "hg38_genes.bed"
  "hg38_exons.bed"
  "ClinGen_haploinsufficiency_gene_GRCh38.bed"
  "ClinGen_triplosensitivity_gene_GRCh38.bed"
  "ClinGen recurrent CNV-hg38.bed"
  "clinvar_SV_clean.bed"
)

# ----------- helper functions ----------- #

usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Options:
  -v  PATH   Path to multi-sample VCF file (default: ${DEFAULT_VCF})
  -o  PATH   Output directory (default: ${DEFAULT_OUT})
  -d  PATH   Databases directory containing BED files (default: ${DEFAULT_DB_DIR})
  -h        Show this help and exit

Example:
  $(basename "$0") -v ../DellyVariation.vcf -o ../out -d ../databases
EOF
}

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

check_command() {
  local cmd="$1"
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    fail "Required command '${cmd}' not found in PATH."
  fi
}

# ----------- parse arguments ----------- #

while getopts ":v:o:d:h" opt; do
  case "${opt}" in
    v) VCF_PATH="${OPTARG}" ;;
    o) OUT_DIR="${OPTARG}" ;;
    d) DB_DIR="${OPTARG}" ;;
    h) usage; exit 0 ;;
    \?) fail "Invalid option: -${OPTARG} (use -h for help)" ;;
    :) fail "Option -${OPTARG} requires an argument." ;;
  esac
done

# ----------- sanity checks ----------- #

echo ">>> Checking required commands..."
check_command python
check_command awk
check_command bedtools
echo "OK: python, awk, bedtools found."
echo

echo ">>> Checking VCF file..."
if [[ ! -f "${VCF_PATH}" ]]; then
  fail "VCF file not found at: ${VCF_PATH}"
fi
echo "OK: VCF file exists: ${VCF_PATH}"
echo

echo ">>> Checking databases directory and BED files..."
if [[ ! -d "${DB_DIR}" ]]; then
  fail "Databases directory not found: ${DB_DIR}"
fi

missing_any=false
for bed in "${REQUIRED_BEDS[@]}"; do
  path="${DB_DIR}/${bed}"
  if [[ ! -f "${path}" ]]; then
    echo "MISSING: ${path}"
    missing_any=true
  else
    echo "FOUND  : ${path}"
  fi
done

if [[ "${missing_any}" == true ]]; then
  fail "One or more required BED files are missing. See README and databases/README.md."
fi
echo "OK: All required BED files found."
echo

echo ">>> Preparing output directory..."
mkdir -p "${OUT_DIR}"
echo "Output directory: ${OUT_DIR}"
echo

# ----------- run the analyzer ----------- #

echo ">>> Running VCF Trio Analyzer"
echo "VCF        : ${VCF_PATH}"
echo "OUT        : ${OUT_DIR}"
echo "DATABASES  : ${DB_DIR}"
echo

# Optionally, you can export DB_DIR so vcf_analyzer.py can use it,
# if you later modify the script to read from this env variable.
export VCF_TRIO_DB_DIR="${DB_DIR}"

python ../src/vcf_analyzer.py \
  --vcf "${VCF_PATH}" \
  --out "${OUT_DIR}"

echo
echo ">>> Done!"
echo "Results written to: ${OUT_DIR}"