#!/usr/bin/env bash
set -euo pipefail

# run_assessments.sh
# Minimal wrapper to run only:
#   scripts/assembly_assessment_1.sh  <ASSEMBLY> <BUSCO_LINEAGE> <ASSEMBLY_PATH>
#   scripts/assembly_assessment_2.sh  <ASSEMBLY> <HiFi_reads>    <ASSEMBLY_PATH>
#   scripts/assembly_assessment_3.sh  <ASSEMBLY> <HiFi_reads>    <ASSEMBLY_PATH>
#
# Usage example:
#   ./run_assessments.sh \
#     --scripts-dir /path/to/scripts \
#     --assembly-path /path/to/polish_random/racon \
#     --assembly tumamoca_random.combined_racon \
#     --hifi /data/hifi_random_1_3.fastq \
#     --busco embryophyta_odb10 \
#     --conda-env myenv            # optional

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

SCRIPTS_DIR="${SCRIPT_DIR}/scripts"
ASSEMBLY_PATH=""
ASSEMBLY=""
HIFI=""
BUSCO_LINEAGE=""
CONDA_ENV=""

# --- parse args ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --scripts-dir)  SCRIPTS_DIR="$2"; shift 2;;
    --assembly-path) ASSEMBLY_PATH="$2"; shift 2;;
    --assembly)     ASSEMBLY="$2"; shift 2;;
    --hifi)         HIFI="$2"; shift 2;;
    --busco)        BUSCO_LINEAGE="$2"; shift 2;;
    --conda-env)    CONDA_ENV="$2"; shift 2;;
    -h|--help)
      cat <<EOF
Usage: $0 --scripts-dir <dir> --assembly-path <dir> --assembly <name> --hifi <reads.fastq> --busco <lineage> [--conda-env <env>]
Runs only the assessment scripts with robust path handling.
EOF
      exit 0;;
    *) echo "Unknown arg: $1" >&2; exit 1;;
  esac
done

# --- required checks ---
[[ -n "$ASSEMBLY_PATH" ]] || { echo "ERROR: --assembly-path is required"; exit 1; }
[[ -n "$ASSEMBLY"      ]] || { echo "ERROR: --assembly is required"; exit 1; }
[[ -n "$HIFI"          ]] || { echo "ERROR: --hifi is required"; exit 1; }
[[ -n "$BUSCO_LINEAGE" ]] || { echo "ERROR: --busco is required"; exit 1; }
[[ -d "$SCRIPTS_DIR"   ]] || { echo "ERROR: scripts dir not found: $SCRIPTS_DIR"; exit 1; }

# Optional: ensure conda env is active even in non-interactive shells
if [[ -n "$CONDA_ENV" ]]; then
  if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
    conda activate "$CONDA_ENV"
  else
    echo "WARNING: conda not found; skipping activation of $CONDA_ENV"
  fi
fi

A1="${SCRIPTS_DIR}/assembly_assessment_1.sh"
A2="${SCRIPTS_DIR}/assembly_assessment_2.sh"
A3="${SCRIPTS_DIR}/assembly_assessment_3.sh"

echo "ASSEMBLY_PATH : $ASSEMBLY_PATH"
echo "ASSEMBLY      : $ASSEMBLY"
echo "HiFi reads    : $HIFI"
echo "BUSCO lineage : $BUSCO_LINEAGE"
echo "Scripts dir   : $SCRIPTS_DIR"
echo

# Helper to run or emit a clear skip message
run_or_skip () {
  local script="$1"; shift
  if [[ -x "$script" ]]; then
    echo ">>> Running: $script $*"
    sh "$script" "$@"
  else
    echo "SKIP: $script not found or not executable (chmod +x ?)"
  fi
  echo
}

run_or_skip "$A1" "$ASSEMBLY" "$BUSCO_LINEAGE" "$ASSEMBLY_PATH"
run_or_skip "$A2" "$ASSEMBLY" "$HIFI"          "$ASSEMBLY_PATH"
run_or_skip "$A3" "$ASSEMBLY" "$HIFI"          "$ASSEMBLY_PATH"

echo "All assessment steps attempted."

