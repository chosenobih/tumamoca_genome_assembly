#!/usr/bin/env bash
set -euo pipefail

########## USER SETTINGS ##########
# Raw input paths
HIFI_READS="/data2/service_center_analyses/tumamoca_genome/raw_data/hifi/tumamoca_hifi.fastq"
HIC_R1="/data2/service_center_analyses/tumamoca_genome/raw_data/hic/tumamoca_1.fastq.gz"
HIC_R2="/data2/service_center_analyses/tumamoca_genome/raw_data/hic/tumamoca_2.fastq.gz"

# Run tag (output folder)
RUN_TAG="subset_run_results_11_06_25"

# BUSCO lineage for assessments
BUSCO_LINEAGE="embryophyta_odb10"

# If you want the script to activate a conda env itself, set this (else leave empty)
CONDA_ENV=""   # e.g., "genome-env" or leave "" to skip

# Threads / memory knobs
CPU_BUILD=150          # hifiasm & heavy steps
CPU_MAP_SHORT=100       # minimap2/samtools for smaller set
CPU_MAP_LONG=150       # minimap2/samtools for larger set
CPU_KMC=100
CPU_KMC_TOOLS=50
CPU_RACON=150
MM2_INDEX_SPLIT="128G" # minimap2 -I (index split size)
KMC_K=21
KMC_MEM=150             # GB for KMC (-m)

# Slurm toggle: 1 uses srun; 0 runs directly
USE_SRUN=0
###################################

# Resolve script location and helper scripts dir (works even after cd)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS_DIR="${SCRIPT_DIR}/scripts"   # change if your helpers live elsewhere

# Optional: ensure conda env is active even in non-interactive shells
if [[ -n "$CONDA_ENV" ]]; then
  if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
    conda activate "$CONDA_ENV"
  else
    echo "WARNING: conda not found; skipping activation of $CONDA_ENV"
  fi
fi

# srun wrapper
run() {
  local cpus="$1"; shift
  if [[ "$USE_SRUN" -eq 1 ]]; then
    srun -c "${cpus}" "$@"
  else
    "$@"
  fi
}

echo "=== Sanity checks ==="
for p in "$HIFI_READS" "$HIC_R1" "$HIC_R2"; do
  [[ -s "$p" ]] || { echo "Missing input: $p" >&2; exit 1; }
done

# Create/enter working dir
mkdir -p "${RUN_TAG}"
cd "${RUN_TAG}"

echo "=== Step 1: Randomly sample one-third of HiFi reads ==="
if [[ ! -s hifi_random_1_3.fastq ]]; then
  seqtk sample -s100 "$HIFI_READS" 0.33 > hifi_random_1_3.fastq
fi

echo "=== Step 2: Filter HiFi reads >11kb and sample two-thirds ==="
if [[ ! -s hifi_over11kb.fastq ]]; then
  seqkit seq -m 11000 "$HIFI_READS" > hifi_over11kb.fastq
fi
if [[ ! -s hifi_long_11kb_2_3.fastq ]]; then
  seqtk sample -s100 hifi_over11kb.fastq 0.66 > hifi_long_11kb_2_3.fastq
fi

echo "=== Step 3: Assemble random subset with Hi-C and -l3 ==="
mkdir -p asm_random
if [[ ! -s asm_random/asm_random.hic.hap1.p_ctg.fa ]]; then
  run "$CPU_BUILD" hifiasm -o asm_random/asm_random -t "$CPU_BUILD" -l3 \
    --h1 "$HIC_R1" --h2 "$HIC_R2" hifi_random_1_3.fastq > asm_random/asm_random.log
  awk '/^S/{print ">"$2;print $3}' asm_random/asm_random.hic.hap1.p_ctg.gfa > asm_random/asm_random.hic.hap1.p_ctg.fa
  awk '/^S/{print ">"$2;print $3}' asm_random/asm_random.hic.hap2.p_ctg.gfa > asm_random/asm_random.hic.hap2.p_ctg.fa
  cat asm_random/asm_random.hic.hap1.p_ctg.fa asm_random/asm_random.hic.hap2.p_ctg.fa > asm_random/asm_random.combined.p_ctg.fa
fi

echo "=== Step 4: Assemble >11kb subset with Hi-C and -l3 ==="
mkdir -p asm_long
if [[ ! -s asm_long/asm_long.hic.hap1.p_ctg.fa ]]; then
  run "$CPU_BUILD" hifiasm -o asm_long/asm_long -t "$CPU_BUILD" -l3 \
    --h1 "$HIC_R1" --h2 "$HIC_R2" hifi_long_11kb_2_3.fastq > asm_long/asm_long.log
  awk '/^S/{print ">"$2;print $3}' asm_long/asm_long.hic.hap1.p_ctg.gfa > asm_long/asm_long.hic.hap1.p_ctg.fa
  awk '/^S/{print ">"$2;print $3}' asm_long/asm_long.hic.hap2.p_ctg.gfa > asm_long/asm_long.hic.hap2.p_ctg.fa
  cat asm_long/asm_long.hic.hap1.p_ctg.fa asm_long/asm_long.hic.hap2.p_ctg.fa > asm_long/asm_long.combined.p_ctg.fa
fi

echo "=== Step 5: Assembly stats (N50, size, etc.) ==="
assembly-stats asm_random/asm_random.hic.hap1.p_ctg.fa > asm_random/contig_stats.hic.hap1.txt
assembly-stats asm_random/asm_random.hic.hap2.p_ctg.fa > asm_random/contig_stats.hic.hap2.txt
assembly-stats asm_random/asm_random.combined.p_ctg.fa > asm_random/contig_stats.combined.txt
assembly-stats asm_long/asm_long.hic.hap1.p_ctg.fa   > asm_long/contig_stats.hic.hap1.txt
assembly-stats asm_long/asm_long.hic.hap2.p_ctg.fa   > asm_long/contig_stats.hic.hap2.txt
assembly-stats asm_long/asm_long.combined.p_ctg.fa   > asm_long/contig_stats.combined.txt

###############################################################################
# Step 6–7: COVERAGE (hap1, hap2, combined) for random and >11kb assemblies
###############################################################################
mm2_depth() {
  local threads="$1"
  local fasta="$2"
  local reads="$3"
  local bam_out="$4"
  local depth_out="$5"

  run "$threads" bash -c \
    'minimap2 -ax map-hifi -t '"$threads"' '"$fasta"' '"$reads"' \
    | samtools sort -@'"$threads"' -o '"$bam_out"''
  samtools index -@ "$threads" "$bam_out"
  samtools depth -a "$bam_out" > "$depth_out"
}

echo "=== Step 6: Coverage for RANDOM assembly (hap1, hap2, combined) ==="
mm2_depth "$CPU_MAP_SHORT" "asm_random/asm_random.hic.hap1.p_ctg.fa" hifi_random_1_3.fastq \
          "asm_random/aln_random.hic.hap1.bam" "asm_random/coverage_hic.hap1.txt"
mm2_depth "$CPU_MAP_SHORT" "asm_random/asm_random.hic.hap2.p_ctg.fa" hifi_random_1_3.fastq \
          "asm_random/aln_random.hic.hap2.bam" "asm_random/coverage_hic.hap2.txt"
mm2_depth "$CPU_MAP_SHORT" "asm_random/asm_random.combined.p_ctg.fa" hifi_random_1_3.fastq \
          "asm_random/aln_random.combined.bam" "asm_random/coverage_combined.txt"

echo "=== Step 7: Coverage for >11kb assembly (hap1, hap2, combined) ==="
mm2_depth "$CPU_MAP_LONG" "asm_long/asm_long.hic.hap1.p_ctg.fa" hifi_long_11kb_2_3.fastq \
          "asm_long/aln_long.hic.hap1.bam" "asm_long/coverage_hic.hap1.txt"
mm2_depth "$CPU_MAP_LONG" "asm_long/asm_long.hic.hap2.p_ctg.fa" hifi_long_11kb_2_3.fastq \
          "asm_long/aln_long.hic.hap2.bam" "asm_long/coverage_hic.hap2.txt"
mm2_depth "$CPU_MAP_LONG" "asm_long/asm_long.combined.p_ctg.fa" hifi_long_11kb_2_3.fastq \
          "asm_long/aln_long.combined.bam" "asm_long/coverage_combined.txt"

echo "=== Step 8: Run KAT comp on random assembly ==="
run "$CPU_BUILD" kat comp -o kat_random -t "$CPU_BUILD" hifi_random_1_3.fastq asm_random/asm_random.hic.hap1.p_ctg.fa
run "$CPU_BUILD" kat comp -o kat_random -t "$CPU_BUILD" hifi_random_1_3.fastq asm_random/asm_random.hic.hap2.p_ctg.fa
run "$CPU_BUILD" kat comp -o kat_random -t "$CPU_BUILD" hifi_random_1_3.fastq asm_random/asm_random.combined.p_ctg.fa

echo "=== Step 9: Run KAT comp on >11kb assembly ==="
run "$CPU_BUILD" kat comp -o kat_long -t "$CPU_BUILD" hifi_long_11kb_2_3.fastq asm_long/asm_long.hic.hap1.p_ctg.fa
run "$CPU_BUILD" kat comp -o kat_long -t "$CPU_BUILD" hifi_long_11kb_2_3.fastq asm_long/asm_long.hic.hap2.p_ctg.fa
run "$CPU_BUILD" kat comp -o kat_long -t "$CPU_BUILD" hifi_long_11kb_2_3.fastq asm_long/asm_long.combined.p_ctg.fa

################################################################################
# NEW STEPS: KMC → minimap2 (map-pb) → racon polishing → assessment (both sets)
# (now for hap1, hap2, and combined; plus polished coverage)
################################################################################

# Map reads to FASTA using minimap2 map-pb (SAM output for racon)
mm2_map_pb_to_sam() {
  local threads="$1"
  local fasta="$2"
  local reads="$3"
  local sam_out="$4"
  run "$threads" bash -c \
    'minimap2 -x map-pb -I '"$MM2_INDEX_SPLIT"' -a -t '"$threads"' -c --secondary=no '"$fasta"' '"$reads"' > '"$sam_out"''
}

# Coverage helper using map-pb
mm2pb_depth() {
  local threads="$1"
  local fasta="$2"
  local reads="$3"
  local bam_out="$4"
  local depth_out="$5"
  run "$threads" bash -c \
    'minimap2 -x map-pb -t '"$threads"' '"$fasta"' '"$reads"' \
    | samtools sort -@'"$threads"' -o '"$bam_out"''
  samtools index -@ "$threads" "$bam_out"
  samtools depth -a "$bam_out" > "$depth_out"
}

# For each label (random/long), KMC once, then polish+assess hap1/hap2/combined
polish_and_assess_all() {
  local label="$1"              # random | long
  local reads="$2"              # hifi_random_1_3.fastq | hifi_long_11kb_2_3.fastq
  local asm_dir="$3"            # asm_random | asm_long

  # Input FASTAs
  local fasta_h1="${asm_dir}/${asm_dir}.hic.hap1.p_ctg.fa"
  local fasta_h2="${asm_dir}/${asm_dir}.hic.hap2.p_ctg.fa"
  local fasta_cmb="${asm_dir}/${asm_dir}.combined.p_ctg.fa"

  # Base outdir + KMC (run once per subset)
  local outbase="polish_${label}"
  mkdir -p "${outbase}" "${outbase}/kmc"
  echo "=== [${label}] KMC histogram on read subset (${reads}) ==="
  run "$CPU_KMC" kmc -k"${KMC_K}" -m"${KMC_MEM}" -t"${CPU_KMC}" "${reads}" "${outbase}/kmc/${label}" "${outbase}/kmc/tmp"
  run "$CPU_KMC_TOOLS" kmc_tools -t"${CPU_KMC_TOOLS}" transform "${outbase}/kmc/${label}" histogram "${outbase}/kmc/${label}_hist.txt"

  # Array of targets: name | input_fasta
  local targets=(
    "hap1|${fasta_h1}"
    "hap2|${fasta_h2}"
    "combined|${fasta_cmb}"
  )

  for spec in "${targets[@]}"; do
    IFS='|' read -r tname tfasta <<< "$spec"
    local name="tumamoca_${label}.${tname}"
    local tdir="${outbase}/${tname}"

    echo "=== [${label} :: ${tname}] Polishing workflow ==="
    mkdir -p "${tdir}"

    # 1) Map (map-pb) to SAM for racon
    mm2_map_pb_to_sam "$CPU_RACON" "$tfasta" "$reads" "${tdir}/${name}.HiFi.sam"

    # 2) Racon polishing (1 round)
    mkdir -p "${tdir}/racon"
    run "$CPU_RACON" racon -t "$CPU_RACON" "$reads" "${tdir}/${name}.HiFi.sam" "$tfasta" \
      > "${tdir}/racon/${name}.racon.fa"

    # 3) Coverage on polished FASTA (map-pb)
    mm2pb_depth "$CPU_RACON" "${tdir}/racon/${name}.racon.fa" "$reads" \
                "${tdir}/racon/${name}.racon.bam" \
                "${tdir}/racon/${name}.racon.coverage.txt"

    # 4) Assessments (if helper scripts exist)
    local ASSEMBLY="${name}.racon"
    local ASSEMBLY_PATH="${tdir}/racon"

    if [[ -x "${SCRIPTS_DIR}/assembly_assessment_1.sh" ]]; then
      sh "${SCRIPTS_DIR}/assembly_assessment_1.sh" "${ASSEMBLY}" "${BUSCO_LINEAGE}" "${ASSEMBLY_PATH}"
    else
      echo "  (skip) ${SCRIPTS_DIR}/assembly_assessment_1.sh not found or not executable"
    fi

    if [[ -x "${SCRIPTS_DIR}/assembly_assessment_2.sh" ]]; then
      sh "${SCRIPTS_DIR}/assembly_assessment_2.sh" "${ASSEMBLY}" "${reads}" "${ASSEMBLY_PATH}"
    else
      echo "  (skip) ${SCRIPTS_DIR}/assembly_assessment_2.sh not found or not executable"
    fi

    if [[ -x "${SCRIPTS_DIR}/assembly_assessment_3.sh" ]]; then
      sh "${SCRIPTS_DIR}/assembly_assessment_3.sh" "${ASSEMBLY}" "${reads}" "${ASSEMBLY_PATH}"
    else
      echo "  (skip) ${SCRIPTS_DIR}/assembly_assessment_3.sh not found or not executable"
    fi

    echo "=== [${label} :: ${tname}] Polishing, coverage, and assessments complete ==="
  done
}

polish_and_assess_all "random" "hifi_random_1_3.fastq" "asm_random"
polish_and_assess_all "long"   "hifi_long_11kb_2_3.fastq" "asm_long"

echo "=== ALL DONE: coverage (hap1/hap2/combined) + KMC + racon + polished coverage + assessments for both assemblies. ==="

