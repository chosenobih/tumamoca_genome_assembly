#!/usr/bin/env bash
set -euo pipefail

########## USER SETTINGS ##########
# Raw input paths
HIFI_READS="/data2/service_center_analyses/tumamoca_genome/raw_data/hifi/tumamoca_hifi.fastq"
HIC_R1="/data2/service_center_analyses/tumamoca_genome/raw_data/hic/tumamoca_1.fastq.gz"
HIC_R2="/data2/service_center_analyses/tumamoca_genome/raw_data/hic/tumamoca_2.fastq.gz"

# Run name & lineage for BUSCO/Merqury helper scripts
RUN_TAG="subset_run_results_11_03_25"
BUSCO_LINEAGE="embryophyta_odb10"

# Threads / memory knobs
CPU_BUILD=150          # hifiasm & heavy steps
CPU_MAP_SHORT=100       # minimap2/samtools for smaller set
CPU_MAP_LONG=150       # minimap2/samtools for larger set
CPU_KMC=100
CPU_KMC_TOOLS=50
CPU_RACON=150
MM2_INDEX_SPLIT="150G" # minimap2 -I
KMC_K=21
KMC_MEM=100             # GB for KMC (-m)
USE_SRUN=0             # 1 to prefix heavy jobs with srun -c <CPUS>, 0 to run directly
###################################

# Helper for optional srun
run() {
  local cpus="$1"; shift
  if [[ "$USE_SRUN" -eq 1 ]]; then
    srun -c "${cpus}" "$@"
  else
    "$@"
  fi
}

# Create and enter working directory
mkdir -p "${RUN_TAG}"
cd "${RUN_TAG}"

echo "=== Sanity checks ==="
for p in "$HIFI_READS" "$HIC_R1" "$HIC_R2"; do
  [[ -s "$p" ]] || { echo "Missing input: $p" >&2; exit 1; }
done

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
  run "$CPU_BUILD" hifiasm -o asm_random/asm_random -t "$CPU_BUILD" --dual-scaf --telo-m TTTAGGG -l3 \
    --h1 "$HIC_R1" --h2 "$HIC_R2" hifi_random_1_3.fastq > asm_random/asm_random.log
  awk '/^S/{print ">"$2;print $3}' asm_random/asm_random.hic.hap1.p_ctg.gfa > asm_random/asm_random.hic.hap1.p_ctg.fa
  awk '/^S/{print ">"$2;print $3}' asm_random/asm_random.hic.hap2.p_ctg.gfa > asm_random/asm_random.hic.hap2.p_ctg.fa
  cat asm_random/asm_random.hic.hap1.p_ctg.fa asm_random/asm_random.hic.hap2.p_ctg.fa > asm_random/asm_random.combined.p_ctg.fa
fi

echo "=== Step 4: Assemble >11kb subset with Hi-C and -l3 ==="
mkdir -p asm_long
if [[ ! -s asm_long/asm_long.hic.hap1.p_ctg.fa ]]; then
  run "$CPU_BUILD" hifiasm -o asm_long/asm_long -t "$CPU_BUILD" --dual-scaf --telo-m TTTAGGG -l3 \
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

echo "=== Step 6: Coverage plot prep for random subset ==="
run "$CPU_MAP_SHORT" bash -c \
'minimap2 -ax map-hifi -t '"$CPU_MAP_SHORT"' asm_random/asm_random.hic.hap1.p_ctg.fa hifi_random_1_3.fastq | samtools sort -@'"$CPU_MAP_SHORT"' -o asm_random/aln_random.hic.hap1.bam'
samtools depth -a asm_random/aln_random.hic.hap1.bam > asm_random/coverage_hic.hap1.txt

run "$CPU_MAP_SHORT" bash -c \
'minimap2 -ax map-hifi -t '"$CPU_MAP_SHORT"' asm_random/asm_random.hic.hap2.p_ctg.fa hifi_random_1_3.fastq | samtools sort -@'"$CPU_MAP_SHORT"' -o asm_random/aln_random.hic.hap2.bam'
samtools depth -a asm_random/aln_random.hic.hap2.bam > asm_random/coverage_hic.hap2.txt

run "$CPU_MAP_SHORT" bash -c \
'minimap2 -ax map-hifi -t '"$CPU_MAP_SHORT"' asm_random/asm_random.combined.p_ctg.fa hifi_random_1_3.fastq | samtools sort -@'"$CPU_MAP_SHORT"' -o asm_random/aln_random.combined.bam'
samtools depth -a asm_random/aln_random.combined.bam > asm_random/coverage_combined.txt

echo "=== Step 7: Coverage plot prep for >11kb subset ==="
run "$CPU_MAP_LONG" bash -c \
'minimap2 -ax map-hifi -t '"$CPU_MAP_LONG"' asm_long/asm_long.hic.hap1.p_ctg.fa hifi_long_11kb_2_3.fastq | samtools sort -@'"$CPU_MAP_LONG"' -o asm_long/aln_long.hic.hap1.bam'
samtools depth -a asm_long/aln_long.hic.hap1.bam > asm_long/coverage_hic.hap1.txt

run "$CPU_MAP_LONG" bash -c \
'minimap2 -ax map-hifi -t '"$CPU_MAP_LONG"' asm_long/asm_long.hic.hap2.p_ctg.fa hifi_long_11kb_2_3.fastq | samtools sort -@'"$CPU_MAP_LONG"' -o asm_long/aln_long.hic.hap2.bam'
samtools depth -a asm_long/aln_long.hic.hap2.bam > asm_long/coverage_hic.hap2.txt

run "$CPU_MAP_LONG" bash -c \
'minimap2 -ax map-hifi -t '"$CPU_MAP_LONG"' asm_long/asm_long.combined.p_ctg.fa hifi_long_11kb_2_3.fastq | samtools sort -@'"$CPU_MAP_LONG"' -o asm_long/aln_long.combined.bam'
samtools depth -a asm_long/aln_long.combined.bam > asm_long/coverage_combined.txt

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
################################################################################

polish_one() {
  local label="$1"              # random | long
  local reads="$2"              # hifi_random_1_3.fastq | hifi_long_11kb_2_3.fastq
  local asm_dir="$3"            # asm_random | asm_long
  local combined="${asm_dir}/${asm_dir}.combined.p_ctg.fa"
  local outdir="polish_${label}"
  local name="tumamoca_${label}"  # NAME for KMC files, histogram labels, etc.

  echo "=== [${label}] KMC histogram on read subset (${reads}) ==="
  mkdir -p "${outdir}" "${outdir}/${name}_kmc"
  run "$CPU_KMC" kmc -k"${KMC_K}" -m"${KMC_MEM}" -t"${CPU_KMC}" "${reads}" "${outdir}/${name}" "${outdir}/${name}_kmc"
  run "$CPU_KMC_TOOLS" kmc_tools -t"${CPU_KMC_TOOLS}" transform "${outdir}/${name}" histogram "${outdir}/${name}_hist.txt"

  echo "=== [${label}] Map reads to combined assembly for polishing ==="
  run "$CPU_RACON" bash -c \
  'minimap2 -x map-pb -I '"$MM2_INDEX_SPLIT"' -a -t '"$CPU_RACON"' -c --secondary=no '"$combined"' '"$reads"' > '"${outdir}/${name}_HiFi_combined.sam"''

  echo "=== [${label}] Run racon polishing (1 round) ==="
  mkdir -p "${outdir}/racon"
  run "$CPU_RACON" racon -t "$CPU_RACON" "$reads" "${outdir}/${name}_HiFi_combined.sam" "$combined" > "${outdir}/racon/${name}.combined_racon.fa"

  echo "=== [${label}] Assembly assessment wrappers (if present) ==="
  # Expected helper scripts path: ./scripts/*.sh (relative to this run dir)
  local ASSEMBLY="${name}.combined_racon"
  local ASSEMBLY_PATH="${outdir}/racon"

  if [[ -x ../scripts/assembly_assessment_1.sh ]]; then
    sh ../scripts/assembly_assessment_1.sh "${ASSEMBLY}" "${BUSCO_LINEAGE}" "${ASSEMBLY_PATH}"
  else
    echo "  (skip) scripts/assembly_assessment_1.sh not found or not executable"
  fi

  if [[ -x ../scripts/assembly_assessment_2.sh ]]; then
    sh ../scripts/assembly_assessment_2.sh "${ASSEMBLY}" "${reads}" "${ASSEMBLY_PATH}"
  else
    echo "  (skip) scripts/assembly_assessment_2.sh not found or not executable"
  fi

  if [[ -x ../scripts/assembly_assessment_3.sh ]]; then
    sh ../scripts/assembly_assessment_3.sh "${ASSEMBLY}" "${reads}" "${ASSEMBLY_PATH}"
  else
    echo "  (skip) scripts/assembly_assessment_3.sh not found or not executable"
  fi

  echo "=== [${label}] Polishing & assessments complete ==="
}

polish_one "random" "hifi_random_1_3.fastq" "asm_random"
polish_one "long"   "hifi_long_11kb_2_3.fastq" "asm_long"

echo "=== ALL DONE: Coverage, k-mers (KAT & KMC hist), assembly stats, polishing, and assessments complete. ==="

