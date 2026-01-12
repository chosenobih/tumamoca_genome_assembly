# Tumamoca (GlobeBerry) Genome Assembly

## Overview
This repository contains scripts, configuration files, and summary results used to assemble and evaluate the **Tumamoca (GlobeBerry)** genome using **PacBio HiFi** and **Hi-C** data. The goal is to produce a high-quality, reproducible reference assembly suitable for downstream comparative genomics, conservation genetics, and functional studies.

## What’s in this repo
- Reproducible command logs and scripts for each assembly stage
- Assembly evaluation outputs (e.g., BUSCO/QUAST/Merqury summaries)
- Documentation describing design choices and known issues

> Note: Large raw data and bulky outputs (FASTQ/BAM/FASTA intermediates) are intentionally **not** tracked in git. This repo focuses on **code + configs + lightweight summaries**.

## Data types
- PacBio HiFi reads (primary assembly input)
- Hi-C read pairs (scaffolding / validation)
- Ultima reads (validation)

## High-level workflow
1. **QC & read stats**
2. **Primary HiFi assembly** (hifiasm)
3. **Organelle/contaminant filtering** (using oatk)
4. **Hi-C scaffolding** (HapHiC)
5. **Polishing / curation** (racon)
6. **Assembly evaluation**
   - BUSCO completeness
   - QUAST assembly metrics
   - Merqury k-mer completeness
