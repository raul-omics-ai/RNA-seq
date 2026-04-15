# 🧬 RNA-seq Analysis

This repository contains a **complete pipeline for RNA-se data analysis**, from raw reads preparation to the biological interpretation of results and many cool plots.

The goal is to provide a **reproducible, automated and modular** workflow, organized into three main blocks:

<pre> ``` 
  RNA-seq-analysis/ 
  │── Upstream Analysis/ # Scripts en bash para preprocesamiento y mapeo 
  │ ├── 00_downloead_fastq_files_sra_v5.sh
  │ ├── 01_FASTQC_MULTIQUC_Quality_Control_Parallel_v3.sh
  │ ├── 02_rnaseq_trim_filtering_v4.sh
  │ ├── 03_Alig_Quant_RSEM_STAR_v4.sh
  │ ├── rename_fastq.sh
  │ └── cutadapt_batch_pe.sh
  |
  │── Differential Expression Analysis/ # Script en R para análisis de expresión diferencial 
  │ ├── Differential_Expression_Analysis_Limma_voom.R
  | ├── Differencial_Expression_Analysis_EdgeR_QL.R
  | ├── QC_Count_Matrix_RNAseq.R
  │ └── full_rnaseq_automate_analysis.R
  |
  │── Downstream Analysis/ # Scripts en R para análisis biológico 
  | ├── CIBERSORT.R
  | ├── cibersort_automatic_function.R
  │ ├── crear_kegg_database_V1.R
  │ ├── enrichment_analysis.R 
  │ ├── volcano_gsea.R
  │ └── spia_topology_based_enrichment_analysis.R
  |
  └── README.md 
  ``` </pre>

---
## 🔬 Upstream Analysis
Includes steps for
- Quality control (FastQC, MultiQC)  
- Adapter trimming and filtering (fastp, cutadapt)  
- Mapping to a reference genome (RSEM + STAR)  
- BONUS: Automatic download of public data from SRA (fasterq-dump)  

---

## 📊 Differential Expression Analysis
R Script (`full_rnaseq_automate_analysis.R`) that:
- Evaluates count matrix quality  
- Detects outliers and problematic samples  
- Runs differential expression analysis (edgeR, limma-voom)  
- Generates all diagnostic plots and tabulated results  
- Integrates enrichment analysis  

---

## 🧩 Downstream Analysis
Additional functions for:
- **Over-Representation Analysis (ORA)**  
- **Gene Set Enrichment Analysis (GSEA)**  
- **Signaling Pathway Impact Analysis (SPIA)**
- **Immune profile for bulk RNA-seq with CIBERSORT**
---

## 🚀 How to use this repository

### Requirements
- **Bash** (v5+)
- **R** (4.1)

### Environments
- rnaseq_preprocessing.yml
- r_env_4.1.0.yml
- rsem_1.3.0.yml

### Installation with Conda
``` bash 
conda env create -f rnaseq_preprocessing.yml
conda activate rnaseq_preprocessing
