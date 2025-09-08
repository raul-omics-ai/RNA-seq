# 🧬 RNA-seq Analysis

Este repositorio contiene un **pipeline completo para el análisis de datos RNA-seq**, 
desde la preparación de las lecturas crudas hasta la interpretación biológica de los resultados y muchos gráficos chulos.  

El objetivo es proporcionar un flujo **reproducible, automatizado y modular**, 
organizado en tres grandes bloques:

<pre> ``` 
  RNA-seq-analysis/ 
  │── Upstream Analysis/ # Scripts en bash para preprocesamiento y mapeo 
  │ ├── 00_downloead_fastq_files_sra_v5.sh
  │ ├── 01_FASTQC_MULTIQUC_Quality_Control_Parallel_v3.sh
  │ ├── 02_rnaseq_trim_filtering_v4.sh
  │ ├── cutadapt_batch_pe.sh
  │ └── 03_Alig_Quant_RSEM_STAR_v4.sh
  │── Differential Expression Analysis/ # Script en R para análisis de expresión diferencial 
  │ ├── Differential_Expression_Analysis_Limma_voom_v4.R
  | ├── Differencial_Expression_Analysis_EdgeR_QL_v3.R
  | ├── QC_Count_Matrix_RNAseq_v6.R
  │ └── full_rnaseq_automate_analysis_v4.R
  │── Downstream Analysis/ # Scripts en R para análisis biológico 
  │ ├── crear_kegg_database_V1.R
  │ ├── enrichment_analysis_v4.R 
  │ ├── volcano_gsea_v1.R
  │ └── spia_topology_based_enrichment_analysis_V1.R 
  └── README.md 
  ``` </pre>

---
## 🔬 Upstream Analysis
Incluye los pasos para:
- Control de calidad (FastQC, MultiQC)  
- Filtrado y recorte de adaptadores (fastp, cutadapt)  
- Mapeo a genoma de referencia (RSEM + STAR)  
- Descarga automática de datos públicos desde SRA  

---

## 📊 Differential Expression Analysis
Script en R (`full_rnaseq_automate_analysis_v4.R`) que:
- Evalúa calidad de la matriz de cuentas  
- Detecta outliers y muestras problemáticas  
- Corre análisis de expresión diferencial (edgeR, limma-voom)  
- Genera todos los gráficos de diagnóstico y resultados tabulados  
- Integra análisis de enriquecimiento  

---

## 🧩 Downstream Analysis
Funciones adicionales para:
- **Over-Representation Analysis (ORA)**  
- **Gene Set Enrichment Analysis (GSEA)**  
- **Signaling Pathway Impact Analysis (SPIA)** 
---

## 🚀 Cómo usar este repositorio

### Requisitos
- **Bash** (v5+)
- **R** (4.1)

### Environments
- rnaseq_preprocessing.yml
- r_env_4.1.0.yml
- rsem_1.3.0.yml

### Instalación con Conda
``` bash 
conda env create -f rnaseq_preprocessing.yml
conda activate rnaseq_preprocessing
