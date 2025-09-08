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
  │ └── DESeq2_pipeline.R 
  │── Downstream Analysis/ # Scripts en R para análisis biológico 
  │ ├── enrichment_overrepresentation.R 
  │ ├── pathway_enrichment.R 
  │ └── topology_based_enrichment.R 
  └── README.md 
  ``` </pre>

---
## Upstream Analysis

En esta carpeta se encuentran los pasos para realizar el control de calidad (QC), filtrado de bases de baja calidad y recorte de adaptadores y mapeo de las secuencias dividido en 3 scripts. 

También hay un script para descargar automáticamente archivos de secuenciación públicos desde el repositorio del Sequence Read Archive (SRA).

---
## Differential Expression Analysis

La función principal es `full_rnaseq_automate_analysis_v4.R` que realiza el control de calidad de las matriz de cuentas, la detección de muestras de baja calidad, el análisis de expresión diferencial (limma-voom y edgeR), así como todas las gráficas para comprobar las asunciones de los modelos.

También están incorporadas las funciones para realizar el análisis de enriquecimiento.

---
## Downstream Analysis

En esta carpeta se encuentran todas las funciones sobre análisis de sobrerrepresentación (over-representation analysis, ORA), análisis de enriquecimiento de conjuntos de genes (gene set enrichment analysis, GSEA) y el análisis de impacto para tomar en cuenta la topología de la ruta de señalización (Signaling Pathway Impact Analysis, SPIA).

Todas estas funciones están implementadas dentro de la función principal `full_rnaeq_automate_analysis_v4.R` pero aquí se encuentran más organizadas.

---

## 🚀 Cómo usar este repositorio

### Requisitos
- **Bash** (v5+)
- **R** (4.1)

### Environments
- rnaseq_preprocessing.yml
- r_env_4.1.0.yml
- rsem_1.3.0.yml
