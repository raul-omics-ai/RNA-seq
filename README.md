# 🧬 RNA-seq Analysis

Este repositorio contiene un **pipeline completo para el análisis de datos RNA-seq**, 
desde la preparación de las lecturas crudas hasta la interpretación biológica de los resultados y muchos gráficos chulos.  

El objetivo es proporcionar un flujo **reproducible, automatizado y modular**, 
organizado en tres grandes bloques:

<pre> ``` RNA-seq-analysis/ 
  │── Upstream Analysis/ # Scripts en bash para preprocesamiento y mapeo 
  │ ├── 01_qc.sh 
  │ ├── 02_filtering.sh 
  │ ├── 03_trimming.sh 
  │ └── 04_mapping_counts.sh 
  │── Differential Expression Analysis/ # Script en R para análisis de expresión diferencial 
  │ └── DESeq2_pipeline.R 
  │── Downstream Analysis/ # Scripts en R para análisis biológico 
  │ ├── enrichment_overrepresentation.R 
  │ ├── pathway_enrichment.R 
  │ └── topology_based_enrichment.R 
  └── README.md ``` </pre>

---
## 🚀 Cómo usar este repositorio

### Requisitos
- **Bash** (v5+)
- **R** (4.1)

### Environments
- rnaseq_preprocessing.yml
- r_env_4.1.0.yml
- rsem_1.3.0.yml
