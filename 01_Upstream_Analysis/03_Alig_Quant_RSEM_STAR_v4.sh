#!/bin/bash
START_TIME=$(date +%s)

usage() {
    echo "Uso: $0 [-g GTF_FILE] [-r RSEM_REF] [-i FASTQ_DIR] [-f FASTQ_LIST_FILE] [-o OUTPUT_DIR] [-t THREADS]"
    echo "  -g  Ruta al archivo GTF"
    echo "  -r  Ruta al directorio de referencia de RSEM"
    echo "  -i  Directorio que contiene los archivos FASTQ (*.fastq.gz)"
    echo "  -f  Archivo de texto con rutas absolutas a archivos *_R1.fastq.gz"
    echo "  -o  Directorio de salida"
    echo "  -t  Número de hilos a usar (por defecto 32)"
    echo "  -h  Mostrar esta ayuda"
    exit 1
}

while getopts "g:r:i:f:o:t:h" opt; do
  case $opt in
    g) GTF_FILE="$OPTARG" ;;
    r) RSEM_REF="$OPTARG" ;;
    i) FASTQ_DIR="$OPTARG" ;;
    f) FASTQ_LIST_FILE="$OPTARG" ;;
    o) OUTPUT_DIR="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    h) usage ;;
    *) usage ;;
  esac
done

# Validación básica
if [[ -z "$FASTQ_DIR" && -z "$FASTQ_LIST_FILE" ]]; then
    echo "Error: Debe especificar un directorio de FASTQ (-i) o un archivo con rutas (-f)"
    usage
fi

THREADS=${THREADS:-32}
MEMORY="64G"
ulimit -n 4096

mkdir -p "$OUTPUT_DIR"
BAM_DIR="${OUTPUT_DIR}/01_BAM_FILES"
mkdir -p "${BAM_DIR}"
RSEM_OUTPUT_DIR="${OUTPUT_DIR}/02_RSEM_OUTS"
mkdir -p "${RSEM_OUTPUT_DIR}"
QUALIMAP_REPORT_DIR="${OUTPUT_DIR}/03_Qualimap_report"
mkdir -p "${QUALIMAP_REPORT_DIR}"
COUNT_MATRIX_DIR="${OUTPUT_DIR}/04_Count_Matrices"
mkdir -p "${COUNT_MATRIX_DIR}"
TMP_DIR="/media/rfernandez/Windows/tmp_mapping"
mkdir -p "${TMP_DIR}"

SCRIPT_PATH="${OUTPUT_DIR}/Mapping_quantification_script_used.txt"
cat "$0" > "$SCRIPT_PATH"
MULTI_BAMQC_FILE="${QUALIMAP_REPORT_DIR}/Sample_bam_paths.txt"
ERROR_LOG="${OUTPUT_DIR}/errores.log"
INPUT_FILE_LIST="${OUTPUT_DIR}/input_fastq_list.txt"

# Generar lista de archivos FASTQ
if [[ -n "$FASTQ_LIST_FILE" ]]; then
    echo "Usando archivo de entrada proporcionado: $FASTQ_LIST_FILE"
    cp "$FASTQ_LIST_FILE" "$INPUT_FILE_LIST"
else
    echo "Generando lista de archivos FASTQ en ${INPUT_FILE_LIST}"
    find "${FASTQ_DIR}" -type f -name "*_R1.fastq.gz" | sort > "${INPUT_FILE_LIST}"
fi

# Procesamiento
while IFS= read -r R1; do
    SAMPLE=$(basename "${R1}" _R1.fastq.gz)
    R2="${R1/_R1.fastq.gz/_R2.fastq.gz}"

    if [ ! -f "${R2}" ]; then
        echo "Advertencia: Archivo R2 no encontrado para la muestra ${SAMPLE}. Saltando esta muestra." >> "$ERROR_LOG"
        continue
    fi

    echo "Procesando la muestra: ${SAMPLE}"

    rsem-calculate-expression \
        -p ${THREADS} \
        --star \
        --star-path ~/miniconda3/envs/rnaseq_preprocessing/bin/ \
        --star-gzipped-read-file \
        --temporary-folder "${TMP_DIR}" \
        --paired-end "${R1}" "${R2}" \
        "${RSEM_REF}" \
        "${RSEM_OUTPUT_DIR}/${SAMPLE}_RSEM"

    if [ $? -ne 0 ]; then
        echo "Error: RSEM falló para la muestra ${SAMPLE}. Saltando esta muestra." >> "$ERROR_LOG"
        continue
    fi

    echo "Moviendo BAMs y matrices de conteo para ${SAMPLE}"
    #mv "${RSEM_OUTPUT_DIR}/${SAMPLE}_RSEM.genome.sorted.bam" "${BAM_DIR}/" 2>/dev/null
    #mv "${RSEM_OUTPUT_DIR}/${SAMPLE}_RSEM.genome.sorted.bam.bai" "${BAM_DIR}/" 2>/dev/null
    mv "${RSEM_OUTPUT_DIR}/${SAMPLE}_RSEM.genes.results" "${COUNT_MATRIX_DIR}/" 2>/dev/null

    #echo "Ejecutando Qualimap para la muestra ${SAMPLE}"
    #qualimap bamqc \
    #    -bam "${BAM_DIR}/${SAMPLE}_RSEM.genome.sorted.bam" \
    #    -c \
    #    -gd MOUSE \
    #    -gff "${GTF_FILE}" \
    #    -nt ${THREADS} \
    #    --java-mem-size=${MEMORY} \
    #    -p strand-specific-forward \
    #    -outdir "${QUALIMAP_REPORT_DIR}/${SAMPLE}_qualimap_report"

    if [ $? -ne 0 ]; then
        echo "Error: Qualimap falló para la muestra ${SAMPLE}. Saltando esta muestra." >> "$ERROR_LOG"
    else
        echo "Qualimap finalizado correctamente para ${SAMPLE}"
    fi

done < "${INPUT_FILE_LIST}"

# Multi-bamqc
#echo "Generando archivo ${MULTI_BAMQC_FILE}"
#for sample_dir in "${QUALIMAP_REPORT_DIR}"/*_qualimap_report/; do
#    sample_name=$(basename "$sample_dir" "_qualimap_report")
#    echo -e "${sample_name}\t${sample_dir}" >> "${MULTI_BAMQC_FILE}"
#done

#echo "Ejecutando qualimap multi-bamqc"
#qualimap multi-bamqc --data "${MULTI_BAMQC_FILE}"
#if [ $? -ne 0 ]; then
#    echo "Error: Qualimap Multi-bamqc falló." >> "$ERROR_LOG"
#else
#    echo "Qualimap Multi-bamqc finalizado correctamente."
#fi

# Limpieza
echo "Eliminando carpeta temporal: ${TMP_DIR}"
rm -rf "${TMP_DIR}"

# Resumen
echo "Generando resumen final..."
TOTAL_MUESTRAS=$(wc -l < "${INPUT_FILE_LIST}")
ERRORES=0
if [ -f "${ERROR_LOG}" ]; then
    ERRORES=$(grep -c "Error" "${ERROR_LOG}")
fi
PROCESADAS=$((TOTAL_MUESTRAS - ERRORES))

RESUMEN="${OUTPUT_DIR}/resumen_final.txt"
{
    echo "=============================="
    echo "     RESUMEN DEL ANÁLISIS"
    echo "=============================="
    echo "Total de muestras detectadas : ${TOTAL_MUESTRAS}"
    echo "Muestras procesadas con éxito: ${PROCESADAS}"
    echo "Muestras con errores         : ${ERRORES}"
    echo ""

    if [ -f "${ERROR_LOG}" ]; then
        echo "==== MUESTRAS CON ERRORES ===="
        grep "Error" "${ERROR_LOG}" | sed 's/Error: //' | sort | uniq
    else
        echo "No se detectaron errores."
    fi

    echo ""
    echo "==== TIEMPO DE EJECUCIÓN ===="
    printf "Duración total: %02d:%02d:%02d (hh:mm:ss)\n" $((DURATION / 3600)) $(((DURATION % 3600) / 60)) $((DURATION % 60))
} > "${RESUMEN}"

echo "Resumen final guardado en: ${RESUMEN}"
echo "Script finalizado. Revisa ${ERROR_LOG} para ver muestras con errores."

