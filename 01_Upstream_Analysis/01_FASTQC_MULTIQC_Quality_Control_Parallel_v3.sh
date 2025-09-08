#!/bin/bash

# Script para realizar el análisis de calidad en secuencias de RNA-seq de Illumina
# con FastQC, fastq_screen y MultiQC.

# Función para mostrar el uso correcto del script
mostrar_uso() {
    echo "Uso: $0 -i <directorio_de_entrada> -o <directorio_de_salida> [-t <hilos>]"
    echo "  -i: Directorio de entrada donde se encuentran los archivos .fastq"
    echo "  -o: Directorio de salida donde se generarán los reportes"
    echo "  -t: (Opcional) Número de hilos a usar. Por defecto 8."
    exit 1
}

# Valores por defecto
threads=8  # Número de hilos por defecto

# Procesar las opciones con getopts
while getopts "i:o:t:" opt; do
  case $opt in
    i) input_dir="$OPTARG" ;;
    o) output_dir="$OPTARG" ;;
    t) threads="$OPTARG" ;;
    *) mostrar_uso ;;
  esac
done

# Verificar que se hayan proporcionado los argumentos obligatorios
if [ -z "$input_dir" ] || [ -z "$output_dir" ]; then
    echo "Error: Debes proporcionar los directorios de entrada y salida."
    mostrar_uso
fi

# Verificar que el directorio de entrada exista
if [ ! -d "$input_dir" ]; then
    echo "El directorio de entrada $input_dir no existe."
    exit 1
fi

# Crear un directorio general para los resultados si no existe
mkdir -p "$output_dir"

# Crear subdirectorios para los resultados de FastQC y fastq_screen
fastqc_dir="$output_dir/fastqc_results"
fastq_screen_dir="$output_dir/fastq_screen_results"
mkdir -p "$fastqc_dir" "$fastq_screen_dir"

# Crear un archivo para registrar los archivos con errores
error_log="$output_dir/archivos_con_error.txt"
> "$error_log"  # Vaciar el archivo en caso de que ya exista

echo "Iniciando análisis de calidad con FastQC utilizando $threads hilos..."

# Función para ejecutar FastQC y manejar errores
process_fastq_file() {
    fastq_file="$1"
    echo "Procesando $fastq_file con FastQC..."
    fastqc -o "$fastqc_dir" "$fastq_file" -t "$threads"  # Usar el número de hilos especificado

    if [ $? -ne 0 ]; then
        echo "Error procesando $fastq_file con FastQC. Continuando con el siguiente archivo..."
        echo "$fastq_file" >> "$error_log"  # Registrar el archivo con error en el log
    else
        echo "Finalizado: $fastq_file"
    fi
}

# Exportar la función y variables necesarias para GNU parallel
export -f process_fastq_file
export fastqc_dir
export error_log
export threads

# Usar GNU parallel para procesar múltiples archivos en paralelo con progreso
find "$input_dir" -name "*.fastq*" | parallel -j "$threads" --progress --joblog "$output_dir/parallel_joblog.txt" process_fastq_file {}

echo "Análisis de FastQC finalizado."

# Ejecución de fastq_screen de manera secuencial para cada archivo .fastq
echo "Iniciando análisis de calidad con fastq_screen..."

for fastq_file in "$input_dir"/*.fastq*; do
    echo "Procesando $fastq_file con fastq_screen..."
    fastq_screen --outdir "$fastq_screen_dir" "$fastq_file"

    if [ $? -ne 0 ]; then
        echo "Error procesando $fastq_file con fastq_screen. Continuando con el siguiente archivo..."
        echo "$fastq_file" >> "$error_log"  # Registrar el archivo con error en el log
    else
        echo "Finalizado: $fastq_file"
    fi
done

echo "Análisis de fastq_screen finalizado."

# Crear un reporte consolidado con MultiQC dentro del directorio de salida
echo "Generando reporte consolidado con MultiQC..."
multiqc "$output_dir" -o "$output_dir"

if [ $? -eq 0 ]; then
    echo "Reporte de MultiQC generado exitosamente en $output_dir."
else
    echo "Error generando el reporte de MultiQC"
    exit 1
fi

echo "Proceso de control de calidad completado."

# Guardar el código del script actual en un archivo de texto dentro del directorio de salida
script_name=$(basename "$0")
cp "$0" "$output_dir/script_utilizado.txt"

echo "El script utilizado ha sido guardado."
echo "Los archivos que generaron error se han registrado en $error_log."

