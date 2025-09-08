#!/bin/bash

# Script para procesar múltiples archivos FASTQ (paired-end y single-end) con fastp utilizando getopts

# Variables por defecto
input_dir=""
output_dir=""
file_type=""
CUT_FRONT=26                    # Cortar por el extremo 5' las bases con un PhredScore < 25
CUT_TAIL=26                     # Cortar por el extremo 3' las bases con un PhredScore < 25
CUT_MEAN_QUALITY=26             # Cortar según el valor de calidad Phred medio
DETECT_ADAPTER_FOR_PE="ENABLE"  # Detectar los adaptadores y cortarlos
TRIM_POLY_G="ENABLE"            # Por defecto para NovaSeq/NextSeq
TRIM_POLY_X="ENABLE"            # Por defecto para eliminar colas de polyA de mRNA-seq data
MIN_LENGTH=75                   # Longitud mínima de las reads
THREADS=16                       # Número de hilos a utilizar (por defecto, todos los disponibles)

# Función de uso para mostrar cómo utilizar el script
usage() {
    echo "Uso: $0 -i <directorio_de_entrada> -o <directorio_de_salida> -p <tipo_de_archivo> [-f CUT_FRONT] [-t CUT_TAIL] [-m CUT_MEAN_QUALITY] [-l MIN_LENGTH] [-h THREADS]"
    echo "  -i: Directorio de entrada donde se encuentran los archivos .fastq"
    echo "  -o: Directorio de salida donde se generarán los reportes"
    echo "  -p: Tipo de archivo: 'paired' para paired-end o 'single' para single-end"
    echo "  -f: PhredScore mínimo para cortar por el extremo 5'"
    echo "  -t: PhredScore mínimo para cortar por el extremo 3'"
    echo "  -m: Cortar según el valor de calidad Phred medio"
    echo "  -l: Longitud mínima de las lecturas"
    echo "  -h: (Opcional) Número de hilos a usar. Por defecto 16"
    exit 1
}

# Parsear las opciones utilizando getopts
while getopts ":i:o:p:f:t:m:l:h:" opt; do
    case ${opt} in
        i) input_dir="$OPTARG" ;;
        o) output_dir="$OPTARG" ;;
        p) file_type="$OPTARG" ;;
        f) CUT_FRONT="$OPTARG" ;;
        g) CUT_TAIL="$OPTARG" ;;
        m) CUT_MEAN_QUALITY="$OPTARG" ;;
        l) MIN_LENGTH="$OPTARG" ;;
        h) THREADS="$OPTARG" ;;
        \?) echo "Opción inválida: -$OPTARG" >&2; usage ;;
        :) echo "La opción -$OPTARG requiere un argumento." >&2; usage ;;
    esac
done

# Verifica que las rutas de input, output y tipo de archivo hayan sido proporcionadas
if [ -z "$input_dir" ] || [ -z "$output_dir" ] || [ -z "$file_type" ]; then
    echo "Error: Debes proporcionar los directorios de entrada (-i), salida (-o) y el tipo de archivo (-p)."
    usage
fi

# Validar tipo de archivo
if [[ "$file_type" != "paired" && "$file_type" != "single" ]]; then
    echo "Error: El tipo de archivo debe ser 'paired' o 'single'."
    usage
fi

# Crear directorios de salida
clean_fastq_dir="$output_dir/01_Clean_fastq/"
qc_postprocessing_dir="$output_dir/02_fastp_QC"
mkdir -p "$clean_fastq_dir" "$qc_postprocessing_dir"

# Guardar el script actual en un archivo .txt dentro del directorio de salida
cp "$0" "$output_dir/script_utilizado.txt"

# Guardar las variables configurables en un archivo .txt
cat <<EOL > "$output_dir/variables_utilizadas.txt"
Variables utilizadas en fastp:
ADAPTER TRIMMING OPTIONS:
 --detect_adapter_for_pe: By default, the adapter sequence auto-detection is enabled for SE data only, turn on this option to enable it for PE data ($DETECT_ADAPTER_FOR_PE).

POLY-G TAIL TRIMMING:
 --trim_poly_g: Force polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data ($TRIM_POLY_G).

POLY-X TAIL TRIMMING:
 --trim_poly_x: PolyA trimming in 3' ends ($TRIM_POLY_X).

PER READ CUTTING BY QUALITY OPTIONS:
 --cut_front: Move a sliding window from front (5') to tail, drop the bases in the window if its mean quality < $CUT_FRONT.
 --cut_tail: Move a sliding window from tail (3') to front, drop the bases in the window if its mean quality < $CUT_TAIL.
 --cut_mean_quality: Drop the bases if its mean quality < $CUT_MEAN_QUALITY.
LENGTH FILTERING OPTIONS:
 --length_required: Reads shorter than $MIN_LENGTH will be discarded.

THREADING OPTIONS:
 --thread: $THREADS worker thread number.
EOL

echo "El script utilizado ha sido guardado en $output_dir/script_utilizado.txt"
echo "Las variables configurables han sido guardadas en $output_dir/variables_utilizadas.txt"

########## FASTP FILTRADO Y RECORTE ##########
if [ "$file_type" == "paired" ]; then
    for R1_file in "$input_dir"/*_R1_001.fastq.gz; do
        # Crea el nombre correspondiente del archivo R2 reemplazando R1 con R2
        R2_file="${R1_file/_R1_/_R2_}"

        # Comprueba que el archivo R2 existe
        if [[ -f "$R2_file" ]]; then
            # Extrae el prefijo del nombre del archivo para usarlo en la salida
            base_name=$(basename "$R1_file" _R1_001.fastq.gz)

            echo "Procesando la muestra (paired-end): $base_name"
            # Ejecuta fastp con los archivos R1 y R2
            fastp -i "$R1_file" -I "$R2_file" \
                  -o "$clean_fastq_dir/${base_name}_filtered_R1.fastq.gz" \
                  -O "$clean_fastq_dir/${base_name}_filtered_R2.fastq.gz" \
                  --cut_front "$CUT_FRONT" \
                  --cut_tail "$CUT_TAIL" \
                  --cut_mean_quality "$CUT_MEAN_QUALITY" \
                  --detect_adapter_for_pe \
                  --trim_poly_g \
                  --trim_poly_x \
                  --length_required "$MIN_LENGTH" \
                  --thread "$THREADS" \
                  --html "$qc_postprocessing_dir/${base_name}_fastp_report.html" \
                  --json "$qc_postprocessing_dir/${base_name}_fastp_report.json"
        else
            echo "Archivo R2 no encontrado para: $R1_file"
        fi
    done
elif [ "$file_type" == "single" ]; then
    for SE_file in "$input_dir"/*_R1*.fastq.gz; do
        if [[ -f "$SE_file" ]]; then
            base_name=$(basename "$SE_file" _R1_001.fastq.gz)

            echo "Procesando la muestra (single-end): $base_name"
            # Ejecuta fastp con el archivo single-end
            fastp -i "$SE_file" \
                  -o "$clean_fastq_dir/${base_name}_filtered.fastq.gz" \
                  --cut_front "$CUT_FRONT" \
		  --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
                  --cut_tail "$CUT_TAIL" \
                  --cut_mean_quality "$CUT_MEAN_QUALITY" \
                  --trim_poly_g \
                  --trim_poly_x \
                  --length_required "$MIN_LENGTH" \
                  --thread "$THREADS" \
                  --html "$qc_postprocessing_dir/${base_name}_fastp_report.html" \
                  --json "$qc_postprocessing_dir/${base_name}_fastp_report.json"
        else
            echo "Archivo single-end no encontrado: $SE_file"
        fi
    done
fi

echo "Procesamiento con fastp completado"

echo "Comenzando con el procesamiento de la calidad de los archivos limpios"

./01_FASTQC_MULTIQC_Quality_Control_Parallel_v3.sh -i "$clean_fastq_dir" -o "$output_dir/03_MultiQC_QC"

