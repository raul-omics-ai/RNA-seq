#!/bin/bash

# Función para mostrar la ayuda
mostrar_ayuda() {
    echo "Uso: $0 [-o OUTPUT_DIR] [-l SRR_FILE] [-c compress] <-h help>"
    echo "  -o OUTPUT_DIR    Ruta donde se guardarán los archivos .fastq"
    echo "  -l SRR_FILE      Ruta al archivo que contiene la lista de SRRs"
    echo "  -c compress      Comprime los archivos FASTQ usando pigz (solo si se llama)"
    echo "  -h help"
    exit 1
}

# Parseo de opciones usando getoptse
while getopts ":o:l:c:h:" opt; do
    case ${opt} in
        o ) OUT_DIR="$OPTARG" ;;
        l ) SRR_FILE="$OPTARG" ;;
        c ) COMPRESS=true ;;
        h ) mostrar_ayuda ;;
        \? ) echo "Opción inválida: -$OPTARG" 1>&2
             mostrar_ayuda ;;
        : ) echo "Opción sin argumento: -$OPTARG" 1>&2
            mostrar_ayuda ;;
    esac
done
shift $((OPTIND -1))

# Crea los directorios de salida si no existen
mkdir -p "$OUT_DIR"

# Verifica si el archivo de SRR existe
if [ ! -f "$SRR_FILE" ]; then
    echo "El archivo de SRR no se encontró: $SRR_FILE"
    exit 1
fi

# Itera sobre cada SRR en el archivo
for SRR in $(cat "$SRR_FILE"); do
    echo "Descargando $SRR..."

    # Descargar archivo SRA usando prefetch
    prefetch -O "$OUT_DIR" "$SRR" --max-size 50g

    # Verifica si el archivo SRA fue descargado correctamente
    SRA_FILE_PATH="${OUT_DIR}/${SRR}.sra"
    if [ ! -f "$SRA_FILE_PATH" ]; then
        SRA_SUBFOLDER_PATH="${OUT_DIR}/${SRR}/${SRR}.sra"
        if [ ! -f "$SRA_SUBFOLDER_PATH" ]; then
            echo "Error: el archivo SRA para $SRR no se encontró. Saltando..."
            continue
        else
            SRA_FILE_PATH="$SRA_SUBFOLDER_PATH"
        fi
    fi

    echo "Convirtiendo ${SRR}.sra a formato FASTQ..."

    # Convierte a FASTQ usando fasterq-dump con 32 hilos
    fasterq-dump "$SRA_FILE_PATH" --split-files -e 30 -p -O "$OUT_DIR" -x -t "$OUT_DIR"

    # Verifica si fasterq-dump fue exitoso
    if [ $? -ne 0 ]; then
        echo "Error al convertir $SRR con fasterq-dump. Saltando..."
        continue
    fi

    # Determinar si es paired-end o single-end
    if [ -f "${OUT_DIR}/${SRR}_2.fastq" ]; then
        # Caso paired-end
        echo "Detección: Paired-end para $SRR"
        
        if [ "$COMPRESS" = true ]; then
            echo "Comprimiendo archivos FASTQ paired-end de $SRR con pigz..."
            # Comprimir tanto _1.fastq como _2.fastq
            pigz -p 15 "${OUT_DIR}/${SRR}_1.fastq" "${OUT_DIR}/${SRR}_2.fastq"

            # Verifica si la compresión fue exitosa y elimina los archivos originales
            if [ $? -eq 0 ]; then
                echo "Eliminando archivos FASTQ originales no comprimidos para paired-end..."
                rm -f "${OUT_DIR}/${SRR}_1.fastq" "${OUT_DIR}/${SRR}_2.fastq"
            else
                echo "Error durante la compresión de los archivos FASTQ paired-end de $SRR."
            fi
        fi

    else
        # Caso single-end
        echo "Detección: Single-end para $SRR"
        
        if [ "$COMPRESS" = true ]; then
            echo "Comprimiendo archivo FASTQ single-end de $SRR con pigz..."
            # Comprimir solo el archivo _1.fastq
            pigz -p 15 "${OUT_DIR}/${SRR}_1.fastq"

            # Verifica si la compresión fue exitosa y elimina el archivo original
            if [ $? -eq 0 ]; then
                echo "Eliminando archivo FASTQ original no comprimido para single-end..."
                rm -f "${OUT_DIR}/${SRR}_1.fastq"
            else
                echo "Error durante la compresión del archivo FASTQ single-end de $SRR."
            fi
        fi
    fi

    # Limpieza: elimina la carpeta SRR descargada
    echo "Eliminando subcarpetas SRR..."
    if [ -d "${OUT_DIR}/${SRR}" ]; then
        rm -r "${OUT_DIR}/${SRR}"
    fi

    # Elimina el SRR descargado del archivo de lista
    sed -i "/^${SRR}$/d" "$SRR_FILE"

    echo "Proceso completado para $SRR."
done

echo "Proceso completado para todos los SRR en $SRR_FILE."

