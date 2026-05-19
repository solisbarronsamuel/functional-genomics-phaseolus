#!/bin/bash

# ==========================================
# Directorios
# ==========================================

SRA_DIR="/home/samuel/Documents/school/final-project/functional-genomics-phaseolus/01_RawData/sra/16S"

FASTQ_DIR="/home/samuel/Documents/school/final-project/functional-genomics-phaseolus/01_RawData/fastq/16S"

LIST="/home/samuel/Documents/school/final-project/functional-genomics-phaseolus/01_RawData/txt/16S/SRR_Acc_List16S.txt"

# ==========================================
# Crear directorios si no existen
# ==========================================

mkdir -p "$SRA_DIR"
mkdir -p "$FASTQ_DIR"

# ==========================================
# Descargar y convertir
# ==========================================

while read -r SRR
do

    echo "Procesando: $SRR"

    # ----------------------------------
    # Descargar .sra
    # ----------------------------------

    prefetch "$SRR" \
        --output-directory "$SRA_DIR"

    # ----------------------------------
    # Convertir a FASTQ
    # ----------------------------------

    fasterq-dump "$SRA_DIR/$SRR/$SRR.sra" \
        --outdir "$FASTQ_DIR" \
        --threads 8

done < "$LIST"

echo "DESCARGA COMPLETADA"

 