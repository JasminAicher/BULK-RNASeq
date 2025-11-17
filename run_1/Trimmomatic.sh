#!/bin/bash
#SBATCH --account=ag_ukb_irnb_bruestle
#SBATCH --job-name=Trimmoatic
#SBATCH --output=Trimmoatic_%j.log    # Protokoll-Datei
#SBATCH --error=Trimmoatic_%j.err     # Fehler-Protokoll
#SBATCH --time=08:00:00                 # Maximale Laufzeit
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jaicher@uni-bonn.de  # Tipp: E-Mail-Benachrichtigung
#SBATCH --nodes=1 # nodes requested
#SBATCH --tasks-per-node=1 # number of tasks per node
#SBATCH --cpus-per-task=25 # number of allocated cores per process
#SBATCH --mem=250G

# Dateipräfix- und Zielordner definieren
INPUT_DIR=jaicher_hpc@marvin.hpc.uni-bonn.de:/lustre/scratch/data/jaicher_hpc-CNTNAP2/P2025-296-RNA/
OUTPUT_DIR=jaicher_hpc@marvin.hpc.uni-bonn.de:/lustre/scratch/data/jaicher_hpc-CNTNAP2/Trimmomatic
ADAPTERS=jaicher_hpc@marvin.hpc.uni-bonn.de:~/adapters/TruSeq3-PE.fa

# Automatisch alle R1-Dateien finden und Sample-Namen extrahieren
for R1 in ${INPUT_DIR}/*_R1.fq.gz; do
    # Sample-Namen aus Dateinamen extrahieren
    SAMPLE=$(basename ${R1} _R1.fq.gz)
    
    echo "Verarbeite ${SAMPLE}..."
    
    # Input-Dateien
    R2=${INPUT_DIR}/${SAMPLE}_R2.fq.gz

    # Prüfen ob R2 existiert
    if [ ! -f ${R2} ]; then
        echo "WARNUNG: ${R2} nicht gefunden. Überspringe ${SAMPLE}."
        continue
    fi
    
    # Trimmomatic ausführen
    trimmomatic PE -threads 4 -phred33 ${R1} ${R2} \
        ${OUTPUT_DIR}/${SAMPLE}_R1_paired.fq.gz ${OUTPUT_DIR}/${SAMPLE}_R1_unpaired.fq.gz \
        ${OUTPUT_DIR}/${SAMPLE}_R2_paired.fq.gz ${OUTPUT_DIR}/${SAMPLE}_R2_unpaired.fq.gz \
        ILLUMINACLIP:${ADAPTERS}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    
    echo "${SAMPLE} abgeschlossen."
done

echo "Alle Samples verarbeitet!"
