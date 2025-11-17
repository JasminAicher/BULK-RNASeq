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
#SBATCH --partition=intelsr_short

# Dateipräfix- und Zielordner definieren
OUTPUT_DIR=../Trimmomatic
module load Java

# Automatisch alle R1-Dateien finden und Sample-Namen extrahieren
for R1 in *_1.fq.gz; do
    # Sample-Namen aus Dateinamen extrahieren
    SAMPLE=$(basename ${R1} _1.fq.gz)
    
    echo "Verarbeite ${SAMPLE}..."
    
    # Input-Dateien
    R2=${SAMPLE}_2.fq.gz

    # Prüfen ob R2 existiert
    if [ ! -f ${R2} ]; then
        echo "WARNUNG: ${R2} nicht gefunden. Überspringe ${SAMPLE}."
        continue
    fi
    
    # Trimmomatic ausführen
    java -jar trimmomatic-0.40.jar PE ${R1} ${R2} ${OUTPUT_DIR}/${SAMPLE}_1_paired.fq.gz ${OUTPUT_DIR}/${SAMPLE}_1_unpaired.fq.gz ${OUTPUT_DIR}/${SAMPLE}_2_paired.fq.gz ${OUTPUT_DIR}/${SAMPLE}_2_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
    
    echo "${SAMPLE} abgeschlossen."
done

echo "Alle Samples verarbeitet!"


