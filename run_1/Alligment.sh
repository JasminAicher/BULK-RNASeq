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


OUTPUTDIR=jaicher_hpc@marvin.hpc.uni-bonn.de:/lustre/scratch/data/jaicher_hpc-CNTNAP2/Trimmomatic/sams
INDEX=jaicher_hpc@marvin.hpc.uni-bonn.de:/lustre/scratch/data/jaicher_hpc-CNTNAP2/BULK-RNASeq/run_1

# Schleife über alle *_R1_paired.fq.gz Dateien 
for R1 in "$INPUTDIR"/*_1_paired.fq.gz; do
    SAMPLE=$(basename "$R1" _1_paired.fq.gz)
    R2="$INPUTDIR/${SAMPLE}_2_paired.fq.gz"
    
    # Prüfe ob die R2-Datei existiert
    if [ ! -f "$R2" ]; then
        echo "WARNUNG: $R2 fehlt, überspringe $SAMPLE"
        continue
    fi
    
    echo "Starte Alignment für $SAMPLE"
    hisat2 -p 8 \
      -x "$INDEX" \
      -1 "$R1" \
      -2 "$R2" \
      -S "$OUTPUTDIR/${SAMPLE}_combined.sam" \
      --summary-file "$OUTPUTDIR/${SAMPLE}_alignment_summary.txt"
    echo "Fertig mit $SAMPLE"
    
done