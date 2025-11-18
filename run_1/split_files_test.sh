#!/bin/bash
#SBATCH --account=ag_ukb_irnb_bruestle
#SBATCH --job-name=split_files_test
#SBATCH --output=split_files_test_%j.log    # Protokoll-Datei
#SBATCH --error=split_files_test_%j.err     # Fehler-Protokoll
#SBATCH --time=8:00:00                 # Maximale Laufzeit
#SBATCH --partition=intelsr_short
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jaicher@uni-bonn.de  # Tipp: E-Mail-Benachrichtigung
#SBATCH --nodes=1 # nodes requested
#SBATCH --tasks-per-node=1 # number of tasks per node
#SBATCH --cpus-per-task=25 # number of allocated cores per process
#SBATCH --mem=250G

# Variablen anpassen
TABLE="DE07NGSUKBR151852_table.tsv"
FASTQ1="DE07NGSUKBR151852_1_paired.fq.gz"
FASTQ2="DE07NGSUKBR151852_2_paired.fq.gz"
SAMPLE="DE07NGSUKBR151852"

echo "Extrahiere Human Read-IDs"
awk '$4=="human"{print $1}' table.tsv > human_ids.txt

echo "Extrahiere Mouse Read-IDs"
awk '$4=="mouse"{print $1}' table.tsv > mouse_ids.txt


echo "Erzeuge Human FASTQ .gz Dateien"
seqtk subseq "$FASTQ1" human_ids.txt | gzip > "${SAMPLE}_1_human.fq.gz"
seqtk subseq "$FASTQ2" human_ids.txt | gzip > "${SAMPLE}_2_human.fq.gz"

echo "Erzeuge Mouse FASTQ .gz Dateien"
seqtk subseq "$FASTQ1" mouse_ids.txt | gzip > "${SAMPLE}_1_mouse.fq.gz"
seqtk subseq "$FASTQ2" mouse_ids.txt | gzip > "${SAMPLE}_2_mouse.fq.gz"

echo "Fertig! Ergebnis sind direkt komprimierte Files:
- ${SAMPLE}_1_human.fq.gz
- ${SAMPLE}_2_human.fq.gz
- ${SAMPLE}_1_mouse.fq.gz
- ${SAMPLE}_2_mouse.fq.gz
"
