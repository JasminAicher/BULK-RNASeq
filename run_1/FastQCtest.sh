#!/bin/bash
#SBATCH --account=ag_ukb_irnb_bruestle
#SBATCH --job-name=hisat2_index
#SBATCH --output=hisat2_index_%j.log    # Protokoll-Datei
#SBATCH --error=hisat2_index_%j.err     # Fehler-Protokoll
#SBATCH --time=8:00:00                 # Maximale Laufzeit
#SBATCH --partition=intelsr_short
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jaicher@uni-bonn.de  # Tipp: E-Mail-Benachrichtigung
#SBATCH --nodes=1 # nodes requested
#SBATCH --tasks-per-node=1 # number of tasks per node
#SBATCH --cpus-per-task=25 # number of allocated cores per process
#SBATCH --mem=250G

fastqc -t 8 -o qc_results DE07NGSUKBR151852_1_paired.fq.gz DE07NGSUKBR151852_2_paired.fq.gz
