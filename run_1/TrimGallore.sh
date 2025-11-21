#!/bin/bash
#SBATCH --account=ag_ukb_irnb_bruestle
#SBATCH --job-name=Gallore
#SBATCH --output=Gallore_%j.log    # Protokoll-Datei
#SBATCH --error=Gallore_%j.err     # Fehler-Protokoll
#SBATCH --time=8:00:00                 # Maximale Laufzeit
#SBATCH --partition=intelsr_short
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jaicher@uni-bonn.de  # Tipp: E-Mail-Benachrichtigung
#SBATCH --nodes=1 # nodes requested
#SBATCH --tasks-per-node=1 # number of tasks per node
#SBATCH --cpus-per-task=25 # number of allocated cores per process
#SBATCH --mem=250G

module load Python/3.11.3-GCCcore-12.3.0


INPUTDIR=/lustre/scratch/data/jaicher_hpc-CNTNAP2/P2025-296-RNA/

if [ ! -z "${SLURM_CPUS_PER_TASK}" ] ; then
    export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
    export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}
else
    export OMP_NUM_THREADS=1
    export SRUN_CPUS_PER_TASK=1
fi

# Schleife über alle *_R1_paired.fq.gz Dateien 
for R1 in "$INPUTDIR"/*_1_paired.fq.gz; do
    SAMPLE=$(basename "$R1" _1_paired.fq.gz)
    R2="$INPUTDIR/${SAMPLE}_2_paired.fq.gz"
    
    # Prüfe ob die R2-Datei existiert
    if [ ! -f "$R2" ]; then
        echo "WARNUNG: $R2 fehlt, überspringe $SAMPLE"
        continue
    fi
    
 trim_galore --paired --fastqc --length 20 -o trimmed DE07NGSUKBR151852_1_paired.fq.gz DE07NGSUKBR151852_2_paired.fq.gz
    
done

