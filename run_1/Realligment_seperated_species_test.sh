#!/bin/bash
#SBATCH --account=ag_ukb_irnb_bruestle
#SBATCH --job-name=hisat2_realligment
#SBATCH --output=hisat2_realligment_%j.log    # Protokoll-Datei
#SBATCH --error=hisat2_realligment_%j.err     # Fehler-Protokoll
#SBATCH --time=8:00:00                 # Maximale Laufzeit
#SBATCH --partition=intelsr_short
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jaicher@uni-bonn.de  # Tipp: E-Mail-Benachrichtigung
#SBATCH --nodes=1 # nodes requested
#SBATCH --tasks-per-node=1 # number of tasks per node
#SBATCH --cpus-per-task=25 # number of allocated cores per process
#SBATCH --mem=250G

if [ ! -z "${SLURM_CPUS_PER_TASK}" ] ; then
    export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
    export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}
else
    export OMP_NUM_THREADS=1
    export SRUN_CPUS_PER_TASK=1
fi

echo "Running a job"
echo " Using ${SLURM_NNODES} nodes"
echo " Using ${SLURM_NTASKS} tasks"
echo " Using ${OMP_NUM_THREADS} threads per process"
export PATH=$HOME/samtools-1.20:$PATH

# Human Reads mit Human Genome alignieren
#hisat2 -p 8 -x human_index -1 DE07NGSUKBR151852_1_human.fq.gz -2 DE07NGSUKBR151852_2_human.fq.gz -S sample_human_final.sam
#echo "Human reads alligned"
# Mouse Reads mit Mouse Genome alignieren
#hisat2 -p 8 -x mouse_index -1 DE07NGSUKBR151852_1_mouse.fq.gz -2 DE07NGSUKBR151852_2_mouse.fq.gz  -S sample_mouse_final.sam
#echo "Mouse reads alligned"
# SAM zu BAM konvertieren und sortieren
samtools view -bS sample_human_final.sam | samtools sort -o sample_human_final.sorted.bam
samtools view -bS sample_mouse_final.sam | samtools sort -o sample_mouse_final.sorted.bam
echo "BAMzuSAM"
# Index erstellen
samtools index sample_human_final.sorted.bam
samtools index sample_mouse_final.sorted.bam
echo "Indexe erstellt"

wait
echo "All jobs completed.