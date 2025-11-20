#!/bin/bash
#SBATCH --account=ag_ukb_irnb_bruestle
#SBATCH --job-name=FeatureCount
#SBATCH --output=FeatureCount_%j.log    # Protokoll-Datei
#SBATCH --error=FeatureCount_%j.err     # Fehler-Protokoll
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


#counts human reads
featureCounts -T 8 -p -B -a Homo_sapiens.GRCh38.113.gtf -o counts_human.txt -g gene_id -t exon sample_human_final.sorted.bam
echo "human reads sorted"
#counts mouse reads
featureCounts -T 8 -p -B -a Mus_musculus.GRCm39.113.gtf -o counts_mouse.txt -g gene_id -t exon sample_mouse_final.sorted.bam
echo "mouse reads sorted"
echo "done"