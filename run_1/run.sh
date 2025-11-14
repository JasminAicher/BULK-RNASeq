#!/bin/bash
#SBATCH --account=ag_ukb_irnb_bruestle
#SBATCH --job-name=hisat2_index
#SBATCH --output=hisat2_index_%j.log    # Protokoll-Datei
#SBATCH --error=hisat2_index_%j.err     # Fehler-Protokoll
#SBATCH --time=08:00:00                 # Maximale Laufzeit
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

hisat2-build -p 16 combined_genome.fa combined_index

wait
echo "All jobs completed."