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


if [ ! -z "${SLURM_CPUS_PER_TASK}" ] ; then
    export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
    export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}
else
    export OMP_NUM_THREADS=1
    export SRUN_CPUS_PER_TASK=1
fi

trim_galore --paired --length 20 -o trimmed DE34NGSUKBR151851_1.fq.gz DE34NGSUKBR151851_2.fq.gz

