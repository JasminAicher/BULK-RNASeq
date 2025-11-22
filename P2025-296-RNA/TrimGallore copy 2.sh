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

export PATH=/home/jaicher_hpc/TrimGalore-0.6.10:/opt/software/easybuild-INTEL/software/Python/3.11.3-GCCcore-12.3.0/bin:/opt/software/easybuild-INTEL/software/OpenSSL/1.1/bin:/opt/software/easybuild-INTEL/software/XZ/5.4.2-GCCcore-12.3.0/bin:/opt/software/easybuild-INTEL/software/SQLite/3.42.0-GCCcore-12.3.0/bin:/opt/software/easybuild-INTEL/software/Tcl/8.6.13-GCCcore-12.3.0/bin:/opt/software/easybuild-INTEL/software/ncurses/6.4-GCCcore-12.3.0/bin:/opt/software/easybuild-INTEL/software/bzip2/1.0.8-GCCcore-12.3.0/bin:/opt/software/easybuild-INTEL/software/binutils/2.40-GCCcore-12.3.0/bin:/opt/software/easybuild-INTEL/software/GCCcore/12.3.0/bin:/home/jaicher_hpc/.local/bin:/home/jaicher_hpc/bin:/home/jaicher_hpc:/home/jaicher_hpc/subread-2.0.2-Linux-x86_64/bin:/home/jaicher_hpc/samtools-1.20-install/bin:/home/jaicher_hpc/hisat2-2.2.1:/opt/software/easybuild-INTEL/software/Miniforge3/24.1.2-0/bin:/opt/software/easybuild-INTEL/software/Miniforge3/24.1.2-0/condabin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin
export LD_LIBRARY_PATH=/opt/software/easybuild-INTEL/software/Python/3.11.3-GCCcore-12.3.0/lib:/opt/software/easybuild-INTEL/software/OpenSSL/1.1/lib:/opt/software/easybuild-INTEL/software/libffi/3.4.4-GCCcore-12.3.0/lib64:/opt/software/easybuild-INTEL/software/XZ/5.4.2-GCCcore-12.3.0/lib:/opt/software/easybuild-INTEL/software/SQLite/3.42.0-GCCcore-12.3.0/lib:/opt/software/easybuild-INTEL/software/Tcl/8.6.13-GCCcore-12.3.0/lib:/opt/software/easybuild-INTEL/software/libreadline/8.2-GCCcore-12.3.0/lib:/opt/software/easybuild-INTEL/software/ncurses/6.4-GCCcore-12.3.0/lib:/opt/software/easybuild-INTEL/software/bzip2/1.0.8-GCCcore-12.3.0/lib:/opt/software/easybuild-INTEL/software/binutils/2.40-GCCcore-12.3.0/lib:/opt/software/easybuild-INTEL/software/zlib/1.2.13-GCCcore-12.3.0/lib:/opt/software/easybuild-INTEL/software/GCCcore/12.3.0/lib64


if [ ! -z "${SLURM_CPUS_PER_TASK}" ] ; then
    export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
    export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}
else
    export OMP_NUM_THREADS=1
    export SRUN_CPUS_PER_TASK=1
fi

echo "Start byPrim"
cat "combined_sams/DE82NGSUKBR151860_combined.sam" | python byPrim.py -s DE82NGSUKBR151860


