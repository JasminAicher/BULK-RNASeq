#!/bin/bash
#SBATCH --account=ag_ukb_irnb_bruestle
#SBATCH --job-name=byPrim
#SBATCH --output=byPrim_%j.log    # Protokoll-Datei
#SBATCH --error=byPrim_%j.err     # Fehler-Protokoll
#SBATCH --time=8:00:00                 # Maximale Laufzeit
#SBATCH --partition=intelsr_short
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jaicher@uni-bonn.de  # Tipp: E-Mail-Benachrichtigung
#SBATCH --nodes=1 # nodes requested
#SBATCH --tasks-per-node=1 # number of tasks per node
#SBATCH --cpus-per-task=25 # number of allocated cores per process
#SBATCH --mem=250G

module load Python/3.9

export PATH=$HOME/hisat2-2.2.1:$PATH
export PATH=$HOME/samtools-1.20:$PATH
export PATH=$HOME/optimized/byPrim.py:$PATH

python $HOME/optimized/byPrim.py -s sample_combined.sam -r1 DE07NGSUKBR151852_1_paired.fq.gz -r2 DE07NGSUKBR151852_2_paired.fq.gz -prefix DE07NGSUKBR151852
