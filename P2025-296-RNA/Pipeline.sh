#!/bin/bash
#SBATCH --account=ag_ukb_irnb_bruestle
#SBATCH --job-name=Pipeline
#SBATCH --output=Pipeline_%j.log    # Protokoll-Datei
#SBATCH --error=Pipeline_%j.err     # Fehler-Protokoll
#SBATCH --time=24:00:00                 # Maximale Laufzeit
#SBATCH --partition=intelsr_medium
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jaicher@uni-bonn.de  # Tipp: E-Mail-Benachrichtigung
#SBATCH --nodes=5 # nodes requested
#SBATCH --tasks-per-node=4 # number of tasks per node
#SBATCH --cpus-per-task=24 # number of allocated cores per process

module purge

export PATH=/home/jaicher_hpc/TrimGalore-0.6.10:/opt/software/easybuild-INTEL/software/Python/3.11.3-GCCcore-12.3.0/bin:/opt/software/easybuild-INTEL/software/OpenSSL/1.1/bin:/opt/software/easybuild-INTEL/software/XZ/5.4.2-GCCcore-12.3.0/bin:/opt/software/easybuild-INTEL/software/SQLite/3.42.0-GCCcore-12.3.0/bin:/opt/software/easybuild-INTEL/software/Tcl/8.6.13-GCCcore-12.3.0/bin:/opt/software/easybuild-INTEL/software/ncurses/6.4-GCCcore-12.3.0/bin:/opt/software/easybuild-INTEL/software/bzip2/1.0.8-GCCcore-12.3.0/bin:/opt/software/easybuild-INTEL/software/binutils/2.40-GCCcore-12.3.0/bin:/opt/software/easybuild-INTEL/software/GCCcore/12.3.0/bin:/home/jaicher_hpc/.local/bin:/home/jaicher_hpc/bin:/home/jaicher_hpc:/home/jaicher_hpc/subread-2.0.2-Linux-x86_64/bin:/home/jaicher_hpc/samtools-1.20-install/bin:/home/jaicher_hpc/hisat2-2.2.1:/opt/software/easybuild-INTEL/software/Miniforge3/24.1.2-0/bin:/opt/software/easybuild-INTEL/software/Miniforge3/24.1.2-0/condabin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin
export LD_LIBRARY_PATH=/opt/software/easybuild-INTEL/software/Python/3.11.3-GCCcore-12.3.0/lib:/opt/software/easybuild-INTEL/software/OpenSSL/1.1/lib:/opt/software/easybuild-INTEL/software/libffi/3.4.4-GCCcore-12.3.0/lib64:/opt/software/easybuild-INTEL/software/XZ/5.4.2-GCCcore-12.3.0/lib:/opt/software/easybuild-INTEL/software/SQLite/3.42.0-GCCcore-12.3.0/lib:/opt/software/easybuild-INTEL/software/Tcl/8.6.13-GCCcore-12.3.0/lib:/opt/software/easybuild-INTEL/software/libreadline/8.2-GCCcore-12.3.0/lib:/opt/software/easybuild-INTEL/software/ncurses/6.4-GCCcore-12.3.0/lib:/opt/software/easybuild-INTEL/software/bzip2/1.0.8-GCCcore-12.3.0/lib:/opt/software/easybuild-INTEL/software/binutils/2.40-GCCcore-12.3.0/lib:/opt/software/easybuild-INTEL/software/zlib/1.2.13-GCCcore-12.3.0/lib:/opt/software/easybuild-INTEL/software/GCCcore/12.3.0/lib64

module load Python/3.11.3-GCCcore-12.3.0
module load Biopython/1.84-foss-2024a
module load Java



mkdir combined_sams
mkdir separated_files
mkdir realligned_files
mkdir final_txt
mkdir log_err

Workspacepath=/lustre/scratch/data/jaicher_hpc-CNTNAP2

if [ ! -z "${SLURM_CPUS_PER_TASK}" ] ; then
    export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
    export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}
else
    export OMP_NUM_THREADS=1
    export SRUN_CPUS_PER_TASK=1
fi

pipeline_single_sample() {
    #echo "Start TrimmGalore"
    #trim_galore --paired --length 20 -o trimmed $1 $2 

    echo "Start Mixed Allignment"
    hisat2 -x "${Workspacepath}/indexes/combined_index/combined_index" -1 "trimmed/$3_1_val_1.fq.gz" -2 "trimmed/$3_2_val_2.fq.gz" -S "combined_sams/$3_combined.sam" --summary-file "combined_sams/$3_alignment_summary.txt"

    echo "Start byPrim"
    cat "combined_sams/$3_combined.sam" | python byPrim.py -s $3

    echo "Start Separating Files"
    TABLE="separated_files/$3_table.tsv"
    FASTQ1="trimmed/$3_1_val_1.fq.gz"
    FASTQ2="trimmed/$3_2_val_2.fq.gz"
    echo "Extrahiere Human Read-IDs"
    awk '$4=="human"{print $1}' "$TABLE" > separated_files/human_ids.txt
    echo "Extrahiere Mouse Read-IDs"
    awk '$4=="mouse"{print $1}' "$TABLE" > separated_files/mouse_ids.txt
    echo "Erzeuge Human FASTQ .gz Dateien"
    seqtk subseq "$FASTQ1" separated_files/human_ids.txt | gzip > "separated_files/$3_1_human.fq.gz"
    seqtk subseq "$FASTQ2" separated_files/human_ids.txt | gzip > "separated_files/$3_2_human.fq.gz"
    echo "Erzeuge Mouse FASTQ .gz Dateien"
    seqtk subseq "$FASTQ1" separated_files/mouse_ids.txt | gzip > "separated_files/$3_1_mouse.fq.gz"
    seqtk subseq "$FASTQ2" separated_files/mouse_ids.txt | gzip > "separated_files/$3_2_mouse.fq.gz"

    # Human Reads mit Human Genome alignieren
    echo "Start Allignment of Human Samples"
    hisat2 -p 8 -x "${Workspacepath}/indexes/human_index/human_index" -1 "separated_files/$3_1_human.fq.gz" -2 "separated_files/$3_2_human.fq.gz" -S "realligned_files/$3_human_final.sam"
    echo "Human reads alligned"
    echo "Start Allignment of Human Samples"
    # Mouse Reads mit Mouse Genome alignieren
    hisat2 -p 8 -x "${Workspacepath}/indexes/mouse_index/mouse_index" -1 "separated_files/$3_1_mouse.fq.gz" -2 "separated_files/$3_2_mouse.fq.gz"  -S "realligned_files/$3_mouse_final.sam"
    echo "Mouse reads alligned"
    echo "Start BAM zu SAM"
    # SAM zu BAM konvertieren und sortieren
    samtools view -bS "realligned_files/$3_human_final.sam" | samtools sort -o "realligned_files/$3_human_final.sorted.bam"
    samtools view -bS "realligned_files/$3_mouse_final.sam" | samtools sort -o "realligned_files/$3_mouse_final.sorted.bam"
    echo "Start Index erstellen"
    # Index erstellen
    samtools index "realligned_files/$3_human_final.sorted.bam"
    samtools index "realligned_files/$3_mouse_final.sorted.bam"
    echo "Indexe erstellt"

    #counts human reads
    echo "Start FeatureCountHuman"
    featureCounts -T 8 -p -B -a "${Workspacepath}/genomes/Homo_sapiens.GRCh38.113.gtf" -o "final_txt/$3_counts_human.txt" -g gene_id -t exon "realligned_files/$3_human_final.sorted.bam"
    echo "human reads sorted"
    #counts mouse reads
    echo "Start FeatureCountMouse"
    featureCounts -T 8 -p -B -a "${Workspacepath}/genomes/Mus_musculus.GRCm39.113.gtf" -o "final_txt/$3_counts_mouse.txt" -g gene_id -t exon "realligned_files/$3_mouse_final.sorted.bam"
    echo "mouse reads sorted"
    echo "done"
}

export -f pipeline_single_sample

# Automatisch alle R1-Dateien finden und Sample-Namen extrahieren
for R1 in *_1.fq.gz; do
    # Sample-Namen aus Dateinamen extrahieren
    SAMPLE=$(basename ${R1} _1.fq.gz)
    
    echo "Verarbeite ${SAMPLE}..."
    
    # Input-Dateien
    R2=${SAMPLE}_2.fq.gz

    # Prüfen ob R2 existiert
    if [ ! -f ${R2} ]; then
        echo "WARNUNG: ${R2} nicht gefunden. Überspringe ${SAMPLE}."
        continue
    fi
    
    srun --exclusive --nodes 1 --ntasks 1 --output "log_err/${SAMPLE}.log" --error "log_err/${SAMPLE}.err" bash -c "pipeline_single_sample ${R1} ${R2} ${SAMPLE}" &
    echo "${SAMPLE} pipeline gestarted."

done
wait
echo "Alle Samples verarbeitet!"