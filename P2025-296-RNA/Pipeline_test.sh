#!/bin/bash
#SBATCH --account=ag_ukb_irnb_bruestle
#SBATCH --job-name=Pipeline
#SBATCH --output=Pipeline_%j.log    # Protokoll-Datei
#SBATCH --error=Pipeline_%j.err     # Fehler-Protokoll
#SBATCH --time=01:00:00                 # Maximale Laufzeit
#SBATCH --partition=intelsr_devel
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jaicher@uni-bonn.de  # Tipp: E-Mail-Benachrichtigung
#SBATCH --nodes=1 # nodes requested
#SBATCH --tasks-per-node=4 # number of tasks per node
#SBATCH --cpus-per-task=24 # number of allocated cores per process

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
    echo "Start TrimmGalore"
    trim_galore --paired --length 20 -o trimmed $1 $2 

    echo "Start Mixed Allignment"
    hisat2 -x "${Workspacepath}/indexes/combined_index" -1 "trimmed/$3_1_val_1.fq.gz" -2 "trimmed/$3_2_val_2.fq.gz" -S "combined_sams/$3_combined.sam" --summary-file "combined_sams/$3_alignment_summary.txt"

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
    hisat2 -p 8 -x "${Workspacepath}/indexes/human_index" -1 "separated_files/$3_1_human.fq.gz" -2 "separated_files/$3_2_human.fq.gz" -S "realligned_files/$3_human_final.sam"
    echo "Human reads alligned"
    echo "Start Allignment of Human Samples"
    # Mouse Reads mit Mouse Genome alignieren
    hisat2 -p 8 -x "${Workspacepath}/indexes/mouse_index" -1 "separated_files/$3_1_mouse.fq.gz" -2 "separated_files/$3_2_mouse.fq.gz"  -S "realligned_files/$3_mouse_final.sam"
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
    
    srun --exclusive --nodes 1 --ntasks 1 --output "log_err/${SAMPLE}.log" --error "log_err/${SAMPLE}.err" pipeline_single_sample ${R1} ${R2} ${SAMPLE} &
    echo "${SAMPLE} pipeline gestarted."

done
wait
echo "Alle Samples verarbeitet!"