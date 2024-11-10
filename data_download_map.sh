#!/bin/bash
################################################################################
#                               RNA-seq analysis                               #
#                                                                              #
# Downloads 6 samples from Sequence Read Archive (3 controls, 3 treated).      #
# Consecutively extracts data to FASTQ files, performs subsampling to 10%      #
#   of pair-end reads and produces quality control reports.                    #
# Creates a STAR index for human chromosome X, maps the subsamples to it       #
#   and counts number of unique pairs mapping to each gene.                    #
#                                                                              #
################################################################################

# from R project directory

##### Data download and subsampling ############################################

base_url="https://sra-pub-run-odp.s3.amazonaws.com/sra"
declare -a samples=("SRR5223500" "SRR5223522" "SRR5223543" \
	"SRR5223505" "SRR5223547" "SRR5223570")

# SRA/fastq files take a lot of space: so download, process and delete
cd data

for sample in "${samples[@]}"; do

    url="${base_url}/${sample}/${sample}"
    echo ${url}
    
    # download
    cd data
    echo "Downloading ${sample} from ${base_url}..."
    curl -O ${url}

    # extract and split fastq files (compressed)
    fastq-dump -split-files --gzip ${sample}
    rm ${sample}

    # sample 10% of reads
    for file in ${sample}_*; do
    	seqkit sample -s100 -p 0.1 ${file} -o sample_${file}
    done
    rm ${sample}_*

    # quality control of sequence data
    fastqc sample_${sample}_* -o fastqc_reports

done

cd ..

################################################################################


##### STAR index and mapping, gene counts ######################################

# annotation data
ann_fasta="annotations/Homo_sapiens.GRCh38.dna.chromosome.X.fa"
ann_gtf="annotations/Homo_sapiens.GRCh38.chrX.gtf"

# prepare STAR index - for chrX (chromosome with FOXP3)
mkdir indexes
STAR --runThreadN 4 \
     --runMode genomeGenerate --genomeDir indexes/chrX \
     --genomeFastaFiles ${ann_fasta} \
     --sjdbGTFfile ${ann_gtf} \
     --sjdbOverhang 74 --outFileNamePrefix chrX \
     --genomeSAindexNbases 12


for sample in "${samples[@]}"; do

    # mapping and QC
    STAR --runThreadN 4 \
      --genomeDir indexes/chrX \
      --readFilesIn data/sample_${sample}_1.fastq.gz \
      data/sample_${sample}_2.fastq.gz \
      --readFilesCommand gunzip -c \
      --outSAMtype BAM SortedByCoordinate \
      --quantMode GeneCounts \
      --outFileNamePrefix alignments/${sample}
    
done

# gene counts
featureCounts -T 4 -p -a ${ann_gtf} -o alignments/counts.txt \
alignments/*Aligned.sortedByCoord.out.bam

################################################################################


