## Analysis of human T-cell RNA-seq
====================================

### DATA DOWNLOAD, QC, MAPPING AND READ COUNTING
To replicate the whole analysis, download the human GTF and FASTA files from  Ensembl [![Ensembl](https://img.shields.io/badge/Ensembl-genome%20data-orange)](http://www.ensembl.org/).  
Subset GTF file for chromosome X and extract to `annotations/Homo_sapiens.GRCh38.chrX.gtf` and save FASTA of chromosome X to `annotations/Homo_sapiens.GRCh38.dna.chromosome.X.fa`.

Samples: run the `data_download_map.sh` for sample download, mapping and generating per gene counts.  
  
*Requirements: `curl`, `fastq-dump`, `seqkit`, `fastqc`, `STAR`, `featureCounts`.*

### DIFFERENTIAL EXPRESSION ANALYSIS
Run the `analysis.R` script.  

*Requirements: `data.table`, `magrittr`, `edgeR` and `ggplot2` packages.*
