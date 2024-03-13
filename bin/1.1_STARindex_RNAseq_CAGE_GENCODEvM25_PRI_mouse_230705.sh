#!/usr/bin/env bash

### This is a log file to obtain STAR Index for RNA-seq and CAGE analysis. 
### July 5 2023

### Contact: Murakawa lab at RIKEN Yokohama and Kyoto University
### Written by:		Shruti Bhagat; bhagat.shruti.53p@st.kyoto-u.ac.jp; Raku Son; raku.son@riken.jp 
### Reviewed by:	Akiko Oguchi; akiko.oguchi@riken.jp 
 
### mouse version

### GENCODE v M25
### Primary assembly (PRI) 
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz
### md5sum from gencode 3bc591be24b77f710b6ba5d41022fc5a  GRCm38.primary_assembly.genome.fa.gz
md5sum 3bc591be24b77f710b6ba5d41022fc5a  GRCm38.primary_assembly.genome.fa.gz

### Comprehensive PRI gtf file
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.primary_assembly.annotation.gtf.gz
### md5sum from gencode c5125258a0a2c5250ddb4c192abbf4e8  gencode.vM25.primary_assembly.annotation.gtf.gz
md5sum c5125258a0a2c5250ddb4c192abbf4e8  gencode.vM25.primary_assembly.annotation.gtf.gz

### Get Mouse ribosomal DNA complete repeating unit https://www.ncbi.nlm.nih.gov/nuccore/BK000964.3?report=fasta
### Get Mouse 5S ribosomal RNA https://www.ncbi.nlm.nih.gov/nuccore/NR_030686.1?report=fasta

### Concatenate as BK000964.3_5S_mouse_rRNA_221002.fa

### Make ribosomal index
STAR --runThreadN 12 \ 
  --runMode genomeGenerate \
  --genomeSAindexNbases 6 \
  --genomeDir STARindex_BK000964.3_5S_rDNA_mouse \
  --genomeFastaFiles BK000964.3_5S_mouse_rRNA_221002.fa 

### Make GRCm38 PRI GENCODE vM25 index
STAR --runThreadN 12 \
 --runMode genomeGenerate \
 --genomeDir GRCm38_index \
 --genomeFastaFiles GRCm38.primary_assembly.genome.fa \
 --sjdbGTFfile gencode.vM25.primary_assembly.annotation.gtf \
 --sjdbOverhang 149

