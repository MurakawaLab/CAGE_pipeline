#!/usr/bin/env bash

### This is a log file to obtain STAR Index for RNA-seq and CAGE analysis. 
### March 23 2022

### Contact: Murakawa lab at RIKEN Yokohama and Kyoto University
### Written by:		Shruti Bhagat; bhagat.shruti.53p@st.kyoto-u.ac.jp; Raku Son; raku.son@riken.jp 
### Reviewed by:	Akiko Oguchi; akiko.oguchi@riken.jp 
###					Kazuhiro Takeuchi; takeuchi.kazuhiro.45v@st.kyoto-u.ac.jp 


### GENCODE v 39
### Primary assembly (PRI) 
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh38.primary_assembly.genome.fa.gz
### md5sum from gencode 6e8f6b77ec0c31e6ccf0d0856619d209  GRCh38.primary_assembly.genome.fa.gz
md5sum 6e8f6b77ec0c31e6ccf0d0856619d209  GRCh38.primary_assembly.genome.fa.gz

### Comprehensive PRI gtf file
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/gencode.v39.primary_assembly.annotation.gtf.gz

### Get Human ribosomal DNA complete repeating unit https://www.ncbi.nlm.nih.gov/nuccore/U13369.1?report=fasta
### 18S rRNA; 28S rRNA; 5.8S rRNA; 5'ETS; 3'ETS; ITS1; ITS2; intergenic spacer; cdc27 pseudogene; p53 binding site.

### Get Human 5S ribosomal DNA https://www.ncbi.nlm.nih.gov/nuccore/V00589.1?report=fasta

### Concatenate as U13369.1_5S_human_rRNA_220131.fa

### Make ribosomal index
STAR --runThreadN 12 \
 --runMode genomeGenerate \
 --genomeSAindexNbases 6 \
 --genomeDir STARindex_U13369.1_5S_rDNA \
 --genomeFastaFiles U13369.1_5S_human_rRNA_220131.fa

### Make GRCh38 PRI GENCODE v 39 index
STAR --runThreadN 12 \
 --runMode genomeGenerate \
 --genomeDir STARindex_GENCODEv39_GRCh38PRI \
 --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
 --sjdbGTFfile gencode.v39.primary_assembly.annotation.gtf \
 --sjdbOverhang 149

