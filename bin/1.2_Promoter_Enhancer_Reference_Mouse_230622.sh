#!/usr/bin/env bash

### This is a log file to obtain mouse chromosome sizes and FANTOM5 promoters and enhancers. 
### May 31st 2023

### Contact: Murakawa lab at RIKEN Yokohama and Kyoto University
### Written by:		Raku Son; raku.son@riken.jp; Shruti Bhagat; bhagat.shruti.53p@st.kyoto-u.ac.jp 
### Reviewed by:	Akiko Oguchi; akiko.oguchi@riken.jp 
###					Kazuhiro Takeuchi; takeuchi.kazuhiro.45v@st.kyoto-u.ac.jp 


### Mouse GRCm38 chromosome sizes
### Chromosome size file can be produced from the reference fasta by samtools (http://www.htslib.org/)
### An example of creating chromosome size file from GENCODE reference file is shown below
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz
gunzip GRCm38.primary_assembly.genome.fa.gz
samtools faidx GRCm38.primary_assembly.genome.fa
cut -f1,2 GRCm38.primary_assembly.genome.fa.fai > GRCm38.primary_assembly.chrom.sizes

### If you use references from ucsc reference, it can also be downloaded from https://hgdownload.cse.ucsc.edu/goldenpath/
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes
sort -k1,1 -k2,2n mm10.chrom.sizes > mm10.chrom.sorted.sizes

### Mouse mm10 FANTOM5 promoters
wget https://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/extra/CAGE_peaks/mm10_fair+new_CAGE_peaks_phase1and2.bed.gz
### Mouse mm10 FANTOM5 enhancers
wget https://fantom.gsc.riken.jp/5/datafiles/reprocessed/mm10_latest/extra/enhancer/F5.mm10.enhancers.bed.gz
### Mouse mm9 FANTOM-NET enhancers
wget https://fantom.gsc.riken.jp/5/suppl/Hirabayashi_et_al_2019/data/Supplementary_Data_2_Mouse_FANTOM-NET_enhancers.bed.gz

## Unzip
gunzip *.gz

### Get number of promoters and enhancers
wc -l *bed

### 164672 mm10_fair+new_CAGE_peaks_phase1and2.bed
### 49797 F5.mm10.enhancers.bed
### 57646 Supplementary_Data_2_Mouse_FANTOM-NET_enhancers.bed

### Convert bed12 to bed6 by selecting first 6 columns 
for infile in *.bed
 do
 echo ${infile}
   base=$(basename ${infile} .bed)
       awk 'BEGIN{OFS="\t"}{ print $1, $2, $3, $4, "0", $6 }' ${infile} | sort -k1,1 -k2,2n \
       > ${base}.6.bed
 done

### Edit promoter name to include hg38 coordinates
awk 'BEGIN{OFS="\t"}{ print $1, $2, $3, $1":"$2".."$3","$6"|"$4, $5, $6 }' mm10_fair+new_CAGE_peaks_phase1and2.6.bed > mm10_fair+new_CAGE_peaks_phase1and2.edit.6.bed

### Get number of promoters and enhancers
wc -l *6.bed
### 164672 mm10_fair+new_CAGE_peaks_phase1and2.6.bed
### 49797 F5.mm10.enhancers.6.bed
### 57646 Supplementary_Data_2_Mouse_FANTOM-NET_enhancers.6.bed

### Convert mm9 to mm10 coordinates with default parameters for FANTOM-NET enhancers
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz
gunzip *.gz
liftOver Supplementary_Data_2_Mouse_FANTOM-NET_enhancers.6.bed mm9ToMm10.over.chain Supplementary_Data_2_Mouse_FANTOM-NET_enhancers.mm10.bed Supplementary_Data_2_Mouse_FANTOM-NET_enhancers.6.unmapped.bed

wc -l Supplementary_Data_2_Mouse_FANTOM-NET_enhancers.mm10.bed 
### 57642
wc -l Supplementary_Data_2_Mouse_FANTOM-NET_enhancers.6.unmapped.bed 
### 8

awk 'BEGIN{OFS="\t"}{ print $1, $2, $3, $1":"$2".."$3"|mm9::"$4, $5, $6 }' Supplementary_Data_2_Mouse_FANTOM-NET_enhancers.mm10.bed > Supplementary_Data_2_Mouse_FANTOM-NET_enhancers.mm10.6.bed
