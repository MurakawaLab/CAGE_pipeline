#!/usr/bin/env bash


### prep
function usage()
{
  cat <<EOF
 April 29th 2022
 This script can be used to process raw reads obtained from single-end CAGE-based sequencing.
 The final output is uniquely mapped bam files which can be visualized using any genome browser tool (like IGV) and used for downstream analysis. 
 
 This script uses the following tools: 
 Trimmomatic-0.39 or higher (used to trim reads. http://www.usadellab.org/cms/?page=trimmomatic). 
 STAR-2.7.10a or higher (used to align fastq files. https://github.com/alexdobin/STAR).
 pigz-2.3.4 or higher (parallel implementation of gzip. https://zlib.net/pigz/).
 Samtools-1.14 or higher (http://www.htslib.org/).
  
 The requirements are as follows:
 -r Path to Directory with STAR index for ribosomal RNA.
 -g Path to Directory with STAR index for reference genome.
 -t Thread count.
 -m Multimap count.
 -p Path to trimmomatic Directory.
All fastq input files must be in your working directory and have the suffix ".fastq.gz".

usage: $0 -r rRNA_STARindex -g Genome_STARindex -t integer -m integer -p Trimmomatic

Contact: Murakawa lab at RIKEN Yokohama and Kyoto University
Written by:		Shruti Bhagat; bhagat.shruti.53p@st.kyoto-u.ac.jp
Reviewed by:	Raku Son; raku.son@riken.jp 
				Akiko Oguchi; akiko.oguchi@riken.jp 
				Kazuhiro Takeuchi; takeuchi.kazuhiro.45v@st.kyoto-u.ac.jp 
				Yasuhiro Murakawa; murakawa.yasuhiro.0r@kyoto-u.ac.jp 
EOF
  exit 1;
}

while getopts r:g:t:m:p: opt
do
  case ${opt} in
  r) rRNA=${OPTARG};;
  g) genome=${OPTARG};;
  t) IntegerT=${OPTARG};;
  m) IntegerM=${OPTARG};;
  p) PathTrimmomatic=${OPTARG};;
  *) usage;;
  esac
done


if [ "${rRNA}" = "" ]; then usage; fi
if [ "${genome}" = "" ]; then usage; fi
if [ "${IntegerT}" = "" ]; then usage; fi
if [ "${IntegerM}" = "" ]; then usage; fi
if [ "${PathTrimmomatic}" = "" ]; then usage; fi


### Trimming with default parameters and minimum length 25 
### If trimmomatic does not work with high number of threads, please input "java -Xms4g -Xmx4g -jar"
for infile in *.fastq.gz
 do
   base=$(basename ${infile} _R1_001.fastq.gz)
   java -jar ${PathTrimmomatic}/trimmomatic-0.39.jar SE -threads ${IntegerT} -phred33 \
                ${infile} ${base}_trimmed.fastq.gz \
                ILLUMINACLIP:${PathTrimmomatic}/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:25
 done

### Arrange directories
mkdir Trimmed
mv *_trimmed.fastq.gz Trimmed
cd Trimmed

### Align with to rDNA genome with default parameters and desired multimap number
for infile in *_trimmed.fastq.gz
 do
   base=$(basename ${infile} _trimmed.fastq.gz)
   STAR \
    --runThreadN ${IntegerT} \
    --genomeDir ${rRNA} \
    --readFilesIn ${infile} \
    --readFilesCommand zcat \
    --outFilterMultimapNmax ${IntegerM} \
    --outReadsUnmapped Fastx \
    --outFilterMismatchNmax 10 \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix ${base}_
 done

### Parallel gzip fastq files
pigz -p ${IntegerT} *.mate1

### Arrange directories
mkdir ForAlignment
mv *.mate1.gz ForAlignment

mkdir MappedRibosomal
mv *.out MappedRibosomal
mv *out* MappedRibosomal
mv *_STARtmp MappedRibosomal

cd ForAlignment

### Align with to genome of interest with default parameters and desired multimap number
for infile in *.mate1.gz
 do
   base=$(basename ${infile} .mate1.gz)
   STAR \
    --runThreadN ${IntegerT} \
    --genomeDir ${genome} \
    --readFilesIn ${infile} \
    --readFilesCommand zcat \
    --outFilterMultimapNmax ${IntegerM} \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${base}_
 done

mkdir BamFiles
mv *.bam BamFiles
cd BamFiles

for infile in *.bam; do
echo ${infile}
base=$(basename ${infile} .bam)
samtools view -@ ${IntegerT} -b -q 255 ${infile} > ${base}.uniqmapped.bam
done

mkdir MappedUniq
mv *.uniqmapped.bam MappedUniq/
cd MappedUniq

