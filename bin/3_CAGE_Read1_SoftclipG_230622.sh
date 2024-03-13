#!/usr/bin/env bash

### prep
function usage()
{
  cat <<EOF
 July 12th 2022
 This script filters CAGE sequencing bam files by only retainnig reads starting with a soft-clipped G.
 
 This script uses Samtools-1.14 or higher (http://www.htslib.org/).
 
 The requirement is as follows:
 -t Thread count.
All input BAM files (only Read1 if you have used paired-end CAGE) must be in the working directory and have the suffix ".bam".
 
usage: $0 -t integer 

Contact: Murakawa lab at RIKEN Yokohama and Kyoto University
Written by:	Akiko Oguchi; akiko.oguchi@riken.jp 
                Shruti Bhagat; bhagat.shruti.53p@st.kyoto-u.ac.jp
                Raku Son; raku.son@riken.jp 
Reviewed by:   	Kazuhiro Takeuchi; takeuchi.kazuhiro.45v@st.kyoto-u.ac.jp 
		Yasuhiro Murakawa; murakawa.yasuhiro.0r@kyoto-u.ac.jp 
EOF
  exit 1;
}


while getopts G:t:p:e: opt
do
  case ${opt} in
  t) IntegerT=${OPTARG};;
  *) usage;;
  esac
done

BASE=G
ALT=C

for infile in *.bam
 do
  base=$(basename ${infile} .bam)
  samtools view -@ ${IntegerT} -H ${infile} > ${base}_header.sam

### Gain reads with one-nucleotide soft-clipped on the forward strand
  samtools view -@ ${IntegerT} -F 16 ${infile} | awk -F '\t' -v BASE=${BASE} '
  BEGIN {OFS="\t"} {
  if ($6 ~ /^1S[0-9]/ && $10 ~ /^'${BASE}'/) {print $0} \
  }
 ' \
  > ${base}_SoftclipG_F.sam

 ### Gain reads with one-nucleotide soft-clipped on the reverse strand
  samtools view -@ ${IntegerT} -f 16 ${infile} | awk -F '\t' -v ALT=${ALT} '
  BEGIN {OFS="\t"} {
  if ($6 ~ /[0-9]M1S$/ && $10 ~ /'${ALT}'$/) {print $0} \
    }
   ' \
  > ${base}_SoftclipG_R.sam

### Combine the header, F.sam and R.sam
cat \
  ${base}_header.sam \
  ${base}_SoftclipG_F.sam \
  ${base}_SoftclipG_R.sam \
 | samtools sort -@ ${IntegerT} -O bam -o SoftclipG_${base}.bam

### Delete unnecessary files
rm ${base}*.sam
done
