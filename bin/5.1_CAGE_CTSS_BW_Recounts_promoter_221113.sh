#!/usr/bin/env bash

### prep
function usage()
{
  cat <<EOF
 November 13th 2022
 This script can be used to gain promoter counts from CTSS Bigwig files. 
 
 This script uses bigWigAverageOverBed (http://hgdownload.soe.ucsc.edu/admin/exe/).
 Please edit your PATH environmental variable in bash_profile so that the above tools can be executable from any directory.
 
 And /path/to/files/ need be provided for the requirements given below:
 -p A bed 6 file containing promoters.

 All input raw CTSS Bigwig files must be in the working directory and separated to have the suffix ".fwd.bw" for the forward strand and ".rev.bw" for the reverse strand.
 
 The major output files are as follows:
 PREFIX.fwd.rev.promoter.recount.txt
 
usage: $0 -p promoter.bed 

Contact: Murakawa lab at RIKEN Yokohama and Kyoto University
Written by:	Shruti Bhagat; bhagat.shruti.53p@st.kyoto-u.ac.jp
		Raku Son; raku.son@riken.jp
Reviewed by:	Akiko Oguchi; akiko.oguchi@riken.jp 
		Kazuhiro Takeuchi; takeuchi.kazuhiro.45v@st.kyoto-u.ac.jp 
		Yasuhiro Murakawa; murakawa.yasuhiro.0r@kyoto-u.ac.jp 
EOF
  exit 1;
}


while getopts p: opt
do
  case ${opt} in
  p) promoter=${OPTARG};;
  *) usage;;
  esac
done

if [ "${promoter}" = "" ]; then usage; fi

### Counting CAGE TSS for promoter on the forward strand
for infile in *.fwd.bw
 do
   base=$(basename ${infile} .fwd.bw)
   ### Promoter forward
   cat ${promoter} \
    | awk '{if($6 == "+"){print}}' \
    | bigWigAverageOverBed ${infile} /dev/stdin /dev/stdout \
    | cut -f 1,4 >> ${base}.promoter.recount.fwd.tab
  done

### Counting CAGE TSS for promoter on the reverse strand
for infile in *.rev.bw
 do
   base=$(basename ${infile} .rev.bw)
   ### Promoter reverse
   cat ${promoter} \
    | awk '{if($6 == "-"){print}}' \
    | bigWigAverageOverBed ${infile} /dev/stdin /dev/stdout \
    | cut -f 1,4 >> ${base}.promoter.recount.rev.tab
 done

### Concatenate forward and reverse counts for promoters
for infile in *.promoter.recount.fwd.tab
 do
   base=$(basename ${infile} .promoter.recount.fwd.tab)
   cat ${infile} ${base}.promoter.recount.rev.tab >> ${base}.promoter.recount.fwd.rev.txt
 done

### Delete unnecessary files
rm *.promoter.recount.fwd.tab
rm *.promoter.recount.rev.tab

### Arrange directories
mkdir Promoter_recounts
mv *.promoter.recount.fwd.rev.txt Promoter_recounts
