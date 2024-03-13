#!/usr/bin/env bash

### prep
function usage()
{
  cat <<EOF
 November 13th 2022
 This script can be used to gain enhancer counts from CTSS Bigwig files. 
 
 This script uses bigWigAverageOverBed (http://hgdownload.soe.ucsc.edu/admin/exe/).
 Please edit your PATH environmental variable in bash_profile so that the above tools can be executable from any directory.
 
 And /path/to/files/ need be provided for the requirements given below:
 -e A bed 6 file containing enhancers with 6th column of "."
 
 All input raw CTSS Bigwig files must be in the working directory and separated to have the suffix ".fwd.bw" for the forward strand and ".rev.bw" for the reverse strand.
 
 The output files are as follows:
 PREFIX.fwd.enhancer.recount.txt
 PREFIX.rev.enhancer.recount.txt
 
usage: $0 -e enhancer.bed 

Contact: Murakawa lab at RIKEN Yokohama and Kyoto University
Written by:	Shruti Bhagat; bhagat.shruti.53p@st.kyoto-u.ac.jp
		Raku Son; raku.son@riken.jp
Reviewed by:	Akiko Oguchi; akiko.oguchi@riken.jp 
		Kazuhiro Takeuchi; takeuchi.kazuhiro.45v@st.kyoto-u.ac.jp 
		Yasuhiro Murakawa; murakawa.yasuhiro.0r@kyoto-u.ac.jp 
EOF
  exit 1;
}

while getopts e: opt
do
  case ${opt} in
  e) enhancer=${OPTARG};;
  *) usage;;
  esac
done

if [ "${enhancer}" = "" ]; then usage; fi

### Counting CAGE TSS for enhancer on the forward strand
for infile in *.fwd.bw
 do
   base=$(basename ${infile} .fwd.bw)
    ### Enhancer forward
    cat ${enhancer} \
    | bigWigAverageOverBed ${infile} /dev/stdin /dev/stdout \
    | cut -f 1,4 >> ${base}.enhancer.recount.fwd.txt
 done

### Counting CAGE TSS for enhancer on the reverse strand
for infile in *.rev.bw
 do
   base=$(basename ${infile} .rev.bw)
    ### Enhancer reverse
    cat ${enhancer} \
    | bigWigAverageOverBed ${infile} /dev/stdin /dev/stdout \
    | cut -f 1,4 >> ${base}.enhancer.recount.rev.txt
 done

### Arrange directories
mkdir Enhancer_recounts
mv *.enhancer.recount.fwd.txt Enhancer_recounts
mv *.enhancer.recount.rev.txt Enhancer_recounts
