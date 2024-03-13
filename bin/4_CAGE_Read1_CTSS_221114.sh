#!/usr/bin/env bash


### prep
function usage()
{
  cat <<EOF
 November 14th 2022
 This script can be used to gain an overview of CAGE sequencing based data.
 It provides ctss bed file, promoter and enhancer counts, and depth nomalized BigWig files for visualization on IGV, UCSC. 
 
 This script uses BEDTools and can be installed from here: https://github.com/arq5x/bedtools2.
 Samtools-1.14 or higher (http://www.htslib.org/).
 bedGraphToBigWig can be downloaded from here: http://hgdownload.soe.ucsc.edu/admin/exe/.
 Please edit your PATH environmental variable in bash_profile so that the above tools can be executable from any directory.
 
 And /path/to/files/ need be provided for the requirements given below:
 -g Sorted chrom sizes file (Human hg38 can be found here: https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/latest/ The UCSC goldenpath also contains referance files for all other species).
 -t Thread count.
 All input BAM files (only Read1 if you have used paired-end CAGE) must be in the working directory and have the suffix ".bam".
 
 The major output files are as follows:
 PREFIX.ctss.bed
 PREFIX.ctss.fwd.cpm.bw
 PREFIX.ctss.rev.cpm.bw
 
usage: $0 -G genome -t integer 

Contact: Murakawa lab at RIKEN Yokohama and Kyoto University
Written by:		Shruti Bhagat; bhagat.shruti.53p@st.kyoto-u.ac.jp
Reviewed by:	Raku Son; raku.son@riken.jp 
				Akiko Oguchi; akiko.oguchi@riken.jp 
				Kazuhiro Takeuchi; takeuchi.kazuhiro.45v@st.kyoto-u.ac.jp 
				Yasuhiro Murakawa; murakawa.yasuhiro.0r@kyoto-u.ac.jp 
EOF
  exit 1;
}


while getopts G:t: opt
do
  case ${opt} in
  G) genome=${OPTARG};;
  t) IntegerT=${OPTARG};;
  *) usage;;
  esac
done


if [ "${genome}" = "" ]; then usage; fi
if [ "${IntegerT}" = "" ]; then usage; fi

### Index BAM file
for infile in *.bam; do
echo ${infile}
samtools index -@ ${IntegerT} ${infile}
done

### Convert bam to bed file
for infile in *.bam
 do
   base=$(basename ${infile} .bam)
       bedtools bamtobed -i ${infile} \
       | awk '$1 ~ /chr/ { print }' | sort -k1,1 -k2,2n \
       > ${base}.bed
 done

### Make bedGraph files and sort them
for infile in *.bed
 do
   base=$(basename ${infile} .bed)   
   ### bedGraph forward
       bedtools genomecov -g ${genome} \
       -i ${base}.bed -5 -bg -strand + \
       | sort -k1,1 -k2,2n \
       > ${base}.fwd.bg
       
   ### bedGraph reverse
       bedtools genomecov -g ${genome} \
       -i ${base}.bed -5 -bg -strand - \
       | sort -k1,1 -k2,2n \
       > ${base}.rev.bg
 done
 
### Convert bedGraph to bigWig to be used for other analysis
for infile in *.bg
 do
   base=$(basename ${infile} .bg)   
   bedGraphToBigWig ${infile} ${genome} ${base}.bw
 done

### Make CTSS bed files
for infile in *.bed
 do
   base=$(basename ${infile} .bed)   
   ### CTSS forward
   awk 'BEGIN{OFS="\t"}{if($6=="+"){print $1,$2,$5}}' ${infile} \
   | sort -k1,1 -k2,2n \
   | bedtools groupby -i stdin -g 1,2 -c 3 -o count \
   | awk 'BEGIN{OFS="\t"}{ print $1, $2, $2+1, $1":"$2".."$2+1",+", $3, "+" }' > ${base}.CTSS.fwd.bed
   
   ### CTSS reverse
   awk 'BEGIN{OFS="\t"}{if($6=="-"){print $1,$3,$5}}' ${infile} \
   | sort -k1,1 -k2,2n \
   | bedtools groupby -i stdin -g 1,2 -c 3 -o count \
   | awk 'BEGIN{OFS="\t"}{ print $1, $2-1, $2, $1":"$2-1".."$2",-", $3, "-" }' > ${base}.CTSS.rev.bed
   
   ### Concatenate and sort
   cat ${base}.CTSS.fwd.bed ${base}.CTSS.rev.bed | sort -k1,1 -k2,2n >> ${base}.CTSS.fwd.rev.bed
 done
 
### CPM transform CTSS bed files and convert to bigWig for visualization
for infile in *.CTSS.fwd.rev.bed
 do
   echo ${infile}
   base=$(basename ${infile} .CTSS.fwd.rev.bed)   
   ### Sum CTSS 
   sum=$(awk 'BEGIN{sum=0}{sum=sum+$5}END{print sum}' ${infile} )
   
   ### CPM transform forward
   awk --assign sum=${sum} 'BEGIN{OFS="\t"} { if($6=="+") { printf("%s\t%i\t%i\t%1.2f\n", $1,$2,$3, 1e6 * $5 / sum) } }' ${infile} \
   | sort -k1,1 -k2,2n > ${base}.CTSS.CPM.fwd.bg
   
   ### CPM transform reverse
   awk --assign sum=${sum} 'BEGIN{OFS="\t"} { if($6=="-") { printf("%s\t%i\t%i\t%1.2f\n", $1,$2,$3, 1e6 * $5 / sum) } }' ${infile} \
   | sort -k1,1 -k2,2n > ${base}.CTSS.CPM.rev.bg
   
   ### Convert to bigWig forward
   bedGraphToBigWig ${base}.CTSS.CPM.fwd.bg ${genome} ${base}.CTSS.CPM.fwd.bw
   
   ### Convert to bigWig reverse
   bedGraphToBigWig ${base}.CTSS.CPM.rev.bg ${genome} ${base}.CTSS.CPM.rev.bw
 done

### Delete unnecessary files
rm *.bg
rm *.CTSS.fwd.bed
rm *.CTSS.rev.bed

### Compress bed files
pigz -p ${IntegerT} *.bed

### Arrange directories
mkdir BigWig
mkdir CTSSbed
mkdir BED
mv *.bw BigWig
mv *.CTSS.fwd.rev.bed.gz CTSSbed
mv *.bed.gz BED
cd BigWig
mkdir RAWbw
mkdir CPMbw
mv *.CPM.* CPMbw
mv *.bw RAWbw
