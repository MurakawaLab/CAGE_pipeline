# CAGE_pipeline
A set of scripts to analyse CAGE datatsets.


Most of these scripts are based on the pipeline developed by the FANTOM5 consortium. 
Some key references: Forrest et al Nature 2014; Andersson et al Nature 2014

All scripts are in the bin directory, please check the scripts for detailed usage. 
A brief description is provided below.

1.1: Creating STAR index for human and mouse. This script uses GENCODE gtf files.
1.2: Creating a mask file used during enhancer identification. 
2: Mapping single-end CAGE reads.
3: Identifying reads with an unencodedG at the 5' end. 
4: Identifying CAGE transcription start sites (CTSS). 
5: Counting CAGE reads mapping to known promoters and enhancers.
6: Identifying bidirectionally transcribed enhancers (Original code: https://github.com/anderssonrobin/enhancers)


