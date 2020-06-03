#!/bin/bash
#Parameter file to run ScarMapper Pipeline
#File generated 2020-06-03

python3 /full/path/to/scarmapper.py --options_file /full/path/to/run_ScarMapper_IndelProcessins.sh
exit

--IndelProcessing	True

--FASTQ1	/full/path/to/FASTQ1.gz
--FASTQ2	/full/path/to/FASTQ2.gz

--RefSeq	/full/path/to/reference_sequence.fa
--Master_Index_File	/full/path/to/Master_Index_File
--SampleManifest	/full/path/to/SampleManifest.csv
--TargetFile	/full/path/to/Target_File.txt
--PrimerPhasingFile	/full/path/to/PrimerPhasing.csv

--WorkingFolder	/full/path/to/Working_Folder/<Location of output files>

--Verbose	# INFO or DEBUG
--Job_Name	# No spaces or special characters, prepended to output files.
--Spawn	3 # How many parallel jobs?  Max should be n-1 threads or cpu's.  Minimum is 1.
--Demultiplex	# True or False.  Write demultiplexed FASTQ files?
--HR_Donor	# 10 - 15 nucleotide sequence for HR Donor search.  Can be left blank. 
--Platform	# Illumina, Ramsden

--N_Limit	0.01
--Minimum_Length	100	# Length after trimming
--OutputRawData	False # True or False.  Output raw data files.

