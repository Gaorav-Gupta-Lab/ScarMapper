#!/bin/bash
#Parameter file to run ScarMapper Pipeline
#Version 1

python3.10 /mnt/hgfs/OneDrive_UNC/Projects/Programs/ScarMapper/scarmapper.py --options_file /mnt/hgfs/Drive_D/run_ScarMapper_IndelProcessins.sh
exit

--IndelProcessing	True

--FASTQ1	/mnt/hgfs/Drive_D/C226a-Hek293T-siBRCA2-C_R1_001.fastq.gz
--FASTQ2	/mnt/hgfs/Drive_D/C226a-Hek293T-siBRCA2-C_R2_001.fastq.gz

--RefSeq	/mnt/hgfs/OneDrive/Bioinformatics/RefSeq/GRCh38/GRCh38.p12.fa.bgz
--Master_Index_File	/mnt/hgfs/Drive_D/MSK_Indices.txt
--SampleManifest	/mnt/hgfs/Drive_D/MSK_SampleManifest.txt
--TargetFile	/mnt/hgfs/Drive_D/ScarMapper_Targets.txt
--WorkingFolder	/mnt/hgfs/Drive_D/

--Verbose	INFO# INFO or DEBUG
--Job_Name	MSK_Test# No spaces or special characters, prepended to output files.
--Spawn	3 # How many parallel jobs?  Max should be n-1 threads or cpu's.  Minimum is 1.
--Demultiplex	False# True or False.  Write demultiplexed FASTQ files?
--DeleteConsensusFASTQ	True
--HR_Donor	# 10 - 15 nucleotide sequence for HR Donor search.  Can be left blank. 
--Platform	Illumina# Illumina, TruSeq, Ramsden

--N_Limit	0.01
--Minimum_Length	100	# Length after trimming
--OutputRawData	False # True or False.  Output raw data files.

# PEAR Options.  Leave blank for defaults
--TestMethod	
--PValue	0.05 # Default 0.01
--Memory	20000M # Default 200M.  Recomend >1000M
--MinOverlap	# Default 10
--QualityThreshold	# Default 40
--PhredValue	# Default 33
--MinConsensusLength	# Default 50

# Plot Options
--PatternThreshold	0.0001# Cutoff frequency for patterns to plot such as 0.0001
--FigureType	pdf # svg, jpg, tiff, pdf, png