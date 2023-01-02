#!/bin/bash
#Parameter file to run ScarMapper Pipeline
#Version 2

python3.10 /Full/Path/To/ScarMapper/scarmapper.py --options_file /Full/Path/To/run_ScarMapper_IndelProcessins.sh
exit

--IndelProcessing	True # Always True for Indel Processing

--FASTQ1	/Full/Path/To/R1.fastq.gz # Uncompressed file permitted but not recommended
--FASTQ2	/Full/Path/To/R2.fastq.gz # Uncompressed file permitted but not recommended

--RefSeq	/Full/Path/To/RefSeq/RefSeq.fa.bgz # Index file must be located here as well.
--Master_Index_File	/Full/path/to/index.file
--SampleManifest	/Full/path/to/SampleManifest.file
--TargetFile	/Full/path/to/Target.file
--WorkingFolder	/Full/path/to/working folder/ # This is where the output files will be written.

--Verbose	INFO# INFO or DEBUG
--Job_Name	 # No spaces or special characters, prepended to output files.
--Spawn	3 # How many parallel jobs?  Max should be n-1 threads or cpu's.  Minimum is 1.
--Demultiplex	False # True or False.  Write demultiplexed FASTQ files?
--DeleteConsensusFASTQ	True
--HR_Donor	# 10 - 15 nucleotide sequence for HR Donor search.  Leave blank if not doing HR search.
--Platform	Illumina # Illumina, TruSeq, Ramsden
--Search_KMER	10 # Size of sliding window for junction search.  Recomended size 10.

--N_Limit	0.01
--Minimum_Length	100	# Length after trimming
--OutputRawData	False # True or False.  Output raw data files.

# PEAR Options.  Leave blank for defaults
--TestMethod	
--PValue	0.05 # Default 0.01
--Memory	20000M # Default 200M.  Recommend >3000M.  Program will crash if value > available memory.
--MinOverlap	# Default 10
--QualityThreshold	# Default 40
--PhredValue	# Default 33
--MinConsensusLength	# Default 50

# Plot Options
--PatternThreshold	0.0001 # Cutoff frequency for patterns to plot such as 0.0001
--FigureType	pdf # svg, jpg, tiff, pdf, png