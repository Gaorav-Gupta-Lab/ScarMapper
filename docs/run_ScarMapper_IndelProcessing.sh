#!/bin/bash
#Parameter file to run ScarMapper Pipeline

python3 /full/path/to/scarmapper.py --options_file /full/path/to/run_ScarMapper_IndelProcessins.sh
exit

--IndelProcessing	True

--FASTQ1	/full/path/to/FASTQ1.gz
--FASTQ2	/full/path/to/FASTQ2.gz

--RefSeq	/full/path/to/reference_sequence.fa
--Master_Index_File	/full/path/to/Master_Index_File
--SampleManifest	/full/path/to/SampleManifest.txt
--TargetFile	/full/path/to/Target_File.txt

--Working_Folder	/full/path/to/Working_Folder/

--Verbose	INFO
--Job_Name	# No spaces or special characters
--Spawn	3
--Demultiplex	False
--Species	# Mouse or Human
--Platform	# Illumina, Ion, Ramsden

--N_Limit	0.01
--Minimum_Length	100	# Length after trimming
--OutputRawData	False # True or False.  Output raw data files.

# This is to trim adapter sequences.  Will probably be removed.
--Atropos_Trim	False
--Anchored_Adapters_5p	/full/path/to/5'_anchored_adapters.fa
--Anchored_Adapters_3p	/full/path/to/3'_anchored_adapters.fa
--Atropos_Aligner	adapter
--NextSeq_Trim	1
--Adapter_Mismatch_Fraction	0.15
--Read_Queue_Size	500000
--Result_Queue_Size	100000
