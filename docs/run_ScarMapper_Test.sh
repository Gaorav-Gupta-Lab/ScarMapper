#!/bin/bash
#Parameter file to run ScarMapper Pipeline

python3 /mnt/hgfs/OneDrive_UNC/Projects/ScarMapper/scarmapper.py --options_file /mnt/hgfs/OneDrive_UNC/Projects/ScarMapper/docs/run_ScarMapper_Test.sh
exit

--IndelProcessing	True

--FASTQ1	/mnt/hgfs/OneDrive_UNC/Projects/ScarMapper/FASTQ/ScarMapper_Test_Trim.R1.fastq.gz
--FASTQ2	/mnt/hgfs/OneDrive_UNC/Projects/ScarMapper/FASTQ/ScarMapper_Test_Trim.R2.fastq.gz
#--FASTQ1	/mnt/hgfs/OneDrive_UNC/Projects/ScarMapper/FASTQ/RJ_JCG1_S1_L001_R1_001.fastq.gz
#--FASTQ2	/mnt/hgfs/OneDrive_UNC/Projects/ScarMapper/FASTQ/RJ_JCG1_S1_L001_R2_001.fastq.gz

--Ref_Seq	/mnt/hgfs/OneDrive/Bioinformatics/RefSeq/GRCm38/GRCm38.p6.fa.gz
--Master_Index_File	/mnt/hgfs/OneDrive/Bioinformatics/Indices/Master_Lists/Ramsden_Rosa26a_Oligos.bed
--Index_File	/mnt/hgfs/OneDrive_UNC/Projects/ScarMapper/docs/Ramsden_Index3_TEST.bed
--Target_File	/mnt/hgfs/OneDrive_UNC/Projects/ScarMapper/docs/ScarMapper_Targets.txt

--Working_Folder	/mnt/hgfs/Drive_D/ScarMapper_Testing/

--Verbose	DEBUG
--Job_Name	Testing
--Spawn	3
--Demultiplex	False
--Species	Mouse
--Platform	Ramsden

--Atropos_Trim	False
--Anchored_Adapters_5p	/mnt/hgfs/OneDrive/Bioinformatics/Atropos_Adapters/5'_anchored_adapters.fa
--Anchored_Adapters_3p	/mnt/hgfs/OneDrive/Bioinformatics/Atropos_Adapters/3'_anchored_adapters.fa
--Atropos_Aligner	adapter
--NextSeq_Trim	1
--Adapter_Mismatch_Fraction	0.15
--Read_Queue_Size	500000
--Result_Queue_Size	100000

--N_Limit	0.01
--Minimum_Length	100	# Length after trimming
--OutputRawData	False

