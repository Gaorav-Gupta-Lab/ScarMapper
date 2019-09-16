#!/bin/bash
#Parameter file to run ScarMapper Pipeline

# SLURM Commands
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=07-00:00:00
#SBATCH --mem=60g

python3 /full/path/to/scarmapper.py --options_file /full/path/to/run_ScarMapper_IndelProcessins.sh
exit

--IndelProcessing	True

--FASTQ1	/full/path/to/FASTQ1.gz
--FASTQ2	/full/path/to/FASTQ2.gz

--Ref_Seq	/full/path/to/reference_sequence.fa
--Master_Index_File	/full/path/to/Ramsden_Rosa26a_Oligos.bed
--Index_File	/mnt/hgfs/OneDrive_UNC/Projects/ScarMapper/docs/Ramsden_Index3_TEST.bed
--Target_File	/mnt/hgfs/OneDrive_UNC/Projects/ScarMapper/docs/ScarMapper_Targets.txt

--Working_Folder	/full/path/to/Working_Folder/

--Verbose	INFO
--Job_Name	
--Spawn	3
--Demultiplex	False
--Species	Mouse
--Platform	Ramsden

--Atropos_Trim	False
--Anchored_Adapters_5p	/full/path/to/5'_anchored_adapters.fa
--Anchored_Adapters_3p	/full/path/to/3'_anchored_adapters.fa
--Atropos_Aligner	adapter
--NextSeq_Trim	1
--Adapter_Mismatch_Fraction	0.15
--Read_Queue_Size	500000
--Result_Queue_Size	100000

--N_Limit	0.01
--Minimum_Length	100	# Length after trimming
--OutputRawData	False