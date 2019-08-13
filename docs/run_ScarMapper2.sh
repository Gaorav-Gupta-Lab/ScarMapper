#!/bin/bash
#Parameter file to run ScarMapper
#File generated 2019-07-19

# SLURM Commands
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=07-00:00:00
#SBATCH --mem=80g

python3 /nas/longleaf/home/dennis/scripts/ScarMapper/scarmapper.py --options_file /pine/scr/d/e/dennis/run_ScarMapper2.sh
exit

--ScarMapper	True
--ThruPLEX	False

--FASTQ1	/pine/scr/d/e/dennis/ScarMapper_Test_Trim.R1.fastq.gz
--FASTQ2	/pine/scr/d/e/dennis/ScarMapper_Test_Trim.R2.fastq.gz

--Ref_Seq	/proj/guptalab/Reference_Files/GRCm38/GRCm38.p6.fa.gz
--Index_File	/nas/longleaf/home/dennis/Reference_Files/Indicies/ScarMapper/Ramsden_Index2.bed
--Target_File	/nas/longleaf/home/dennis/Reference_Files/Targets/ScarMapper/ScarMapper_Targets.txt
--Working_Folder	/pine/scr/d/e/dennis/Test2/


--Verbose	INFO
--Job_Name	Test2
--Spawn	46
--Demultiplex	False
--Species	Mouse

--Atropos_Trim	False
--Anchored_Adapters_5p	/nas/longleaf/home/dennis/Reference_Files/Atropos_Adapters/5'_anchored_adapters.fa
--Anchored_Adapters_3p	/nas/longleaf/home/dennis/Reference_Files/Atropos_Adapters/3'_anchored_adapters.fa
--Atropos_Aligner	adapter
--NextSeq_Trim	1
--Adapter_Mismatch_Fraction	0.15
--Read_Queue_Size	500000
--Result_Queue_Size	100000

--N_Limit	0.1
--Minimum_Length	100	# Length after trimming

