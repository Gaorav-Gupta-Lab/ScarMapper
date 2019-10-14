#!/bin/bash
#Parameter file to run ScarMapper
#File generated 2019-07-19

# SLURM Commands
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=07-00:00:00
#SBATCH --mem=60g

python3 /full/path/to/scarmapper.py --options_file /full/path/to/run_ScarMapper_Combine.sh
exit


--IndelProcessing	False

--Index_File	/full/path/to/CombineIndex.bed
--Working_Folder	/full/path/to/working/folder/
--DataFiles	/full/path/to/data/files/

--Verbose	INFO
--Job_Name	<Must be the same as the one used for initial IndelProcessing run>
--SampleName	<From sample index>

