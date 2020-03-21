#!/bin/bash
#Parameter file to run ScarMapper Combine module
#File generated 2020-03-21


python3 /full/path/to/scarmapper.py --options_file /full/path/to/run_ScarMapper_Combine.sh
exit


--IndelProcessing	False

--SampleManifest	/full/path/to/SampleManifest_File.csv
--WorkingFolder	/full/path/to/working/folder/
--DataFiles	/full/path/to/data/files/

--Verbose	INFO
--Job_Name	# Use same Job Name as original Indel Processing run.
--SampleName	# This will be part of the output file name

