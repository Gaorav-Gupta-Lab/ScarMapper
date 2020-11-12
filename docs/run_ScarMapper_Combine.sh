#!/bin/bash
#Parameter file to run ScarMapper Combine module
#File generated 11-Nov-2020


python3 /full/path/to/scarmapper.py --options_file /full/path/to/run_ScarMapper_Combine.sh
exit

--IndelProcessing	False

--WorkingFolder	/full/path/to/file/save/location/
--SampleName	# Labels plot
--DataFiles	/full/path/to/frequency/files/

--Verbose	INFO
--Job_Name	# Labels output file

# Plot Options
--FigureType	pdf # svg, jpg, tiff, pdf, png