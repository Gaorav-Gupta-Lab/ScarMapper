"""

@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27599
@copyright: 2019
"""
import datetime
import os
import collections
import subprocess
import argparse
import sys
import time
from distutils.util import strtobool
from scipy.stats import gmean
from scarmapper import TargetMapper as Target_Mapper, INDEL_Processing as Indel_Processing
from Valkyries import Tool_Box, Version_Dependencies as VersionDependencies, FASTQ_Tools


__author__ = 'Dennis A. Simpson'
__version__ = '0.1.0'


def main(command_line_args=None):
    VersionDependencies.python_check()

    if not command_line_args:
        command_line_args = sys.argv
    run_start = datetime.datetime.today().strftime("%a %b %d %H:%M:%S %Y")
    parser = argparse.ArgumentParser(description="A little ditty to manipulate FASTQ files.\n {0} v{1}"
                                     .format(__package__, __version__), formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--options_file', action='store', dest='options_file', required=True,
                        help='File containing program parameters.')

    options_parser = Tool_Box.options_file(parser)
    args = options_parser.parse_args()
    # args, options_parser = string_to_boolean(args, options_parser)
    options_parser.set_defaults(Trim5=0)
    options_parser.set_defaults(Trim3=0)
    options_parser.set_defaults(Minimum_Length=100)
    options_parser.set_defaults(N_Limit=100)
    options_parser.set_defaults(HaloPLEX=False)
    options_parser.set_defaults(ThruPLEX=False)
    options_parser.set_defaults(FASTQ_PreProcess=True)
    args = options_parser.parse_args()

    # Check options file for errors.
    error_checking(args)

    log = Tool_Box.Logger(args)
    Tool_Box.log_environment_info(log, args, command_line_args)
    start_time = time.time()
    module_name = ""

    # Initialize generator to read each FASTQ file
    fastq1 = FASTQ_Tools.FASTQ_Reader(args.FASTQ1, log)
    fastq2 = FASTQ_Tools.FASTQ_Reader(args.FASTQ2, log)
    index1 = FASTQ_Tools.FASTQ_Reader(args.Index1, log)
    index2 = FASTQ_Tools.FASTQ_Reader(args.Index2, log)

    splitter_data = FASTQ_Tools.FastqSplitter(args, log, fastq1, fastq2, index1, index2, paired_end=True)
    new_fastq1, new_fastq2 = splitter_data.file_writer()

    warning = "\033[1;31m **See warnings above**\033[m" if log.warning_occurred else ''
    elapsed_time = int(time.time() - start_time)
    log.info("****FASTQ Preprocessing {0} complete ({1} seconds, {2} Mb peak memory).****"
             .format(module_name, elapsed_time, Tool_Box.peak_memory(), warning))


def error_checking(args):
    """
    Check parameter file for errors.
    :param args:
    """
    if not os.path.exists(args.Working_Folder):
        print("\033[1;31mERROR:\n\tWorking Folder Path: {} Not Found.  Check Options File."
              .format(args.Working_Folder))
        raise SystemExit(1)

    if getattr(args, "FASTQ1", False) and not os.path.isfile(args.FASTQ1):
        print("\033[1;31mERROR:\n\t--FASTQ1: {} Not Found.  Check Options File."
              .format(args.FASTQ1))
        raise SystemExit(1)

    if getattr(args, "FASTQ2", False) and not os.path.isfile(args.FASTQ2):
        print("\033[1;31mERROR:\n\t--FASTQ2: {} Not Found.  Check Options File."
              .format(args.FASTQ2))
        raise SystemExit(1)

    if getattr(args, "Index1", False) and not os.path.isfile(args.Index1):
        print("\033[1;31mERROR:\n\t--Index1: {} Not Found.  Check Options File."
              .format(args.Index1))
        raise SystemExit(1)

    if getattr(args, "Index2", False) and not os.path.isfile(args.Index2):
        print("\033[1;31mERROR:\n\t--Index2: {} Not Found.  Check Options File."
              .format(args.Index1))
        raise SystemExit(1)


if __name__ == '__main__':
    main()