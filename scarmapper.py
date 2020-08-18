"""

@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27599
@copyright: 2020
"""
import csv
import datetime
import glob
import itertools
import os
import collections
import subprocess
import argparse
import sys
import time
import pathos
from scipy import stats
from distutils.util import strtobool
from scipy.stats import gmean
from Valkyries import Tool_Box, Version_Dependencies as VersionDependencies, FASTQ_Tools
import re

# This is a seriously ugly hack to check the existence and age of the compiled file.
folder_content = os.listdir("{0}{1}scarmapper{1}".format(os.path.dirname(__file__), os.sep))
regex_pattern = "SlidingWindow.cpython.*\.so"
regex = re.compile(regex_pattern)
cfile = ""
old_file = False
for f in folder_content:
    if regex.search(f):
        cfile = f
        break
if cfile:
    compiled_time = \
        time.ctime(os.path.getmtime("{0}{1}scarmapper{1}{2}".format(os.path.dirname(__file__), os.sep, cfile)))
    pyx_module_time = \
        time.ctime(os.path.getmtime("{0}{1}scarmapper{1}SlidingWindow.pyx".format(os.path.dirname(__file__), os.sep)))
    if pyx_module_time > compiled_time:
        old_file = True

if not cfile or old_file:
    print("Compiled Module Doesn't Exist or is Old; Compiling New SlidingWindow Module")
    setup_file = "python3 {0}{1}scarmapper{1}setup.py build_ext --inplace".format(os.path.dirname(__file__), os.sep)
    os.chdir(os.path.dirname(__file__))
    subprocess.run([setup_file], shell=True)
    # The sleep is to allow for network or disk latency.
    time.sleep(5.0)

from scarmapper import INDEL_Processing as Indel_Processing, TargetMapper as Target_Mapper

__author__ = 'Dennis A. Simpson'
__version__ = '0.19.0'
__package__ = 'ScarMapper'


def pear_consensus(args, log):
    """
    This will take the input FASTQ files and use PEAR to generate a consensus file.
    :param args:
    :param log:
    :return:
    """
    log.info("Beginning PEAR Consensus")

    fastq_consensus_prefix = "{}{}".format(args.WorkingFolder, args.Job_Name)
    fastq_consensus_file = "{}.assembled.fastq".format(fastq_consensus_prefix)
    discarded_fastq = "{}.discarded.fastq".format(fastq_consensus_prefix)
    r1_unassembled = "{}.unassembled.forward.fastq".format(fastq_consensus_prefix)
    r2_unassembled = "{}.unassembled.reverse.fastq".format(fastq_consensus_prefix)

    y = "-y {} ".format(args.Memory)
    j = "-j {} ".format(int(args.Spawn)-1)
    p_value = ''
    if args.PValue:
        p_value = "-p {} ".format(args.PValue)
    min_overlap = ''
    if args.MinOverlap:
        min_overlap = "-v {} ".format(args.MinOverlap)
    quality_threshold = ""
    if args.QualityThreshold:
        quality_threshold = "-q {} ".format(args.QualityThreshold)
    phred_value = ""
    if args.PhredValue:
        phred_value = "-b {} ".format(args.PhredValue)

    proc = subprocess.run(
        "{}{}Pear{}bin{}./pear -f {} -r {} -o {} {}{}{}{}{}{}"
        .format(os.path.dirname(__file__), os.sep, os.sep, os.sep, args.FASTQ1, args.FASTQ2, fastq_consensus_prefix, y, j, p_value, min_overlap,
                quality_threshold, phred_value), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

    if proc.stderr:
        log.error("{}\n{}\n".format(proc.stderr.decode(), proc.stdout.decode()))
        return
    else:
        log.info(
        "Begin PEAR Output\n"
        "----------------------------------------------------------------------------------------------------------\n{}"
        "\n----------------------------------------------------------------------------------------------------------\n"
        .format(proc.stdout.decode()))

    file_list = [fastq_consensus_file, r1_unassembled, r2_unassembled]
    if os.stat(discarded_fastq).st_size > 0:
        file_list.append(discarded_fastq)
    else:
        Tool_Box.delete([discarded_fastq])

    return file_list


def main(command_line_args=None):
    """
    Let's get this party started.
    :param command_line_args:
    """
    start_time = time.time()
    VersionDependencies.python_check()

    if not command_line_args:
        command_line_args = sys.argv

    run_start = datetime.datetime.today().strftime("%a %b %d %H:%M:%S %Y")
    parser = argparse.ArgumentParser(description="A package to map genomic repair scars at defined loci.\n {} v{}"
                                     .format(__package__, __version__), formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--options_file', action='store', dest='options_file', required=True,
                        help='File containing program parameters.')

    # Check options file for errors and return object.
    args = error_checking(string_to_boolean(parser))

    log = Tool_Box.Logger(args)
    Tool_Box.log_environment_info(log, args, command_line_args)

    module_name = ""
    log.info("{} v{}".format(__package__, __version__))

    if args.IndelProcessing:
        if args.Platform == "Illumina" or args.Platform == "Ramsden":
            log.info("Sending FASTQ files to FASTQ preprocessor.")

            if args.PEAR:
                file_list = pear_consensus(args, log)
                if not file_list:
                    log.error("PEAR failed.  Check logs.")
                    return
                fastq_consensus = file_list[0]
                fq1 = FASTQ_Tools.FASTQ_Reader(fastq_consensus, log)
                fq2 = None
            else:
                fq2 = FASTQ_Tools.FASTQ_Reader(args.FASTQ2, log)
                fq1 = FASTQ_Tools.FASTQ_Reader(args.FASTQ1, log)

            indel_processing = \
                Indel_Processing.DataProcessing(log, args, run_start, __version__,
                                                Target_Mapper.TargetMapper(log, args), fq1, fq2)

            indel_processing.main_loop()
        else:
            log.error("Only 'Illumina' or 'Ramsden' --Platform methods currently allowed.")
            raise SystemExit(1)

    elif not args.IndelProcessing:
        # Run frequency file Combine module
        run_start = datetime.datetime.today().strftime("%a %b %d %H:%M:%S %Y")
        log.info("Process Replicates.")
        data_dict = collections.defaultdict(list)
        file_list = [f for f in glob.glob("{}*ScarMapper_Frequency.txt".format(args.DataFiles, ))]
        file_count = len(file_list)
        page_header = "# ScarMapper File Merge v{}\n# Run: {}\n# Sample Name: {}\n" \
            .format(__version__, run_start, args.SampleName)

        line_num = 0
        index_file = list(csv.reader(open(file_list[0]), delimiter='\t'))
        for line in index_file:
            if not line:
                break
            elif line_num > 3:
                page_header += "{}\n".format(line[0])

            line_num += 1
        page_header += "\n\n"

        for file_name in file_list:
            freq_file_data = Tool_Box.FileParser.indices(log, file_name)

            for row in freq_file_data:
                key = "{}|{}|{}|{}".format(row[3], row[4], row[6], row[8])
                row_data = row[2:]

                if key in data_dict:
                    data_dict[key][0].append(float(row[1]))
                else:
                    data_dict[key] = [[float(row[1])], row_data]

        # Process Data and Write Combined Frequency results file
        freq_results_outstring = \
            "{}# Frequency\tSEM\tScar Type\tLeft Deletions\tRight Deletions\tDeletion Size\tMicrohomology\t" \
            "Microhomology Size\tInsertion\tInsertion Size\tLeft Template\tRight Template\tConsensus Left Junction\t" \
            "Consensus Right Junction\tTarget Left Junction\tTarget Right Junction\tConsensus\tTarget Region\n" \
                .format(page_header)

        for key, row_list in data_dict.items():
            if len(row_list[0]) / file_count >= 0.5:
                row_string = "\t".join(row_list[1])
                freq = gmean(row_list[0])
                sem = stats.sem(row_list[0])

                freq_results_outstring += "{}\t{}\t{}\n".format(freq, sem, row_string)

        freq_results_file = \
            open("{}{}_{}_ScarMapper_Combined_Frequency.txt"
                 .format(args.WorkingFolder, args.Job_Name, args.SampleName), "w")

        freq_results_file.write(freq_results_outstring)
        freq_results_file.close()

    # Compress PEAR files.
    if args.PEAR:
        log.info("Compressing {} FASTQ Files Generated by PEAR.".format(len(file_list)))
        p = pathos.multiprocessing.Pool(int(args.Spawn))
        p.starmap(Tool_Box.compress_files, zip(file_list, itertools.repeat(log)))

    warning = "\033[1;31m **See warnings above**\033[m" if log.warning_occurred else ''
    elapsed_time = int(time.time() - start_time)
    log.info("****ScarMapper {0} complete ({1} seconds, {2} Mb peak memory).****"
             .format(module_name, elapsed_time, Tool_Box.peak_memory(), warning))

    # All done so we need to quit otherwise Python will not release the log file on virtual Linux.
    exit(0)


def error_checking(args):
    """
    Check parameter file for errors.
    :return:
    :param args:
    """

    if not os.path.exists(args.WorkingFolder):
        print("\033[1;31mERROR:\n\tWorking Folder Path: {} Not Found.  Check Options File."
              .format(args.WorkingFolder))
        raise SystemExit(1)

    if getattr(args, "FASTQ1", False) and getattr(args, "ConsensusSequence", False):
        print("\033[1;31mERROR:\n\t--FASTQ1 and --ConsensusSequence both set.  Pick one or the other and try again.")
        raise SystemExit(1)

    if getattr(args, "FASTQ2", False) and getattr(args, "ConsensusSequence", False):
        print("\033[1;31mERROR:\n\t--FASTQ2 and --ConsensusSequence both set.  Pick one or the other and try again.")
        raise SystemExit(1)

    if getattr(args, "FASTQ1", False) and not os.path.isfile(args.FASTQ1):
        print("\033[1;31mERROR:\n\t--FASTQ1: {} Not Found.  Check Options File."
              .format(args.FASTQ1))
        raise SystemExit(1)

    if getattr(args, "FASTQ2", False) and not os.path.isfile(args.FASTQ2):
        print("\033[1;31mERROR:\n\t--FASTQ2: {} Not Found.  Check Options File."
              .format(args.FASTQ2))
        raise SystemExit(1)

    if getattr(args, "ConsensusSequence", False) and not os.path.isfile(args.ConsensusSequence):
        print("\033[1;31mERROR:\n\t--ConsensusSequence: {} Not Found.  Check Options File."
              .format(args.FASTQ2))
        raise SystemExit(1)

    return args


def string_to_boolean(parser):
    """
    Converts strings to boolean.  Done to keep the eval() function out of the code.
    :param parser:
    :return:
    """
    options_parser = Tool_Box.options_file(parser)
    args = options_parser.parse_args()

    if args.IndelProcessing == "True":
        # Tool_Box.debug_messenger("Pear set to false.")
        options_parser.set_defaults(PEAR=True)
        options_parser.set_defaults(Demultiplex=bool(strtobool(args.Demultiplex)))
        options_parser.set_defaults(OutputRawData=bool(strtobool(args.OutputRawData)))

    options_parser.set_defaults(IndelProcessing=bool(strtobool(args.IndelProcessing)))
    options_parser.set_defaults(Verbose=args.Verbose.upper())

    return options_parser.parse_args()


if __name__ == '__main__':
    main()
