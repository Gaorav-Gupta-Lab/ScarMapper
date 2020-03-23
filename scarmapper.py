"""

@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27599
@copyright: 2020
"""
import datetime
import glob
import shutil
import os
import collections
import subprocess
import argparse
import sys
import time
from scipy import stats
from distutils.util import strtobool
from scipy.stats import gmean
from Valkyries import Tool_Box, Version_Dependencies as VersionDependencies, FASTQ_Tools

try:
    # If the cythonized file doesn't exist then create it.
    from scarmapper import INDEL_Processing as Indel_Processing, TargetMapper as Target_Mapper
except ImportError:
    setup_file = "python3 {0}{1}scarmapper{1}setup.py build_ext --inplace".format(os.path.dirname(__file__), os.sep)
    os.chdir(os.path.dirname(__file__))
    subprocess.run([setup_file], shell=True)
    from scarmapper import INDEL_Processing as Indel_Processing, TargetMapper as Target_Mapper

__author__ = 'Dennis A. Simpson'
__version__ = '0.13.0'
__package__ = 'ScarMapper'


def atropos_trim(args, log, method):
    """
    Trim adapters from reads with Atropos.  This creates the Atropos config file and runs Atropos.
    :param args:
    :param log:
    :param method:
    :return:
    """

    fastq1_trimmed = "{}{}_Trim.R1.fastq.gz".format(args.WorkingFolder, args.Job_Name)
    trim_report = "{}Atropos_{}_Trim_Report.txt".format(args.WorkingFolder, args.Job_Name)

    fastq2_trimmed = "{}{}_Trim.R2.fastq.gz".format(args.WorkingFolder, args.Job_Name)

    if args.NextSeq_Trim:
        nextseq_trim = "--nextseq-trim 1"
        op_order = "--op-order GAWCQ"
    else:
        nextseq_trim = ""
        op_order = "--op-order AWCQ"

    additional_adapters = ""
    # for adapter in user_sequences:
    #     additional_adapters += "-a {0}\n-A {0}\n".format(adapter)
    config_block = \
        "trim\n--aligner {0}\n--threads {1}\n{2}\n-G file:{3}\n-g file:{3}\n-A file:{4}\n-a file:{4}\n-o {5}\n-p {6}\n"\
        "-pe1 {7}\n-pe2 {8}\n{14}\n{9}\n--error-rate {10}\n--times 3\n--quality-cutoff 20\n"\
        "--stats bot\n--read-queue-size {12}\n--result-queue-size {13}\n--report-file {11}\n"\
        .format(args.Atropos_Aligner, args.Spawn, additional_adapters, args.Anchored_Adapters_5p,
                args.Anchored_Adapters_3p, fastq1_trimmed, fastq2_trimmed, args.FASTQ1, args.FASTQ2, nextseq_trim,
                args.Adapter_Mismatch_Fraction, trim_report, args.Read_Queue_Size, args.Result_Queue_Size, op_order)
    config_file = open("{}{}_Atropos_Config.txt".format(args.WorkingFolder, args.Job_Name), "w")
    config_file.write(config_block)
    config_file.close()

    log.info("Beginning Atropos Trim of {} library".format(method))
    subprocess.run("atropos --config {}{}_Atropos_Config.txt".format(args.WorkingFolder, args.Job_Name), shell=True)

    return fastq1_trimmed, fastq2_trimmed


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
    parser = argparse.ArgumentParser(description="A package to map genomic repair scars at defined loci.\n {0} v{1}"
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
        log.info("Sending FASTQ files to FASTQ preprocessor.")
        fastq_file1 = args.FASTQ1
        fastq_file2 = args.FASTQ2

        fq1 = FASTQ_Tools.FASTQ_Reader(fastq_file1, log)
        fq2 = FASTQ_Tools.FASTQ_Reader(fastq_file2, log)

        indel_processing = \
            Indel_Processing.DataProcessing(log, args, run_start, Target_Mapper.TargetMapper(log, args), fq1, fq2)

        indel_processing.main_loop()

    elif not args.IndelProcessing:
        # Run frequency file Combine module
        log.info("Process Replicates.")
        data_dict = collections.defaultdict(list)
        file_list = [f for f in glob.glob("{}*ScarMapper_Frequency.txt".format(args.DataFiles, ))]

        for file_name in file_list:
            # index_name = sample_info[0]
            # file_name = "{}{}_{}_ScarMapper_Frequency.txt".format(args.DataFiles, args.Job_Name, index_name)
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
            "# Frequency\tSEM\tScar Type\tLeft Deletions\tRight Deletions\tDeletion Size\tMicrohomology\t" \
            "Microhomology Size\tInsertion\tInsertion Size\tLeft Template\tRight Template\tConsensus Left Junction\t" \
            "Consensus Right Junction\tTarget Left Junction\tTarget Right Junction\tConsensus\tTarget Region\n"
        for key, row_list in data_dict.items():
            row_string = "\t".join(row_list[1])
            freq = gmean(row_list[0])
            sem = stats.sem(row_list[0])

            freq_results_outstring += "{}\t{}\t{}\n".format(freq, sem, row_string)

        freq_results_file = \
            open("{}{}_{}_ScarMapper_Combined_Frequency.txt".format(args.WorkingFolder, args.Job_Name, args.SampleName), "w")

        freq_results_file.write(freq_results_outstring)
        freq_results_file.close()

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
        # options_parser.set_defaults(Atropos_Trim=bool(strtobool(args.Atropos_Trim)))
        options_parser.set_defaults(Demultiplex=bool(strtobool(args.Demultiplex)))
        options_parser.set_defaults(OutputRawData=bool(strtobool(args.OutputRawData)))
        # options_parser.set_defaults(NextSeq_Trim=bool(strtobool(args.NextSeq_Trim)))
    options_parser.set_defaults(IndelProcessing=bool(strtobool(args.IndelProcessing)))
    options_parser.set_defaults(Verbose=args.Verbose.upper())

    return options_parser.parse_args()


if __name__ == '__main__':
    main()
