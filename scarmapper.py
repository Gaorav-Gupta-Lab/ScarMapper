"""

@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27599
@copyright: 2019
"""
import datetime
import os
import subprocess
import argparse
import sys
import time
from distutils.util import strtobool
import scarmapper.INDEL_Processing as Indel_Processing
from Valkyries import Tool_Box, Version_Dependencies as VersionDependencies, FASTQ_Tools


__author__ = 'Dennis A. Simpson'
__version__ = '0.3.0'
__package__ = 'ScarMapper'


def atropos_trim(args, log, fq1, fq2, method):
    """
    Trim adapters from reads with Atropos.  This creates the Atropos config file and runs Atropos.
    :param fq1:
    :param fq2:
    :param args:
    :param log:
    :param method:
    :return:
    """

    fastq1_trimmed = "{}{}_Trim.R1.fastq.gz".format(args.Working_Folder, args.Job_Name)
    trim_report = "{}Atropos_{}_Trim_Report.txt".format(args.Working_Folder, args.Job_Name)

    fastq2_trimmed = "{}{}_Trim.R2.fastq.gz".format(args.Working_Folder, args.Job_Name)

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
                args.Anchored_Adapters_3p, fastq1_trimmed, fastq2_trimmed, fq1, fq2, nextseq_trim,
                args.Adapter_Mismatch_Fraction, trim_report, args.Read_Queue_Size, args.Result_Queue_Size, op_order)
    config_file = open("{}{}_Atropos_Config.txt".format(args.Working_Folder, args.Job_Name), "w")
    config_file.write(config_block)
    config_file.close()

    log.info("Beginning Atropos Trim of {} library".format(method))
    subprocess.run("atropos --config {}{}_Atropos_Config.txt".format(args.Working_Folder, args.Job_Name), shell=True)

    return fastq1_trimmed, fastq2_trimmed


def main(command_line_args=None):
    """
    Let's get this party started.
    :param command_line_args:
    """
    VersionDependencies.python_check()

    if not command_line_args:
        command_line_args = sys.argv
    run_start = datetime.datetime.today().strftime("%a %b %d %H:%M:%S %Y")
    parser = argparse.ArgumentParser(description="A package to map genomic repair scars.\n {0} v{1}"
                                     .format(__package__, __version__), formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--options_file', action='store', dest='options_file', required=True,
                        help='File containing program parameters.')

    options_parser = Tool_Box.options_file(parser)
    args = options_parser.parse_args()
    error_checking(args)
    log = Tool_Box.Logger(args)
    Tool_Box.log_environment_info(log, args, command_line_args)
    start_time = time.time()
    module_name = ""
    args, options_parser = string_to_boolean(args, options_parser)
    log.info("Sending FASTQ files to FASTQ preprocessor.")
    # fastq_data = FASTQ_Tools.FastqSplitter(args, log, FASTQ_Tools.FASTQ_Reader(args.FASTQ1, log), FASTQ_Tools.FASTQ_Reader(args.FASTQ2, log))
    # fq1, fq2 = fastq_data.file_writer()
    fq1 = FASTQ_Tools.FASTQ_Reader(args.FASTQ1, log)
    fq2 = FASTQ_Tools.FASTQ_Reader(args.FASTQ2, log)

    if args.Atropos_Trim:
        fq1, fq2 = atropos_trim(args, log, fq1, fq2, "ScarMapper")

    indel_processing = Indel_Processing.DataProcessing(log, args, fq1, fq2, run_start)
    indel_processing.main_loop()

    warning = "\033[1;31m **See warnings above**\033[m" if log.warning_occurred else ''
    elapsed_time = int(time.time() - start_time)
    log.info("****ScarMapper {0} complete ({1} seconds, {2} Mb peak memory).****"
             .format(module_name, elapsed_time, Tool_Box.peak_memory(), warning))

    # All done so we need to quit otherwise Python will not release the log file on virtual Linux.
    exit(0)


def error_checking(args):
    if not os.path.exists(args.Working_Folder):
        print("\033[1;31mWARNING:\n\tWorking Folder Path: {} Not Found.  Check Options File."
              .format(args.Working_Folder))
        raise SystemExit(1)


def string_to_boolean(args, options_parser):
    """
    Converts strings to boolean.  Done to keep the eval function out of the code.
    :param args:
    :param options_parser:
    :return:
    """

    options_parser.set_defaults(ThruPLEX=bool(strtobool(args.ThruPLEX)))
    options_parser.set_defaults(Atropos_Trim=bool(strtobool(args.Atropos_Trim)))
    options_parser.set_defaults(Demultiplex=bool(strtobool(args.Demultiplex)))
    options_parser.set_defaults(ScarMapper=bool(strtobool(args.ScarMapper)))
    options_parser.set_defaults(OutputRawData=bool(strtobool(args.OutputRawData)))
    options_parser.set_defaults(NextSeq_Trim=bool(strtobool(args.NextSeq_Trim)))

    args = options_parser.parse_args()

    return args, options_parser


if __name__ == '__main__':
    main()
