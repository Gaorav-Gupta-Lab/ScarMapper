"""
@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27599
@copyright: 2021
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
from natsort import natsort
from scipy import stats
from distutils.util import strtobool
from scipy.stats import gmean
from Valkyries import Tool_Box, Version_Dependencies as VersionDependencies, FASTQ_Tools
from scarmapper import ScarMapperPlot
import re
import pathlib
from distutils.version import StrictVersion
import platform

# This is a seriously ugly hack to check the existence and age of the compiled file.
folder_content = os.listdir("{0}{1}scarmapper{1}".format(pathlib.Path(__file__).parent.absolute(), os.sep))
python_ver = StrictVersion(platform.python_version())

subver = 5
if "3.6.0" <= python_ver < "3.7.0":
    subver = 6
elif "3.7.0" <= python_ver < "3.8.0":
    subver = 7
elif "3.8.0" <= python_ver < "3.9.0":
    subver = 8

regex = re.compile("SlidingWindow.cpython-3{}.*.so".format(subver))
cfile = ""
old_file = False
for f in folder_content:
    if regex.search(f):
        cfile = f
        break
if cfile:
    cpath = "{0}{1}scarmapper{1}{2}".format(pathlib.Path(__file__).parent.absolute(), os.sep, cfile)
    pyx_file = "{0}{1}scarmapper{1}SlidingWindow.pyx".format(pathlib.Path(__file__).parent.absolute(), os.sep)
    compiled_time = pathlib.Path(cpath).stat().st_ctime
    pyx_module_time = pathlib.Path(pyx_file).stat().st_ctime

    if pyx_module_time >= compiled_time:
        old_file = True

if not cfile or old_file:
    print("Compiled Module Doesn't Exist or is Old; Compiling New SlidingWindow Module")
    
    setup_file = \
        "python3.{2} setup.py build_ext --inplace"\
        .format(pathlib.Path(__file__).parent.absolute(), os.sep, subver)

    os.chdir(pathlib.Path(__file__).parent.absolute())

    
    os.chdir("scarmapper")
    subprocess.run([setup_file], shell=True)
    os.chdir("..")
    # The sleep is to allow for network or disk latency.
    time.sleep(5.0)

from scarmapper import INDEL_Processing as Indel_Processing, TargetMapper as Target_Mapper

__author__ = 'Dennis A. Simpson'
__version__ = '0.26.1'
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
    test_method = ""
    if args.TestMethod:
        test_method = "-g {}".format(args.TestMethod)
    n = ""
    if args.MinConsensusLength:
        n = "-n {} ".format(args.MinConsensusLength)

    proc = subprocess.run(
        "{}{}Pear{}bin{}./pear -f {} -r {} -o {} {}{}{}{}{}{}{}"
        .format(pathlib.Path(__file__).parent.absolute(), os.sep, os.sep, os.sep, args.FASTQ1, args.FASTQ2,
                fastq_consensus_prefix, y, j, n, p_value, min_overlap, quality_threshold, phred_value, test_method),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

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

    run_start = datetime.datetime.today().strftime("%H:%M:%S %Y  %a %b %d")
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
        file_list = []
        if args.Platform == "Illumina" or args.Platform == "Ramsden" or args.Platform == "TruSeq":
            log.info("Sending FASTQ files to FASTQ preprocessor.")

            if args.PEAR:
                file_list = pear_consensus(args, log)
                if not file_list:
                    log.error("PEAR failed.  Check logs.")
                    raise SystemExit(1)
                fastq_consensus = file_list[0]

                fq1 = FASTQ_Tools.FASTQ_Reader(fastq_consensus, log)
                fq2 = None

            else:
                fq2 = FASTQ_Tools.FASTQ_Reader(args.FASTQ2, log)
                fq1 = FASTQ_Tools.FASTQ_Reader(args.FASTQ1, log)

            sample_manifest = Tool_Box.FileParser.indices(log, args.SampleManifest)
            indel_processing = \
                Indel_Processing.DataProcessing(log, args, run_start, __version__,
                                                Target_Mapper.TargetMapper(log, args, sample_manifest), fq1, fq2)

            indel_processing.main_loop()

            # Compress or delete PEAR files.
            if args.PEAR and file_list:
                if args.DeleteConsensusFASTQ:
                    log.info("Deleting PEAR FASTQ Files.")
                    Tool_Box.delete(file_list)
                else:
                    log.info("Compressing {} FASTQ Files Generated by PEAR.".format(len(file_list)))
                    p = pathos.multiprocessing.Pool(int(args.Spawn))
                    p.starmap(Tool_Box.compress_files, zip(file_list, itertools.repeat(log)))
        else:
            log.error("Only 'Illumina', 'TruSeq' or 'Ramsden' --Platform methods currently allowed.")
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

        plot_data_dict = collections.defaultdict(list)
        label_dict = collections.defaultdict(float)
        output_data_dict = collections.defaultdict(list)
        marker_list = []

        for key, row_list in data_dict.items():
            # Force pattern to be in at least half of the files.
            if len(row_list[0]) / file_count >= 0.5:
                row_string = "\t".join(row_list[1])
                freq = gmean(row_list[0])
                sem = stats.sem(row_list[0])
                freq_results_outstring = "{}\t{}\t{}\n".format(freq, sem, row_string)
                output_key = freq

                # Freq is a 17 digit float so it is very unlikely to be duplicated but if it is this increments it by
                # a small number then checks the uniqueness again.
                if output_key in output_data_dict:
                    output_key = output_key+1e-16
                    if output_key in output_data_dict:
                        output_key = output_key + 1e-16

                scar_type = row_list[1][0]
                label_dict[scar_type] += freq

                # Gather up our data for plotting
                lft_del = int(row_list[1][1])
                rt_del = int(row_list[1][2])
                mh_size = int(row_list[1][5])
                ins_size = int(row_list[1][7])

                output_data_dict[output_key] = \
                    [(freq, lft_del, rt_del, mh_size, ins_size, scar_type), freq_results_outstring]

        freq_results_outstring = \
            "{}# Frequency\tSEM\tScar Type\tLeft Deletions\tRight Deletions\tDeletion Size\tMicrohomology\t" \
            "Microhomology Size\tInsertion\tInsertion Size\tLeft Template\tRight Template\tConsensus Left Junction\t" \
            "Consensus Right Junction\tTarget Left Junction\tTarget Right Junction\tConsensus\tTarget Region\n" \
            .format(page_header)

        # Now draw a pretty graph of the data if we are not dealing with a negative control.
        for k in natsort.natsorted(output_data_dict, reverse=True):
            data_list = output_data_dict[k]
            freq_results_outstring += data_list[1]

            freq = data_list[0][0]
            lft_del = data_list[0][1]
            rt_del = data_list[0][2]
            mh_size = data_list[0][3]
            ins_size = data_list[0][4]
            scar_type = data_list[0][5]

            # Plotting all scar patterns is messy.  This provides a cutoff.
            if freq < 0.00025:
                continue

            y_value = freq * 0.5
            lft_ins_width = freq
            rt_ins_width = freq

            # This is gathered up to find the largest value.  Used to set the x-axis limits.
            marker_list.extend([lft_del + (mh_size * 0.5), rt_del + (mh_size * 0.5), ins_size])

            # Deletion size included half the size of any microhomology present.
            lft_del_plot_value = (lft_del + (mh_size * 0.5)) * -1
            rt_del_plot_value = rt_del + (mh_size * 0.5)

            # Insertions are centered on 0 so we need to take half the value for each side.
            lft_ins_plot_value = (ins_size * 0.5) * -1
            rt_ins_plot_value = ins_size * 0.5

            # Scale the width of bars for insertions inside of deletions
            if lft_del + (mh_size * 0.5) != 0:
                lft_ins_width = freq * 0.5
            if rt_del + (mh_size * 0.5) != 0:
                rt_ins_width = freq * 0.5

            if scar_type not in plot_data_dict:
                plot_data_dict[scar_type] = \
                    [[freq], [lft_del_plot_value], [rt_del_plot_value], [lft_ins_plot_value],
                     [rt_ins_plot_value], [lft_ins_width], [rt_ins_width], [y_value]]
            else:
                # Get some previous plot data
                count = len(plot_data_dict[scar_type][0])
                previous_freq = plot_data_dict[scar_type][0][count - 1]
                previous_y = plot_data_dict[scar_type][7][count - 1]

                plot_data_dict[scar_type][0].append(freq)
                plot_data_dict[scar_type][1].append(lft_del_plot_value)
                plot_data_dict[scar_type][2].append(rt_del_plot_value)
                plot_data_dict[scar_type][3].append(lft_ins_plot_value)
                plot_data_dict[scar_type][4].append(rt_ins_plot_value)
                plot_data_dict[scar_type][5].append(lft_ins_width)
                plot_data_dict[scar_type][6].append(rt_ins_width)

                # Use the previous plot data to find the y-value of the current bar.
                plot_data_dict[scar_type][7] \
                    .append(previous_y + 0.002 + (0.5 * previous_freq) + y_value)

        plot_data_dict['Marker'] = [(max(marker_list)) * -1, max(marker_list)]
        # sample_name = "{}.{}".format(args.Job_Name, args.SampleName)

        ScarMapperPlot.scarmapperplot(args, datafile=None, sample_name=args.SampleName, plot_data_dict=plot_data_dict,
                                      label_dict=label_dict)

        freq_results_file = \
            open("{}{}_ScarMapper_Combined_Frequency.txt".format(args.WorkingFolder, args.SampleName), "w")

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

    # if not os.path.exists(args.WorkingFolder):
    if not pathlib.Path(args.WorkingFolder).exists():
        print("\033[1;31mERROR:\n\tWorking Folder Path: {} Not Found.  Check Options File."
              .format(args.WorkingFolder))
        raise SystemExit(1)

    if getattr(args, "FASTQ1", False) and getattr(args, "ConsensusSequence", False):
        print("\033[1;31mERROR:\n\t--FASTQ1 and --ConsensusSequence both set.  Pick one or the other and try again.")
        raise SystemExit(1)

    if getattr(args, "FASTQ2", False) and getattr(args, "ConsensusSequence", False):
        print("\033[1;31mERROR:\n\t--FASTQ2 and --ConsensusSequence both set.  Pick one or the other and try again.")
        raise SystemExit(1)

    if getattr(args, "FASTQ1", False) and not pathlib.Path(args.FASTQ1).exists():
        print("\033[1;31mERROR:\n\t--FASTQ1: {} Not Found.  Check Options File."
              .format(args.FASTQ1))
        raise SystemExit(1)

    if getattr(args, "FASTQ2", False) and not pathlib.Path(args.FASTQ2).exists():
        print("\033[1;31mERROR:\n\t--FASTQ2: {} Not Found.  Check Options File."
              .format(args.FASTQ2))
        raise SystemExit(1)

    if getattr(args, "ConsensusSequence", False) and not pathlib.Path(args.ConsensusSequence).exists():
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
        # Tool_Box.debug_messenger("Pear set to FALSE.")
        options_parser.set_defaults(PEAR=True)
        options_parser.set_defaults(Demultiplex=bool(strtobool(args.Demultiplex)))
        options_parser.set_defaults(OutputRawData=bool(strtobool(args.OutputRawData)))
        options_parser.set_defaults(DeleteConsensusFASTQ=bool(strtobool(args.DeleteConsensusFASTQ)))

    options_parser.set_defaults(IndelProcessing=bool(strtobool(args.IndelProcessing)))
    options_parser.set_defaults(Verbose=args.Verbose.upper())

    return options_parser.parse_args()


if __name__ == '__main__':
    main()
