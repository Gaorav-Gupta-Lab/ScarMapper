"""

@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27599
@copyright: 2019
"""
import collections
import datetime
import itertools
import os
import statistics
import subprocess
import time
import pathos
import pysam
from Valkyries import Tool_Box, Sequence_Magic, FASTQ_Tools
from scarmapper import TargetMapper as Target_Mapper, SlidingWindow

__author__ = 'Dennis A. Simpson'
__version__ = '0.4.1'
__package__ = 'ScarMapper'


class ScarSearch:
    def __init__(self, log, args, target_dict, fastq):
        self.log = log
        self.args = args
        self.target_dict = target_dict
        self.fastq = fastq
        self.summary_data = None
        self.target_region = None
        self.sgrna = None
        self.cutsite = None
        self.lower_limit = None
        self.upper_limit = None
        self.target_length = None
        self.left_target_windows = []
        self.right_target_windows = []

    def window_mapping(self):
        """
        Predetermine all the sliding window results for the target region.
        """
        self.target_length = len(self.target_region)
        self.lower_limit = int(0.15 * self.target_length)
        self.upper_limit = int(0.85 * self.target_length)
        lft_position = self.cutsite-10
        rt_position = self.cutsite

        while rt_position > self.lower_limit:
            self.left_target_windows.append(self.target_region[lft_position:rt_position])
            lft_position -= 1
            rt_position -= 1

        lft_position = self.cutsite
        rt_position = self.cutsite+10
        while lft_position < self.upper_limit:
            self.right_target_windows.append(self.target_region[lft_position:rt_position])
            lft_position += 1
            rt_position += 1

    def data_processing(self, data_list):
        """
        Generate the consensus sequence and find indels.  Write the frequency file.  Called by pathos pool
        :param data_list:
        :return:
        """

        index_name = data_list[0]
        sequence_list = data_list[1]

        self.log.info("Begin Processing {}".format(index_name))
        """
        Summary List: index_name, total aberrant, left deletions, right deletions, total deletions, left insertions, 
        right insertions, total insertions, microhomology, number filtered
        """

        self.summary_data = [index_name, 0, 0, 0, 0, 0, 0, 0, 0]
        junction_type_data = [0, 0, 0, 0, 0]
        read_results_list = []
        results_freq_dict = collections.defaultdict(list)
        refseq = pysam.FastaFile(self.args.Ref_Seq)
        start = int(self.target_dict["Rosa26a"][1])
        stop = int(self.target_dict["Rosa26a"][2])
        chrm = self.target_dict["Rosa26a"][0]
        self.sgrna = self.target_dict["Rosa26a"][5]
        self.target_region = refseq.fetch(chrm, start, stop)
        self.cutsite_search(int(self.target_dict["Rosa26a"][3]) - int(self.target_dict["Rosa26a"][1]) - 50)
        self.window_mapping()
        loop_count = 0
        start_time = time.time()
        split_time = start_time

        for seq in sequence_list:

            loop_count += 1

            if loop_count % 5000 == 0:
                self.log.info("Processed {} reads of {} for {} in {} seconds. Elapsed time: {} seconds."
                              .format(loop_count, len(sequence_list), index_name, time.time() - split_time,
                                      time.time() - start_time))
                split_time = time.time()

            if self.fastq:
                # Generate a gapped simple consensus from the paired reads.
                left_seq, right_seq = seq

                consensus_seq = \
                    self.gapped_aligner(">left\n{}\n>right\n{}\n"
                                        .format(left_seq, Sequence_Magic.rcomp(right_seq)))
            else:
                consensus_seq = seq

            # No need to attempt an analysis of bad data.
            if consensus_seq.count("N") / len(consensus_seq) > float(self.args.N_Limit):
                self.summary_data[7] += 1
                continue

            # No need to analyze sequences that are too short.
            if len(consensus_seq) <= 4*self.lower_limit:
                self.summary_data[7] += 1
                continue

            '''
            The summary_data list contains information for a single library.  [0] index name; [1] reads passing all 
            filters; [2] reads with a left junction; [3] reads with a right junction; [4] reads with an insertion;
            [5] reads with microhomology; [6] reads with no identifiable cut; [7] filtered reads; [8] junction_type_data

            The junction_type_data list contains the repair type category counts.  [0] TMEJ, del > 3 and 
            mcirohomology > 1; [1] NHEJ, del <=3 and Ins < 5 or microhomology <= 1 and del > 3; [2] insertions >= 5 
            [3] Junctions with scars not represented by the other categories..
            '''

            sub_list, self.summary_data = \
                SlidingWindow.sliding_window(consensus_seq, self.target_region, self.cutsite, self.target_length,
                                             self.lower_limit, self.upper_limit, self.summary_data,
                                             self.left_target_windows, self.right_target_windows)

            '''
            The sub_list holds the data for a single consensus read.  These data are [left deletion, right deletion, 
            insertion, microhomology, consensus sequence].  The list could be empty if nothing was found or the 
            consensus was too short.
            '''

            if sub_list:
                read_results_list.append(sub_list)
                freq_key = "{}|{}|{}|{}".format(sub_list[0], sub_list[1], sub_list[2], sub_list[3])
            else:
                continue

            if freq_key in results_freq_dict:
                results_freq_dict[freq_key][0] += 1
            else:
                results_freq_dict[freq_key] = [1, sub_list]

        # Get the number of unmodified reads for library and remove the key from dictionary.
        no_cut_list = results_freq_dict.pop("|||", [0])
        self.summary_data[6] = no_cut_list[0]

        self.log.info("Finished Processing {}".format(index_name))

        # Write frequency results file
        freq_results_file = \
            open("{}{}_{}_ScarMapper_Frequency.txt"
                 .format(self.args.Working_Folder, self.args.Job_Name, index_name), "w")

        freq_results_outstring = \
            "# Total\tFrequency\tLeft Deletions\tRight Deletions\tDeletion Size\tMicrohomology\tMicrohomology Size\t" \
            "Insertion\tInsertion Size\tLeft Template\tRight Template\tConsensus Left Junction\t" \
            "Consensus Right Junction\tTarget Left Junction\tTarget Right Junction\tConsensus\tTarget Region\n"

        for freq_key in results_freq_dict:
            key_count = results_freq_dict[freq_key][0]
            key_frequency = key_count / (self.summary_data[1] - self.summary_data[6])
            lft_del = len(results_freq_dict[freq_key][1][0])
            rt_del = len(results_freq_dict[freq_key][1][1])
            insertion = results_freq_dict[freq_key][1][2]
            ins_size = len(insertion)
            consensus = results_freq_dict[freq_key][1][4]
            microhomology = results_freq_dict[freq_key][1][3]
            microhomology_size = len(microhomology)
            del_size = lft_del + rt_del + microhomology_size
            consensus_lft_junction = results_freq_dict[freq_key][1][5]
            consensus_rt_junction = results_freq_dict[freq_key][1][6]
            template_lft_junction = results_freq_dict[freq_key][1][7]
            template_rt_junction = results_freq_dict[freq_key][1][8]
            lft_template = ""
            rt_template = ""

            # TMEJ counts
            if del_size >= 4 and microhomology_size >= 2:
                junction_type_data[0] += key_count

            # NHEJ counts
            elif del_size < 4 and ins_size < 5:
                junction_type_data[1] += key_count

            # Non-Microhomology Deletions
            elif del_size >= 4 and microhomology_size < 2 and ins_size < 5:
                junction_type_data[4] += key_count

            # Large Insertions with or without Deletions:
            elif ins_size >= 5:
                junction_type_data[2] += key_count
                lft_template, rt_template = \
                    self.templated_insertion_search(insertion, template_lft_junction, template_rt_junction)

            # Scars not part of the previous four
            else:
                junction_type_data[3] += key_count

            freq_results_outstring += "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n" \
                .format(key_count, key_frequency, lft_del, rt_del, del_size, microhomology, microhomology_size,
                        insertion, ins_size, lft_template, rt_template, consensus_lft_junction, consensus_rt_junction,
                        template_lft_junction, template_rt_junction, consensus, self.target_region)

        freq_results_file.write(freq_results_outstring)
        freq_results_file.close()

        # add the junction list to the summary data
        self.summary_data[8] = junction_type_data

        # Format and output raw data if user has so chosen.
        if self.args.OutputRawData:
            self.raw_data_output(index_name, read_results_list)

        return self.summary_data

    def templated_insertion_search(self, insertion, lft_target_junction, rt_target_junction):
        """
        Search for left and right templates for insertions.
        :param insertion:
        :param lft_target_junction:
        :param rt_target_junction:
        :return:
        """
        lft_query1 = Sequence_Magic.rcomp(insertion[:5])
        lft_query2 = insertion[-5:]
        rt_query1 = insertion[:5]
        rt_query2 = Sequence_Magic.rcomp(insertion[-5:])
        lower_limit = lft_target_junction-50
        upper_limit = rt_target_junction+50
        left_not_found = True
        right_not_found = True
        lft_template = ""
        rt_template = ""

        # Set starting positions and search for left template
        lft_position = lft_target_junction-5
        rt_position = lft_target_junction

        while left_not_found and rt_position > lower_limit:
            target_segment = self.target_region[lft_position:rt_position]

            if lft_query1 == target_segment or lft_query2 == target_segment:
                lft_template = target_segment
                left_not_found = False

            lft_position -= 1
            rt_position -= 1

        # Reset starting positions and search for right template
        lft_position = rt_target_junction
        rt_position = rt_target_junction+5
        while right_not_found and lft_position < upper_limit:
            target_segment = self.target_region[lft_position:rt_position]
            if rt_query1 == target_segment or rt_query2 == target_segment:
                rt_template = target_segment
                right_not_found = False

            lft_position += 1
            rt_position += 1

        return lft_template, rt_template

    def raw_data_output(self, index_name, read_results_list):
        """
        Handle formatting and writing raw data.
        :param index_name:
        :param read_results_list:
        """
        results_file = open("{}{}_{}_ScarMapper.txt"
                            .format(self.args.Working_Folder, self.args.Job_Name, index_name), "w")

        results_outstring = "Left Deletions\tRight Deletions\tDeletion Size\tMicrohomology\tInsertion\t" \
                            "Insertion Size\tConsensus\tTarget Region\n"

        for data_list in read_results_list:
            lft_del = len(data_list[0])
            rt_del = len(data_list[1])
            microhomology = data_list[3]
            del_size = lft_del + rt_del + len(microhomology)
            total_ins = data_list[2]
            ins_size = len(total_ins)
            consensus = data_list[4]

            # skip unaltered reads.
            if del_size == 0 and ins_size == 0:
                continue

            results_outstring += "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n" \
                .format(lft_del, rt_del, del_size, microhomology, total_ins, ins_size, consensus, self.target_region)

        results_file.write(results_outstring)
        results_file.close()

    def cutsite_search(self, expected_cut):
        """
        Find the sgRNA cutsite on the gapped genomic DNA.
        :param expected_cut:
        """
        cutposition = 225
        target_site = self.sgrna[:6]  # found it best to search using only +/-3 nucleotides of the cutsite.
        cutsite_found = False
        new_position = expected_cut

        while not cutsite_found:
            if new_position > len(self.target_region) * 0.6:
                break

            guide_query = self.target_region[new_position:new_position + len(target_site)]

            if guide_query == target_site:
                cutposition = new_position + 3
                cutsite_found = True
            else:
                new_position += 1

        # This is a minor problem.  Insertions at cutsite cause errors.
        if not cutsite_found:
            while self.target_region[cutposition:cutposition + 1] == "-":
                cutposition += 1

        self.cutsite = cutposition

    def gapped_aligner(self, fasta_data):
        """
        Generates and returns a simple consensus from the given FASTA data using Muscle.
        :param: self
        :param: fasta_data
        :return:
        """

        # Create gapped alignment file in FASTA format using MUSCLE
        cmd = ['muscle', "-quiet", "-maxiters", "1", "-diags"]
        muscle = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)

        output, err = muscle.communicate(input=fasta_data)
        if err:
            self.log.error(err)

        cat_line = ""
        first_line = True
        gapped_alignment_dict = collections.defaultdict(str)
        key = ""
        list_output = list(output.splitlines())
        consensus_seq = ""

        for line in list_output:
            if ">" in line:
                if not first_line:
                    gapped_alignment_dict[key] = cat_line
                    cat_line = ""
                first_line = False
                key = line.split(">")[1].strip("\n")
            else:
                cat_line += line.strip("\n")

            gapped_alignment_dict[key] = cat_line

        # Build a simple contig from the gapped alignment of the paired reads
        first_lft = False
        for i, (lft, rt) in enumerate(zip(gapped_alignment_dict["left"], gapped_alignment_dict["right"])):

            if not first_lft:
                if lft is not "-" and gapped_alignment_dict["left"][i + 1] is not "-":
                    first_lft = True

            if not first_lft:
                continue

            if lft == rt:
                consensus_seq += lft
            elif lft is "-":
                consensus_seq += rt
            elif rt is "-":
                consensus_seq += lft
            else:
                consensus_seq += "N"

        return consensus_seq


class DataProcessing:
    def __init__(self, log, args, run_start, fq1=None, fq2=None, target_dict=None):
        self.log = log
        self.args = args
        self.date_format = "%a %b %d %H:%M:%S %Y"
        self.run_start = run_start
        self.fastq_outfile_dict = None
        self.fastq_data_dict = None
        self.target_dict = target_dict
        self.index_dict = self.dictionary_build()
        self.refseq = pysam.FastaFile(args.Ref_Seq)
        self.results_dict = collections.defaultdict(list)
        self.sequence_dict = collections.defaultdict(list)
        self.read_count_dict = collections.defaultdict()
        self.fastq1 = fq1
        self.fastq2 = fq2
        self.read_count = 0

    def demultiplex(self):
        """
        Finds reads by index.  Handles writing demultiplexed FASTQ if user desired.
        """
        self.log.info("Index Search")
        # fastq1_short_count = 0
        # fastq2_short_count = 0
        eof = False
        start_time = time.time()
        split_time = time.time()
        fastq_file_name_list = []

        while not eof:
            # Debugging Code Block
            if self.args.Verbose == "DEBUG":
                read_limit = 100000
                if self.read_count > read_limit:
                    if self.args.Demultiplex:
                        for index_name in self.fastq_data_dict:
                            r1_data = self.fastq_data_dict[index_name]["R1"]
                            r2_data = self.fastq_data_dict[index_name]["R2"]
                            r1, r2 = self.fastq_outfile_dict[index_name]
                            r1.write(r1_data)
                            r2.write(r2_data)
                            r1.close()
                            r2.close()
                    Tool_Box.debug_messenger("Limiting Reads Here to {}".format(read_limit))
                    eof = True

            try:
                fastq1_read = next(self.fastq1.seq_read())
                fastq2_read = next(self.fastq2.seq_read())

            except StopIteration:
                if self.args.Demultiplex:
                    for index_name in self.fastq_data_dict:
                        r1_data = self.fastq_data_dict[index_name]["R1"]
                        r2_data = self.fastq_data_dict[index_name]["R2"]
                        r1, r2 = self.fastq_outfile_dict[index_name]
                        r1.write(r1_data)
                        r2.write(r2_data)
                        r1.close()
                        r2.close()
                eof = True
                continue

            # ToDo: Short read scores should be incorporated into the summary file.
            if len(fastq1_read.seq) < int(self.args.Minimum_Length):
                # fastq1_short_count += 1
                continue
            elif len(fastq2_read.seq) < int(self.args.Minimum_Length):
                # fastq2_short_count += 1
                continue

            self.read_count += 1
            if self.read_count % 100000 == 0:
                elapsed_time = int(time.time() - start_time)
                block_time = int(time.time() - split_time)
                split_time = time.time()
                self.log.info("Processed {} reads in {} seconds.  Total elapsed time: {} seconds."
                              .format(self.read_count, block_time, elapsed_time))

            # Match read with library index.
            match_found, left_seq, right_seq, index_name, fastq1_read, fastq2_read = \
                self.index_matching(fastq1_read, fastq2_read)

            if match_found:
                self.sequence_dict[index_name].append([left_seq, right_seq])
                if self.args.Demultiplex:
                    self.fastq_data_dict[index_name]["R1"].append(fastq1_read)
                    self.fastq_data_dict[index_name]["R2"].append(fastq2_read)
                    fastq_file_name_list.append("{}{}_{}_R1.fastq"
                                                .format(self.args.Working_Folder, self.args.Job_Name, index_name))
                    fastq_file_name_list.append("{}{}_{}_R2.fastq"
                                                .format(self.args.Working_Folder, self.args.Job_Name, index_name))
        if self.args.Demultiplex:
            fastq_file_name_list = list(set(fastq_file_name_list))
            self.log.info("Spawning {} Jobs to Compress {} Files.".format(self.args.Spawn, len(fastq_file_name_list)))

            p = pathos.multiprocessing.Pool(int(self.args.Spawn))
            p.starmap(Tool_Box.compress_files, zip(fastq_file_name_list, itertools.repeat(self.log)))

            self.log.info("All Files Compressed")

    def main_loop(self):
        """
        Main entry point for repair scar search and processing.
        """

        self.log.info("Beginning main loop|Demultiplexing FASTQ")

        self.demultiplex()

        self.log.info("Spawning {} Jobs to Process {} Libraries".format(self.args.Spawn, len(self.sequence_dict)))
        p = pathos.multiprocessing.Pool(int(self.args.Spawn))

        # My solution for passing key:value groups to through the multiprocessor.  Largest value group goes first.
        data_list = []
        for key in sorted(self.sequence_dict, key=lambda k: len(self.sequence_dict[k]), reverse=True):
            data_list.append((key, self.sequence_dict[key]))

        # Not sure if clearing this is really necessary but it is not used again so why keep the RAM tied up.
        self.sequence_dict.clear()

        # Initialize the class for doing the search.  Copies of this initialized class are passed to each job.
        scar_search = ScarSearch(self.log, self.args, self.target_dict, fastq=True)

        self.data_output(p.starmap(scar_search.data_processing, zip(data_list)))

        self.log.info("Main Loop Finished")

    def dictionary_build(self):
        """
        Build the index dictionary from the index list.
        :return:
        """

        self.log.info("Building DataFrames.")

        # If we are demultiplexing the input FASTQ then setup the output files and dataframe.
        if self.args.Demultiplex:
            self.fastq_outfile_dict = collections.defaultdict(list)
            self.fastq_data_dict = collections.defaultdict(lambda: collections.defaultdict(list))

        master_index_dict = {}
        with open(self.args.Master_Index_File) as f:
            for l in f:
                if "#" in l or not l:
                    continue
                l_list = [x for x in l.strip("\n").split("\t")]
                master_index_dict[l_list[0]] = [l_list[1], l_list[2]]

        sample_index_list = Tool_Box.FileParser.indices(self.log, self.args.Index_File)
        index_dict = collections.defaultdict(list)

        for sample in sample_index_list:
            index_name = sample[0]
            if index_name in index_dict:
                self.log.error("The index {0} is duplicated.  Correct the error in {1} and try again."
                               .format(sample[0], self.args.Index_File))
                raise SystemExit(1)
            sample_name = sample[1]
            sample_replicate = sample[2]
            if self.args.Platform == "Illumina":
                self.log.error("The Illumina --Platform method is not yet implemented.")
                raise SystemExit(1)

            elif self.args.Platform == "Ramsden":
                # This is for the Ramsden indexing primers.
                left_index_len = 6
                right_index_len = 6
                left_index_sequence, right_index_sequence = master_index_dict[index_name]

                for left, right in zip(left_index_sequence, right_index_sequence):
                    if left.islower():
                        left_index_len += 1
                    if right.islower():
                        right_index_len += 1

                index_dict[index_name] = \
                    [Sequence_Magic.rcomp(right_index_sequence.upper()), left_index_len, left_index_sequence.upper(),
                     right_index_len, index_name, sample_name, sample_replicate]
            else:
                self.log.error("Only 'Illumina' or 'Ramsden' --Platform methods allowed.")
                raise SystemExit(1)

            if self.args.Demultiplex:
                r1 = FASTQ_Tools.Writer(self.log, "{}{}_{}_R1.fastq"
                                        .format(self.args.Working_Folder, self.args.Job_Name, index_name))
                r2 = FASTQ_Tools.Writer(self.log, "{}{}_{}_R2.fastq"
                                        .format(self.args.Working_Folder, self.args.Job_Name, index_name))
                self.fastq_outfile_dict[index_name] = [r1, r2]

        return index_dict

    @staticmethod
    def gapped_aligner(log, fasta_data, consensus=False):
        """
        This approach is depreciated.
        Generates and returns a gapped alignment from the given FASTA data using Muscle.
        :param: log
        :param: fasta_data
        :return:
        """
        # Create gapped alignment file in FASTA format using MUSCLE , "-gapopen", "-12"
        # ToDo: would like to control the gap penalties to better find insertions.
        cmd = ['muscle', "-quiet", "-maxiters", "1", "-diags"]
        if consensus:
            cmd = ['muscle', "-quiet", "-maxiters", "1", "-diags"]

        muscle = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)

        output, err = muscle.communicate(input=fasta_data)
        if err:
            log.error(err)
        cat_line = ""
        first_line = True
        gapped_alignment_dict = collections.defaultdict(str)
        key = ""
        list_output = list(output.splitlines())

        for line in list_output:
            if ">" in line:
                if not first_line:
                    gapped_alignment_dict[key] = cat_line
                    cat_line = ""
                first_line = False
                key = line.split(">")[1].strip("\n")
            else:
                cat_line += line.strip("\n")

            gapped_alignment_dict[key] = cat_line

        return gapped_alignment_dict

    def index_matching(self, fastq1_read, fastq2_read):
        """
        This matches an index sequence with the index found in the sequence reads.
        :param fastq1_read:
        :param fastq2_read:
        :return:
        """

        match_found = False
        # left_found = False
        # right_found = False
        left_seq = ""
        right_seq = ""
        index_key = 'unidentified'

        for index_key in self.index_dict:
            left_index = self.index_dict[index_key][0]
            right_index = self.index_dict[index_key][2]
            left_trim = self.index_dict[index_key][1]
            right_trim = self.index_dict[index_key][3]

            left_match = Sequence_Magic.match_maker(left_index, fastq2_read.seq[:len(left_index)])
            right_match = Sequence_Magic.match_maker(right_index, fastq1_read.seq[:len(right_index)])

            if index_key not in self.read_count_dict:
                self.read_count_dict[index_key] = 0

            if left_match <= 2 and right_match <= 2:
                self.read_count_dict[index_key] += 1
                match_found = True

            if match_found:
                # Trim the staggered nucleotides from the ends of the reads and pull the left and right sequence.
                FASTQ_Tools.read_trim(fastq2_read, trim5=left_trim)
                FASTQ_Tools.read_trim(fastq1_read, trim5=right_trim)
                left_seq = fastq2_read.seq
                right_seq = fastq1_read.seq
                break

        if not match_found:
            if 'unidentified' not in self.read_count_dict:
                self.read_count_dict['unidentified'] = 0
            self.read_count_dict['unidentified'] += 1

        return match_found, left_seq, right_seq, index_key, fastq1_read, fastq2_read

    def data_output(self, summary_data_list):
        """
        Format data and write the summary file.
        :param summary_data_list:
        """

        self.log.info("Formatting data and writing summary file")

        summary_file = open("{}{}_ScarMapper_Summary.txt".format(self.args.Working_Folder, self.args.Job_Name), "w")
        sub_header = \
            "Cut\tCut Fraction\tLeft Deletion Count\tRight Deletion Count\tInsertion Count\tMicrohomology Count\t" \
            "Normalized Microhomology"
        run_stop = datetime.datetime.today().strftime(self.date_format)
        summary_outstring = "ScarMapper {}\nStart: {}\nEnd: {}\nFASTQ1: {}\nFASTQ2: {}\nReads Analyzed: {}\n\n"\
            .format(__version__, self.run_start, run_stop, self.args.FASTQ1, self.args.FASTQ2, self.read_count, )
        summary_outstring += \
            "Index Name\tSample Name\tSample Replicate\tTotal Found\tFraction Total\tPassing Filters\t" \
            "Fraction Passing Filters\tNumber Filtered\t{}\tTMEJ\tNormalized TMEJ\tNHEJ\tNormalized NHEJ\t" \
            "Non-Microhomology Deletions\tNormalized Non-MH Del\tInsertion >=5 +/- Deletions\t" \
            "Normalized Insertion >=5+/- Deletions\tOther Scar Type\n"\
            .format(sub_header)

        '''
        The data_list contains information for each library.  [0] index name; [1] reads passing all 
        filters; [2] reads with a left junction; [3] reads with a right junction; [4] reads with an insertion;
        [5] reads with microhomology; [6] reads with no identifiable cut; [7] filtered reads [8] scar type list.
        '''
        for data_list in summary_data_list:
            index_name = data_list[0]
            sample_name = self.index_dict[index_name][5]
            sample_replicate = self.index_dict[index_name][6]
            library_read_count = self.read_count_dict[index_name]
            fraction_all_reads = library_read_count/self.read_count
            passing_filters = data_list[1]
            fraction_passing = passing_filters/library_read_count
            left_del = data_list[2]
            right_del = data_list[3]
            total_ins = data_list[4]
            microhomology = data_list[5]
            cut = passing_filters-data_list[6]
            cut_fraction = cut/passing_filters
            microhomology_fraction = microhomology/cut
            filtered = data_list[7]
            tmej = data_list[8][0]
            tmej_fraction = tmej/cut
            nhej = data_list[8][1]
            nhej_fraction = nhej/cut
            non_microhomology_del = data_list[8][4]
            non_mh_del_fraction = non_microhomology_del/cut
            large_ins = data_list[8][2]
            large_ins_fraction = large_ins/cut
            other_scar = data_list[8][3]

            summary_outstring += \
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n"\
                .format(index_name, sample_name, sample_replicate, library_read_count,
                        fraction_all_reads, passing_filters, fraction_passing, filtered, cut, cut_fraction,
                        left_del, right_del, total_ins, microhomology, microhomology_fraction, tmej, tmej_fraction,
                        nhej, nhej_fraction, non_microhomology_del, non_mh_del_fraction, large_ins, large_ins_fraction,
                        other_scar)

        summary_outstring += "Unidentified\t{}\t{}\n" \
            .format(self.read_count_dict["unidentified"], self.read_count_dict["unidentified"] / self.read_count)

        summary_file.write(summary_outstring)
        summary_file.close()
