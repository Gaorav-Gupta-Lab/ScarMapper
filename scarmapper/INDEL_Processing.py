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
__version__ = '0.3.0'
__package__ = 'ScarMapper'


class ScarSearch:
    def __init__(self, log, args, target_dict):
        self.log = log
        self.args = args
        self.target_dict = target_dict
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
        Generate the consensus sequence and find indels
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
        junction_type_data = [0, 0, 0, 0]
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

        time_list = []

        for left_seq, right_seq in sequence_list:
            loop_count += 1

            if loop_count % 5000 == 0:
                self.log.info("Processed {} reads of {} for {} in {} seconds. Elapsed time: {} seconds."
                              .format(loop_count, len(sequence_list), index_name, time.time() - split_time,
                                      time.time() - start_time))
                split_time = time.time()

            # Generate a gapped alignment of the paired reads.
            read_alignment_dict = \
                self.gapped_aligner(">left\n{}\n>right\n{}\n"
                                    .format(left_seq, Sequence_Magic.rcomp(right_seq)))

            # Build a simple contig from the gapped alignment of the paired reads
            consensus_seq = ""
            first_lft = False
            for i, (lft, rt) in enumerate(zip(read_alignment_dict["left"], read_alignment_dict["right"])):
                if not first_lft:
                    if lft is not "-" and read_alignment_dict["left"][i + 1] is not "-":
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
            # start_time = time.time()

            sub_list, self.summary_data = \
                SlidingWindow.sliding_window(consensus_seq, self.target_region, self.cutsite, self.target_length,
                                             self.lower_limit, self.upper_limit, self.summary_data, self.left_target_windows, self.right_target_windows)
            # time_list.append(time.time() - start_time)

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

        # Tool_Box.debug_messenger(statistics.mean(time_list))

        # Get the number of unmodified reads for library and remove the key from dictionary.
        no_cut_list = results_freq_dict.pop("|||", [0])
        self.summary_data[6] = no_cut_list[0]

        self.log.info("Finished Processing {}".format(index_name))

        # Write frequency results file
        freq_results_file = \
            open("{}{}_{}_ScarMapper_Frequency.txt"
                 .format(self.args.Working_Folder, self.args.Job_Name, index_name), "w")

        freq_results_outstring = \
            "Total\tFrequency\tLeft Deletions\tRight Deletions\tDeletion Size\tMicrohomology\tMicrohomology Size\t" \
            "Insertion\tInsertion Size\tLeft Template\tRight Template\tConsensus\tTarget Region\n"

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
            lft_template = ""
            rt_template = ""

            # TMEJ counts
            if del_size > 3 and microhomology_size > 1:
                junction_type_data[0] += key_count
            # NHEJ counts
            elif del_size <= 3 and ins_size < 5 or del_size > 3 and microhomology_size <= 1:
                junction_type_data[1] += key_count
            # Large Insertions:
            elif ins_size >= 5:
                junction_type_data[2] += key_count
                lft_template, rt_template = \
                    self.templated_insertion_search(insertion, results_freq_dict[freq_key][1][5],
                                                    results_freq_dict[freq_key][1][6])
            # Scars not part of the previous three
            else:
                junction_type_data[3] += key_count

            freq_results_outstring += "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n" \
                .format(key_count, key_frequency, lft_del, rt_del, del_size, microhomology, microhomology_size,
                        insertion, ins_size, lft_template, rt_template, consensus, self.target_region)

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
        lft_query2 = Sequence_Magic.rcomp(insertion[-5:])
        rt_query1 = insertion[:5]
        rt_query2 = insertion[-5:]
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
        Generates and returns a gapped alignment from the given FASTA data using Muscle.
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

    def sliding_window(self, consensus):
        """
        The Cython version of this method is 17% faster.
        Find junctions using consensus read.  All numbering is relative to the 5' end of the reference sequence,
        in this case it is the target region.
        :param consensus:
        """
        self.log.warning("This method is depreciated.  Cython SlidingWindow is 17% faster")

        def window(seq, start, stop, n, fill=None, keep=0, step=1):
            # "Returns a sliding window generator (of width n) over data from the iterable"

            it = iter(seq)

            result = tuple(itertools.islice(it, start, n))
            # result_window = "".join(result)

            if len(result) == n:
                yield "".join(result)

            while True:
                #         for elem in it:
                elem = tuple(next(it, fill) for _ in range(step))
                result = result[step:] + elem
                if elem[-1] is fill:
                    if keep:
                        yield "".join(result)
                    break

                yield "".join(result)

        lower_limit = int(0.15 * len(self.target_region))
        upper_limit = len(self.target_region)*0.85
        upper_consensus_limit = 0.9*len(consensus)
        read_results_list = []
        ldel = ""
        rdel = ""

        left_found = False
        right_found = False

        target_lft_junction = self.cutsite
        target_rt_junction = self.cutsite
        consensus_lft_junction = 0
        consensus_rt_junction = 0

        # No need to analyze sequences that are too short.
        if len(consensus) <= 4*lower_limit:
            self.summary_data[7] += 1
            return read_results_list

        '''
        Find the 5' junction.  Start at the cut position, derived from the target region, and move toward the 5'
        end of the read one nucleotide at a time using a 10 nucleotide sliding window.  The 3' position of
        the first perfect match of the window from the query and target defines the 5' junction.
        '''
        rt_position = self.cutsite
        lft_position = rt_position-10
        consensus_rt_position = self.cutsite
        consensus_lft_position = consensus_rt_position-10
        # query_fragment = "".join(reversed(consensus[:self.cutsite]))
        # query_slider = window(query_fragment, start=0, stop=len(query_fragment) - lower_limit, n=10)
        while not left_found and consensus_lft_position > lower_limit:
            query_segment = consensus[consensus_lft_position:consensus_rt_position]
            # target_fragment = "".join(reversed(self.target_region[:self.cutsite]))
            # target_slider = window(target_fragment, start=0, stop=len(target_fragment)-lower_limit, n=10)
            # query_segment = next(query_slider)

            while not left_found and lft_position > lower_limit:
                # target_segment = next(target_slider)
                target_segment = self.target_region[lft_position:rt_position]

                if query_segment == target_segment:
                    left_found = True
                    target_lft_junction = rt_position
                    consensus_lft_junction = rt_position
                    ldel = self.target_region[target_lft_junction:self.cutsite]

                lft_position -= 1
                rt_position -= 1

            # reset the the target region window.
            lft_position = rt_position - 10
            rt_position = self.cutsite

            # Move the consensus window one nucleotide 5'
            consensus_lft_position -= 1
            consensus_rt_position -= 1

        '''
        Find the 3' junction.  Start at the cut position, derived from the target region, and move toward the 5'
        end of the read one nucleotide at a time using a 10 nucleotide sliding window.  One plus the 3' position of
        the first perfect match of the window from the query (FASTQ read 2) and target defines the 5' junction.  The
        query and targets have different numbering and the windows move in opposite directions.
        '''

        lft_position = self.cutsite
        rt_position = lft_position+10

        # Move to the expected cutsite position on the consensus from the 3' end.
        consensus_lft_position = len(consensus)-(len(self.target_region)-self.cutsite)
        consensus_rt_position = consensus_lft_position+10

        while not right_found and consensus_rt_position < upper_consensus_limit:
            query_segment = consensus[consensus_lft_position:consensus_rt_position]

            while not right_found and rt_position < upper_limit:
                target_segment = self.target_region[lft_position:rt_position]

                if query_segment == target_segment:
                    right_found = True
                    target_rt_junction = lft_position
                    consensus_rt_junction = consensus_lft_position
                    rdel = self.target_region[self.cutsite:target_rt_junction]

                lft_position += 1
                rt_position += 1

            # reset target window
            lft_position = self.cutsite
            rt_position = lft_position + 10

            # increment consensus window
            consensus_lft_position += 1
            consensus_rt_position += 1

        # extract the insertion
        consensus_insertion = ""
        if 0 < consensus_lft_junction < consensus_rt_junction:
            consensus_insertion = consensus[consensus_lft_junction:consensus_rt_junction]

            # If there is an N in the insertion then don't include read in the analysis.
            if "N" in consensus_insertion:
                self.summary_data[7] += 1
                return read_results_list

            self.summary_data[4] += 1

        # Count left deletions
        if target_lft_junction < self.cutsite:
            self.summary_data[2] += 1

        # Count right deletions
        if target_rt_junction > self.cutsite:
            self.summary_data[3] += 1

        # extract the microhomology
        consensus_microhomology = ""
        if consensus_lft_junction > consensus_rt_junction > 0:
            consensus_microhomology = consensus[consensus_rt_junction:consensus_lft_junction]
            if consensus_microhomology:
                self.summary_data[5] += 1

        # Count number of reads passing all filters
        self.summary_data[1] += 1

        read_results_list = [ldel, rdel, consensus_insertion, consensus_microhomology, consensus]

        return read_results_list


class DataProcessing:
    def __init__(self, log, args, fq1, fq2, run_start):
        self.log = log
        self.args = args
        self.date_format = "%a %b %d %H:%M:%S %Y"
        self.run_start = run_start
        self.fastq_outfile_dict = None
        self.fastq_data_dict = None
        self.index_dict, self.target_dict = self.dictionary_build()
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

    @staticmethod
    def microhomology_search(lft_junction, rt_junction, target_region, consensus):
        """
        Method no longer used.
        :param lft_junction:
        :param rt_junction:
        :param target_region:
        :param consensus:
        :return:
        """

        lower_limit = int(0.10 * len(target_region))
        upper_limit = int(0.90 * len(target_region))
        rt_microhomology = ""
        lft_microhomology = ""
        loop_stop = False
        left_step = lft_junction
        right_step = rt_junction

        # find the 3' microhomology
        while not loop_stop:
            if right_step >= upper_limit:
                loop_stop = True
            elif target_region[left_step] == target_region[right_step]:
                rt_microhomology += target_region[left_step]
            else:
                loop_stop = True

            right_step += 1
            left_step += 1

        # Find the 5' microhomology
        loop_stop = False
        left_step = lft_junction-1
        right_step = rt_junction-1
        while not loop_stop:
            if left_step <= lower_limit:
                loop_stop = True
            elif target_region[left_step] == target_region[right_step]:
                lft_microhomology += target_region[left_step]
            else:
                loop_stop = True

            right_step -= 1
            left_step -= 1

        return lft_microhomology, rt_microhomology

    def main_loop(self):
        self.log.info("Beginning main loop")

        self.demultiplex()

        self.log.info("Spawning {} Jobs to Process {} Libraries".format(self.args.Spawn, len(self.sequence_dict)))
        p = pathos.multiprocessing.Pool(int(self.args.Spawn))

        # My solution for passing key:value groups to through the multiprocessor.  Largest value ste goes first.
        data_list = []
        for key in sorted(self.sequence_dict, key=lambda k: len(self.sequence_dict[k]), reverse=True):
            data_list.append((key, self.sequence_dict[key]))

        # Not sure if clearing this is really necessary but it is not used again so why keep the RAM tied up.
        self.sequence_dict.clear()

        scar_search = ScarSearch(self.log, self.args, self.target_dict)

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

        index_list = Tool_Box.FileParser.indices(self.log, self.args.Index_File)
        index_dict = collections.defaultdict(list)

        for l in index_list:
            left_index_len = 6
            right_index_len = 6
            index_name = l[2]
            for left, right in zip(l[1], l[0]):
                if left.islower():
                    left_index_len += 1
                if right.islower():
                    right_index_len += 1

            index_dict[index_name] = \
                [Sequence_Magic.rcomp(l[1].upper()), left_index_len, l[0].upper(), right_index_len, l[3]]
            if self.args.Demultiplex:
                r1 = FASTQ_Tools.Writer(self.log, "{}{}_{}_R1.fastq"
                                        .format(self.args.Working_Folder, self.args.Job_Name, index_name))
                r2 = FASTQ_Tools.Writer(self.log, "{}{}_{}_R2.fastq"
                                        .format(self.args.Working_Folder, self.args.Job_Name, index_name))
                self.fastq_outfile_dict[index_name] = [r1, r2]

        targetmapper = Target_Mapper.TargetMapper(self.log, self.args)
        target_dict = targetmapper.targets

        return index_dict, target_dict

    @staticmethod
    def gapped_aligner(log, fasta_data, consensus=False):
        """
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

            # elif left_match <= 2:
            #     left_found = True
            # elif right_match <= 2:
            #     right_found = True

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
        #     if right_found:
        #         # No left index found.
        #         self.read_count_dict['unidentified'][0] += 1
        #     if left_found:
        #         # No right index found.
        #         self.read_count_dict['unidentified'][1] += 1
        #     if not left_found and not right_found:
        #         # Neither index found.
        #         self.read_count_dict['unidentified'][2] += 1

        return match_found, left_seq, right_seq, index_key, fastq1_read, fastq2_read

    def data_output(self, summary_data_list):
        """
        Format data and write summary file.
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
            "Index Name\tTotal Found\tFraction Total\tPassing Filters\tFraction Passing Filters\t" \
            "Number Filtered\t{}\tTMEJ\tNormalized TMEJ\tNHEJ\tNormalized NHEJ\tInsertion >=5\t" \
            "Normalized Insertion >=5\tOther Scar Type\n"\
            .format(sub_header)

        '''
        The data_list contains information for each library.  [0] index name; [1] reads passing all 
        filters; [2] reads with a left junction; [3] reads with a right junction; [4] reads with an insertion;
        [5] reads with microhomology; [6] reads with no identifiable cut; [7] filtered reads [8] scar type list.
        '''
        for data_list in summary_data_list:
            index_name = data_list[0]
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
            large_ins = data_list[8][2]
            large_ins_fraction = large_ins/cut
            other_scar = data_list[8][3]

            summary_outstring += \
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n"\
                .format(index_name, library_read_count,
                        fraction_all_reads, passing_filters, fraction_passing, filtered, cut, cut_fraction,
                        left_del, right_del, total_ins, microhomology, microhomology_fraction, tmej, tmej_fraction,
                        nhej, nhej_fraction, large_ins, large_ins_fraction, other_scar)

        summary_outstring += "Unidentified\t{}\t{}\n" \
            .format(self.read_count_dict["unidentified"], self.read_count_dict["unidentified"] / self.read_count)

        summary_file.write(summary_outstring)
        summary_file.close()
