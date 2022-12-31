"""

@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27599
@copyright: 2022
"""

import collections
from Valkyries import Tool_Box

__author__ = 'Dennis A. Simpson'
__version__ = '1.0.0'
__package__ = 'ScarMapper'


class TargetMapper:
    def __init__(self, log, args, sample_manifest):
        """
        :param sample_manifest:
        :param log:
        :param args:
        """
        self.target_data = Tool_Box.FileParser.indices(log, args.TargetFile)
        # self.target_data = [x[0] for x in Tool_Box.FileParser.indices(log, args.TargetFile)]
        self.sample_manifest = sample_manifest
        self.log = log
        self.args = args

    @property
    def phasing(self):
        """
        This will build and return a dictionary that contains the primer phasing information.
        :return:
        """

        phasing_dict = collections.defaultdict(lambda: collections.defaultdict(list))

        # Current TruSeq does not use primer phasing.
        if self.args.Platform == "TruSeq":
            return phasing_dict

        for sample in self.sample_manifest:
            locus = sample[6]

            # Only need to process each locus 1x
            if locus in phasing_dict:
                continue

            '''
            While I believe it unlikely to ever be needed, this code allows the forward and reverse reads to have 
            different phasing lengths.
            '''
            forward_phase_seq = sample[4]
            reverse_phase_seq = sample[5]
            forward_seq_length = len(forward_phase_seq)
            reverse_seq_length = len(reverse_phase_seq)
            forward_phase_length = int(forward_seq_length*0.5)
            reverse_phase_length = int(reverse_seq_length*0.5)
            f_left_position = forward_seq_length - forward_phase_length
            r_left_position = reverse_seq_length - reverse_phase_length

            for i in range(forward_phase_length+1):
                f_phase = "F{}".format(i)
                r_phase = "R{}".format(i)
                forward_sequence = forward_phase_seq[f_left_position-i:forward_seq_length-i]
                reverse_sequence = reverse_phase_seq[r_left_position-i:reverse_seq_length-i]
                phasing_dict[locus]["R1"].append([forward_sequence, f_phase])
                phasing_dict[locus]["R2"].append([reverse_sequence, r_phase])

        return phasing_dict

    @property
    def targets(self):
        """
        This will build and return the dictionary that contains the sgRNAs
        :return:
        """
        target_dict = collections.defaultdict(tuple)
        for line in self.target_data:
            locus_name = line[0]
            chromosome = line[1]

            start_pos = int(line[2])
            stop_pos = int(line[3])

            sgrna_seq = line[4].upper()
            rcomp = line[5].upper()

            target_dict[locus_name] = \
                (locus_name, chromosome, start_pos, stop_pos, sgrna_seq, rcomp)

        return target_dict
