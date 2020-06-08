"""

@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27599
@copyright: 2020
"""

import collections
from Valkyries import Tool_Box, Sequence_Magic
import pysam

__author__ = 'Dennis A. Simpson'
__version__ = '0.10.0'
__package__ = 'ScarMapper'


class TargetMapper:
    def __init__(self, log, args):
        """

        :param log:
        :param args:
        """
        self.target_data = Tool_Box.FileParser.indices(log, args.TargetFile)
        self.phasing_data = Tool_Box.FileParser.indices(log, args.PrimerPhasingFile)
        # self.target_data = [x[0] for x in Tool_Box.FileParser.indices(log, args.TargetFile)]
        self.refseq = pysam.FastaFile(args.RefSeq)
        self.log = log
        self.args = args

    @property
    def phasing(self):
        """
        This will build and return a dictionary that contains the primer phasing information.
        :return:
        """

        phasing_dict = collections.defaultdict(lambda: collections.defaultdict(list))
        for row in self.phasing_data:
            phase = row[0]
            sequence = row[1]
            read = row[2]
            locus = row[3]
            # key = "{}+{}".format(locus, phase)

            phasing_dict[locus][read].append([sequence, phase])
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


