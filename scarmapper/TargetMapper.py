"""

@author: Dennis A. Simpson
         University of North Carolina at Chapel Hill
         Chapel Hill, NC  27599
@copyright: 2020
"""

import collections
import Valkyries.Tool_Box as Tool_Box
import pysam

__author__ = 'Dennis A. Simpson'
__version__ = '0.9.0'
__package__ = 'ScarMapper'


class TargetMapper:
    def __init__(self, log, args):
        """

        :param log:
        :param args:
        """
        self.target_data = Tool_Box.FileParser.indices(log, args.Target_File)
        # self.target_data = [x[0] for x in Tool_Box.FileParser.indices(log, args.Target_File)]
        self.refseq = pysam.FastaFile(args.Ref_Seq)
        self.log = log
        self.args = args

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

            # Remove the first 10 nt of the primer.  Allows the experimental reads to be trimmed and anchored to the
            # ends of target region easier.
            start_pos = int(line[2])+10
            stop_pos = int(line[3])-10

            sgrna_seq = line[4].upper()
            rcomp = line[5].upper()

            target_dict[locus_name] = \
                (locus_name, chromosome, start_pos, stop_pos, sgrna_seq, rcomp)

        return target_dict


