import collections

import Valkyries.Tool_Box as Tool_Box
import pysam


class TargetMapper:
    def __init__(self, log, args):

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
            chromosome = line[0]
            start_pos = int(line[1])
            stop_pos = int(line[2])
            locus_name = line[3]
            sgrna_seq = line[4].upper()
            rcomp = line[5].upper()

            target_dict[line[3]] = \
                (chromosome, start_pos, stop_pos, locus_name, sgrna_seq, rcomp)

        return target_dict


