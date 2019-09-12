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
            genomic_cut_position = int(line[3])
            locus_name = line[4]
            sgrna_seq = line[5].upper()
            rcomp = line[6].upper()

            target_dict[line[4]] = \
                (chromosome, start_pos, stop_pos, genomic_cut_position, locus_name, sgrna_seq, rcomp)

        return target_dict


