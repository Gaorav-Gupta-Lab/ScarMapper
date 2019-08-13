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
        target_dict = collections.defaultdict(list)
        for line in self.target_data:
            target_dict[line[4]] = line

        return target_dict


