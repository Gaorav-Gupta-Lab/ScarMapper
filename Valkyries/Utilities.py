"""
Odin specific utility functions
"""

import collections
from collections import Counter
import Valkyries.Tool_Box as Tool_Box
import natsort


class CountingGenerator:
    def __init__(self):
        self.item_count = 0

    def count(self, generator):
        for i in generator:
            self.item_count += 1
            yield i


def family_counts(region, ctg, family_dict, previous_exon_end, targeted_counts, index):

    for alignment in region:
        # Count reads on reverse strand.  Includes pairs when in a proper pair.
        if alignment.is_reverse:
            targeted_counts["reverse"] += 1
            # if alignment.is_proper_pair:
            #     targeted_counts["reverse"] += 1

        umt = "{}.{}|{}".format(ctg, alignment.reference_start, alignment.query_name.split(":")[0])

        if alignment.reference_start > previous_exon_end:
            family_dict[index].append(umt)
            targeted_counts["total"] += 1

    return family_dict, targeted_counts


class FamilyDistribution:
    """
    This gives an approximation of the duplication level of the reads.  Not a true count yet because it does not use
    read position.
    """
    def __init__(self, log, args, data_source):

        self._family_data = collections.defaultdict(list)
        if getattr(args, "Index_File", None) is not None:
            self._indices = Tool_Box.FileParser.indices(log, args.Index_File)

        self._args = args
        self._log = log
        self._data_source = data_source

    def data_collection(self, library_index, umi):
        self._family_data[library_index] += umi

    def data_processing(self):
        self._log.info("Begin Family Size and UMT Analysis")
        umt_stats_outstring = "UMT\tCount"
        family_size_outstring = "Family\tCount"

        for index_key in self._family_data:
            count_list = []

            for k, v in sorted(Counter(self._family_data[index_key]).items(), key=lambda x: x[1], reverse=True):
                umt_stats_outstring += "\n{0}\t{1}".format(k, v)
                count_list.append(str(v))

            c = dict(Counter(count_list))
            for k in natsort.natsorted(c):
                family_size_outstring += "\n{0}\t{1}".format(k, c[k])

        stats_filename = "{0}{1}_UMT_Stats.txt".format(self._args.Working_Folder, self._data_source)
        size_filename = "{0}{1}_Family_Size.txt".format(self._args.Working_Folder, self._data_source)

        # Deleting the files if they exist prevents a random text file busy OSError I am getting using VBox on Windows.
        Tool_Box.delete([stats_filename, size_filename])

        umt_stats_file = open("{0}{1}_UMT_Stats.txt".format(self._args.Working_Folder, self._data_source), 'w')

        family_size_file = open(size_filename, "w")

        umt_stats_file.write(umt_stats_outstring)
        family_size_file.write(family_size_outstring)

        umt_stats_file.close()
        family_size_file.close()

        self._log.info("{0} {1} UMT Family Size and Stats Files Written"
                       .format(self._args.Job_Name, self._data_source))


class FilteredGenerator:
    """Applies filters to a base collection yielding the item and its filter"""
    def __init__(self, filter_dict):
        """
        Args:
            filter_dict (dict): key = filter name, value = function that
                that accepts an item and returns true is that item should
                be excluded. For example: {"div by 2": lambda x: x % 2 == 0}
        """
        self._filters = sorted(filter_dict.items(), key=lambda x: x[0])

    def filter(self, base_collection):
        """Yields subset of base_collection/generator based on filters."""
        for item in base_collection:
            excluded = []
            for (name, exclude) in self._filters:
                if exclude(item):
                    excluded.append(name)
            if excluded:
                filter_value = "; ".join(excluded)
            else:
                filter_value = None
            yield item, filter_value


def byte_array_to_string(value):
    if isinstance(value, str):
        return value
    else:
        return str(value.decode("utf-8"))
