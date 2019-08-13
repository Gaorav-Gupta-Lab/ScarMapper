"""
Options_Parser.py
Version 0.2
    by Dennis A. Simpson
    July 1, 2015
    Convert output to namedtuple from dictionary.

Version 0.1
    by Dennis A. Simpson
    March 31, 2015
    Python class to parse options_file used by all duplex sequencing files.
"""

import re
import csv
import os
from collections import namedtuple

__author__ = 'Dennis A. Simpson'
__package__ = 'Mimir Sequence Analysis'
__version__ = '0.2'


class OptionsFileParser:
    """
    Class to parse options file for duplex sequence analysis package.
    """
    def __init__(self, o):
        self.options_file = o.options_file

    def options_file(self):
        """
        This function parses the file and returns an object.
        :return:
        """
        count = 0
        options_dictionary = {}  # Define dictionary for options file.
        if not os.path.isfile(self.options_file):  # Make sure this really is a file.
            exit("WARNING:\n\tOptions_File Not Found.  Check File Name and Path.")

        options = csv.reader(open(self.options_file), delimiter='\t')

        for line in options:
            count += 1
            if len(line) < 2 or "#" in str(line[0]):  # Strip out any lines that are comments or blank.
                line.clear()
            else:
                try:
                    value = re.sub(" |#.*| ", '', line[1]).strip()  # Strip out any end of line comments and whitespace.

                except IndexError:
                    exit("There is a syntax error in the options file on line "+str(count))

                key = line[0].strip('--')
                if key == 'index':
                    value = value.split(',')
                options_dictionary[key] = value

        return namedtuple('options_file', options_dictionary.keys())(**options_dictionary)
