"""
Version_Dependencies.py
Version 0.1.0
    Dec. 1, 2016
    Dennis A. Simpson
    This file contains the version numbers for the min and max dependencies of various modules.
"""
import platform
import pysam

from distutils.version import StrictVersion

__author__ = 'Dennis A. Simpson'
__version__ = '0.1.0'


def pysam_check():
    pysam_max = '0.16.0'
    pysam_min = '0.14.0'

    if StrictVersion(pysam.__version__) < pysam_min or StrictVersion(pysam.__version__) >= pysam_max:
        raise SystemExit("\033[1;31mERROR:\033[m Pysam Version MUST be > v{} and < v{}.  Current Pysam Version is {}"
                         .format(pysam_min, pysam_max, pysam.__version__))


def python_check():
    python_max = '3.8.0'
    python_min = '3.5.0'
    if StrictVersion(platform.python_version()) < python_min or StrictVersion(platform.python_version()) > python_max:
        raise SystemExit("\033[1;31mERROR:\033[m Python version must be > v{} and < v{}.  Current Python version is {}"
                         .format(python_min, python_max, platform.python_version()))
