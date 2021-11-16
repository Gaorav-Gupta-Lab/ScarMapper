"""
Setup file to Cythonize Alignment Processing using "python3 setup.py build_ext --inplace"
"""
import os
import pathlib
from distutils.core import setup
from Cython.Build import cythonize

slidingwindow_file = '{}{}SlidingWindow.pyx'.format(pathlib.Path(__file__).parent.absolute(), os.sep)

setup(
    name="ScarMapper Sliding Window",
    author='Dennis Simpson',
    author_email='dennis@email.unc.edu',
    ext_modules=cythonize(slidingwindow_file, annotate=False)
)
