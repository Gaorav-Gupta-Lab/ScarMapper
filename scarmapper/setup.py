"""
Setup file to Cythonize Alignment Processing using "python3 setup.py build_ext --inplace"
"""
import os
from distutils.core import setup
from Cython.Build import cythonize

slidingwindow_file = '{0}{1}SlidingWindow.pyx'.format(os.path.dirname(__file__), os.sep)

setup(
    name="ScarMapper Sliding Window",
    author='Dennis Simpson',
    author_email='dennis@email.unc.edu',
    ext_modules=cythonize(slidingwindow_file, annotate=True)
)
