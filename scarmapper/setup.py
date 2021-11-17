"""
Setup file to Cythonize Alignment Processing using "python3 setup.py build_ext --inplace"
"""
import os
import pathlib
from shutil import rmtree, copyfile
from distutils.core import setup
from Cython.Build import cythonize

setup(
    name="scarmapper",
    author='Dennis Simpson',
    author_email='dennis@email.unc.edu',
    ext_modules=cythonize("SlidingWindow.pyx", annotate=False)
)

dst_dir = '{}{}'.format(pathlib.Path(__file__).parent.absolute(), os.sep)
src_dir = '{0}{1}scarmapper{1}'.format(pathlib.Path(__file__).parent.absolute(), os.sep)
files = os.listdir(src_dir)

for file in files:
    file_src = "{}{}".format(src_dir, file)
    file_dst = "{}{}".format(dst_dir, file)
    copyfile(file_src, file_dst)

rmtree(src_dir)
