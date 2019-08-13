"""
Setup file to Cythonize Alignment Processing using "python3 setup.py build_ext --inplace"
"""
from distutils.core import setup
from Cython.Build import cythonize

setup(
    name="ScarMapper Sliding Window",
    ext_modules=cythonize("SlidingWindow.pyx", annotate=True)
)
