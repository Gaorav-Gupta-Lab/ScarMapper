from setuptools import setup, find_namespace_packages
import sys

if sys.version_info < (3, 5):
    sys.stdout.write("At least Python 3.5 is required.\n")
    sys.exit(1)

setup(
    name='ScarMapper',
    version='0.1.0',
    packages=find_namespace_packages(include=['ScarMapper.*']),

    url='',
    license='MIT',
    author='Dennis Simpson',
    author_email='dennis@email.unc.edu',
    long_description=open('README.md').read(),

    description='Package for mapping repair scars in genomic DNA',
    install_requires=['scipy', 'natsort', 'pysam', 'python-magic', 'pathos', 'numpy']
)