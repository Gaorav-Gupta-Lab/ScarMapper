# ScarMapper
## A Python-encoded algorithm that uses an iterative break-associated alignment strategy to classify individual double-strand DNA break repair products based on deletion size, microhomology usage, and insertions

### Getting Started

These instructions should allow a functional copy of this package to be installed on your local computer.  It is highly recomended
that you read the user guide before attempting to use this program.

### Prerequisites
Linux OS, tested on RHEL 7.x, Scientific Linux 7.x, and CentOS 7.x.  Will possibly run on a Mac OS, although it has not been tested.  Will not run on Windows because ScarMapper uses Pysam to parse the reference FASTA file.

```
Minimum System Requirements:
    4 CPUs or threads
    20 Gb RAM
    ~5x the size of the compressed FASTQ file for disk space.

Required Python Libraries:
    Python ≥v3.5, ≥v3.6 recomended
    numpy
    scipy
    python-levenshtein
    python-magic
    natsort
    pathos
    pysam
    cython
    setuptools
```
### Installation

The quickest way to install is to clone or download this repository into a location that accessible to the user.
Test installation by moving to the ScarMapper directory and executing ```python3 scarmapper.py```.  You should get the error message ```usage: scarmapper.py [-h] --options_file OPTIONS_FILE
 scarmapper.py: error: the following arguments are required: --options_file```.
 
### Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

### Authors

* **Dennis Simpson** - *Initial work* 

### Cite

[![PubMed](img/2318832.png)](https://www.ncbi.nlm.nih.gov/pubmed/xxx)

### License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
