<p align="center">
  <img width="2000" src="/figures/logo.png">
</p>

# DeCatCounter

This is the README document for DeCatCounter, a Python pipeline for processing concatenated PacBio reads from _in vitro_ selection experiments. The pipeline can demultiplex the amplicons, deconcatenatenate into sequence and count their count reads. DeCatCounter can be used to process nucleotides or amino acids sequencing data.

# Usage

The pipeline tool is written in Python and it can be run using a Python interpreter. To run it from the Command Line Interface (aka the Terminal), type:

`python DeCatCounter.py input_file barcodes.txt bc_tol_f bc_tol_r constant.txt ct_tol_f ct_tol_r translation(y/n) low_len hi_len`

![sequences](/figures/sequences.png)

* input_file: name of input file (must include the full path to the directory where it's located).
* barcodes.txt: text file with 3 columns: 1) sample name, 2)corresponding forward barcode, 3) reverse barcode.
* bc_tol_f: error tolerance for forward barcode search in the forward reads(for the reverse reads, this is the error tolerance for the rC of the reverse barcode).
* bc_tol_r: error tolerance for reverse barcode search.
* constant.txt: text file with 2 lines: 1) forward constant region, 2) reverse constant region.
* ct_tol_f: error tolerance for forward constant region search.
* ct_tol_r: error tolerance for reverse constant region search.
* translation(y/n): whether translation to amino acids should be performed, value should be either y or n. 
* low_len: minimum length for final DNA variants.
* hi_len: maximum length for final DNA variants.

# Environment setup

Python must be installed to run DeCatCounter. We recommend using the Anaconda distribution of Python, and adding the Bioconda channel to Anaconda's package manager, Conda. See the [Anaconda documentation](https://docs.anaconda.com/anaconda/install/) for installation. 

The following dependencies are needed: [biopython](https://biopython.org/), [python-Levenshtein](https://pypi.org/project/python-Levenshtein/), [tabulate](https://pypi.org/project/tabulate/) and [pandas](https://pandas.pydata.org/). You can either install them one by one, or you can download the file `requirements.txt` and run the following command from the directory it's located:

`conda install --file requirements.txt` 

Although not necessary, we recommend [using Anaconda to create a virtual environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) where you can install the required dependencies. To do this:

`conda create -n py3 python=3`
`source activate py3`

# Input

All input files (sequencing reads, barcodes and constant regions text files) must be located in the same directory.
Sequencing reads are assumed to be in FASTA format. If your sequencing reads file is in FASTQ format, please convert it to fasta using the following command (convert input.fastq to output.fasta):

`cat input.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > output.fasta`

The barcodes files should be a text file with 3 columns: 1) sample name, 2)corresponding forward barcode, 3) reverse barcode. For example:

<p align="center">
  <img width="400" src="/figures/barcodes.png">
</p>

The constant regions files should have 2 lines: 1) forward constant region, 2) reverse constant region. For example:

<p align="center">
  <img width="350" src="/figures/constant.png">
</p>

# Output

The pipeline will generate an output directory, called `output+date&time`. Inside this directory, there will be a subdirectory for the count files and, if translation has been requested, a subdirectory for the amino acid count files. 

There pipeline also creates a log file with the parameters used and a table listing the number of amplicons recovered after demultiplexing, and the number of sequences recovered after deconcatenation and after trimming and filtering by length.

# Test dataset

A mock, test dataset (test_input.fasta) is provided, together with barcodes and constant regions text files (barcodes.txt, constant.txt).
To run the test dataset, place all files (DeCatCounter.py, test_input.fasta, barcodes.txt and constant.txt) in the same folder and type:

`python DeCatCounter.py test_input.fasta barcodes.txt 0 0 constant.txt 0 0 y 5 50`

If everything is working, your terminal should look like this:

<p align="center">
  <img width="400" src="/figures/test_terminal.png">
</p>

and you should have a new folder in your directory:

<p align="center">
  <img width="1000" src="/figures/output.png">
</p>
   
# Reporting bugs

Please report any bugs to Celia Blanco (celiablanco@ucla.edu). 

When reporting bugs, please include the full output printed in the terminal when running the pipeline. 

# Citation

Nisha Kanwar<sup>\*</sup>, Celia Blanco<sup>\*</sup>, Irene A. Chen and Burckhard Seelig. PacBio sequencing output increased through uniform and directional 5-fold concatenation. *Submitted.*

