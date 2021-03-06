<p align="center">
  <img width="2000" src="/figures/logo.png">
</p>

# DeCatCounter

This is the README document for DeCatCounter, a python pipeline for processing concatenated PacBio reads from _in vitro_ selection experiments. The pipeline can demultiplex the amplicons, deconcatenatenate into sequence and count their count reads. DeCatCounter can be used to process nucleotides or amino acids sequencing data.

# Usage

To use the pipeline, type in the terminal:

`python DeCatCounter.py input_file barcodes.txt barc_tol_f barc_tol_r adapters.txt adap_tol_f adap_tol_r translation(y/n) low_len hi_len`

![sequences](/figures/sequences.png)

* input_file: name of input file (must include the full path to the directory where it's located).
* barcodes.txt: text file with 3 columns: 1) sample name, 2)corresponding forward barcode, 3) reverse barcode.
* bc_tol_f: error tolerance for forward barcode search in the forward reads(for the reverse reads, this is the error tolerance for the rC of the reverse barcode).
* bc_tol_r: error tolerance for reverse barcode search.
* adapters.txt: text file with 2 lines: 1) forward adapter, 2) reverse adapter.
* ad_tol_f: error tolerance for forward adapter search.
* ad_tol_r: error tolerance for reverse adapter search.
* translation(y/n): whether translation to amino acids should be performed, value should be either y or n. 
* low_len: minimum length for final DNA variants.
* hi_len: maximum length for final DNA variants.

# Dependencies
The pipeline script was written to run on Unix-based systems, like Linux, Ubuntu, and MacOS. Windows 10 also has a [Linux subsystem](https://docs.microsoft.com/en-us/windows/wsl/faq).

To use the pipeline, first install [Python](https://www.python.org/downloads/). We recommend using the Anaconda distribution of python, and adding the Bioconda channel to Anaconda's package manager, conda. See the [Anaconda documentation](https://docs.anaconda.com/anaconda/install/) for installation. 

The following dependencies are needed: [biopython](https://biopython.org/), [futures](https://docs.python.org/3/library/asyncio-future.html), [python-Levenshtein](https://pypi.org/project/python-Levenshtein/), [tabulate](https://pypi.org/project/tabulate/) and [pandas](https://pandas.pydata.org/). You can either install them one by one, or you can download the file `requirements.txt` and run the following command from the directory it's located:

`pip install -r requirements/txt` 

The pipeline will not be found unless it is stored in the working directory or in a directory that is in the user's PATH environment. For example, for Unix/Linux users, scripts could be placed in `/usr/local/bin/` upon download. These files can be placed in that directory with the command:

`cp /path/to/DeCatCounter.py /usr/local/bin/`

# Input

All input files (sequencing reads, barcodes and adapters text files) must be located in the same directory.
Sequencing reads are assumed to be in FASTA format. If your sequencing reads file is in FASTQ format, please convert it to fasta using the following command (convert input.fastq to output.fasta):

`cat input.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > output.fasta`

The barcodes files should be a text file with 3 columns: 1) sample name, 2)corresponding forward barcode, 3) reverse barcode. For example:

<p align="center">
  <img width="400" src="/figures/barcodes.png">
</p>

The adapter files should have 2 lines: 1) forward adapter, 2) reverse adapter. For example:

<p align="center">
  <img width="400" src="/figures/adapters.png">
</p>

# Output

The pipeline will generate an output directory, called `output+date&time`. Inside this directory, there will be a subdirectory for the count files and, if translation has been requested, a subdirectory for the amino acid count files. 

There pipeline also creates a log file with the parameters used and a table listing the number of amplicons recovered after demultiplexing, and the number of sequences recovered after deconcatenation and after trimming and filtering by length.

# Test dataset

A mock, test dataset (test_input.fasta) is provided, together with barcodes and adapter text files (barcodes.txt, adapters.txt).
To run the test dataset, type:

`python DeCatCounter.py test_input.fasta barcodes.txt 0 0 adapters.txt 0 0 y 5 50`

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

