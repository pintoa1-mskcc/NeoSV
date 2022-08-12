# NeoSV
[![Latest Version](https://pypip.in/version/<NeoSV>/badge.svg)](https://pypi.python.org/pypi/<NeoSV>/)
NeoSV is a computational workflow to identify neoantigens from structural variations (SVs). NeoSV has been sucessfully applied to 2,528 tumors from the Pan Cancer Analysis of Whole Genomes (PCAWG) project. For more detailed information, please see:. 

# Installation
### Prerequisites
* [Python (>3.6)](https://www.python.org/downloads/). NeoSV should work well with all versions of Python3, but has been only tested on Python > 3.6
* [netMHCpan (>4.0)](https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1). After you sign up and get the link for downloading, there will be a accompanied guidance on how to configure netMHCpan.
### Download
  We provide three ways for downloading and installing NeoSV:
  1. PyPI (recommended): if you already have Python, you can directly install it via `pip install NeoSV`<br>
  2. Anaconda (recommended): or you can install NeoSV and all dependencies via: `conda install -c bioconda NeoSV`<br>
  3. From source: download the code from github `python setup.py install`. If you do this way, you also need to install some dependencies, including: biopython, pyensembl.


# User guide
### Input
NeoSV requires 3 types of inputs:
* **_Mutation:_** a file in VCF format which lists all SVs you want to analyze. Please see sv.vcf as a template.
* **_HLA:_** a file listing the HLA types line-by-line. Usually this includes six HLA alleles for a patient. HLA should be in 4 digit format like: HLA-A*02:01. Please see hla.txt as a template.
* **_Reference:_** NeoSV utilizes pyensembl for SV annotation, thus a reference for pyensembl is needed. There are 3 ways to prepare it: <br>
  - **Pre-download by pyensembl (recommended):** When you install NeoSV using pip or conda, pyensembl will be automatically installed as well. Then you can download the reference:<br>
    ```
    export PYENSEMBL_CACHE_DIR=/custom/cache/dir # specify the location for storing reference
    pyensembl install --release <list of Ensembl release numbers> --species <species-name> # download, for hg19 please use release 75, for hg38 please used release 96
    ```
  - **Automatically download by NeoSV:** If NeoSV did not detect a valid reference in , it will automatically download one to that folder. However, please make sure that your server/computer can connect to the internet, because most high performance computing nodes are disconnected.
  - **Prepare the reference file manually:** This would be useful if your data is not from human or mouse. Then you need to prepare the reference by yourself. A FASTA file and a GTF file will be enough. For more details please see the guidance in pyensembl.
### Run

### Output

# Citation
