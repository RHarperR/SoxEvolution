# General description
The custom scripts here aim to retrieve gene cluster involved in dissimilatory sulfur oxidation from the GTDB genome database (https://gtdb.ecogenomic.org/). The hidden markov models (HMMs) targeting each of the sulfur oxidizing genes are compiled in the `00_hmm_model` folder. A subset of sequences (~30 M unzipped; GTDB_r214_demo.faa.gz) from the full GTDB r214 database (~51 G) was used as a demo dataset.


# System requirement
## Hardware requirements

The `hmmsearch` script and python script in finding the gene cluster require a standard computer with >4G RAM and >5 CPUs to run the Demo dataset, and a server with >10G RAM and >20 CPUs to run the full GTDB r214 dataset. 

## Software requirements

### OS Requirements
This package is supported for macOS and Linux. The script has been tested on the following systems:

- macOS: Mojave (10.14.6)
- Linux: Ubuntu 18.04.3 LTS

### Dependencies

The `hmmsearch` command mainly depends on the HMMER package (v3.2; http://hmmer.org/). The python scripts in searching of sulfur-oxidizing gene clusters mainly depend on the pandas package of python.

## Installation Guide:
### install HMMER and BLAST+ using conda

The installation of hmmer/BLAST+ takes less than 5 minutes with a good connection to the conda mirror.

conda install hmmer
conda install blast
pip install pandas




## 1 Identifying functional oxidative sulfur metabolic gene clusters
### 1.1 Homology Search with HMMER
To find Sox, rDsr, and sHdr homologs, we first search against the GTDB v214 using HMM profiles listed in 00_hmm_model with following command. Threshold are listed in Table S3.
	
	hmmsearch --tblout <output_file> --noali <Threshold> <HMM_profile> <protein_database>

The hmmsearch results were then processed into a table containing protein ID and name in 01_gene_cluster
	·01_sox_all.txt
	·0

### 1.2 Gene Cluster Identification
The sox/rDsr/sHdr gene clusters are identified using the following command.
	
	python code_soxcluster_find.py
	python code_rdsrcluster_find.py
	python code_shdrcluster_find.py

The retrieved sox results are listed as follows.
	1. An overview of sox gene clusters
	·01_sox_all.overview


