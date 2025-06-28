# General description
The custom scripts here aim to retrieve gene clusters involved in dissimilatory sulfur oxidation from the GTDB genome database (https://gtdb.ecogenomic.org/). The hidden markov models (HMMs) targeting each of the sulfur oxidizing genes are compiled in the `00_hmm_model` folder. A subset of sequences (~22.2 M unzipped; GTDB_r214_demo.faa) from the full GTDB r214 database (~51 G) was used as a demo dataset.


# System requirement
## Hardware requirements

The `hmmsearch` script and python script in finding the gene cluster require a standard computer with >4G RAM and >5 CPUs to run the Demo dataset, and a server with >10G RAM and >20 CPUs to run the full GTDB r214 dataset. 

## Software requirements

### OS Requirements
This package is supported for macOS and Linux. The script has been tested on the following systems:

- macOS: Mojave (10.14.6)
- Linux: Ubuntu 18.04.3 LTS

### Dependencies

The `hmmsearch` command mainly depends on the HMMER package (v3.2; http://hmmer.org/). The python scripts in searching of sulfur-oxidizing gene clusters mainly depend on build-in modules including `os`, `re`, `sys`, and `collections`.

## Installation Guide:
### install HMMER and python3 using conda

The installation of hmmer/BLAST+ takes less than 5 minutes with a good connection to the conda mirror.
```
conda install hmmer=3.2 python=3.9
```


# Run Demo
### Step 1: Search for individual genes involved in dissimilatory sulfur oxidation using HMMER
To find Sox, rDsr, and sHdr homologs, we first search the HMM profiles listed in `00_hmm_model` against the demo dataset `GTDB_r214_demo.faa` using the `code_hmmsearch.sh` command. The results of hmmsearch were further parsed using the `code_parsehmm.py`. The output of `code_hmmsearch.sh` can be found in the directory `01_hmmsearch_res/` and the output of `code_parsehmm.py` can be found in the directory `02_gene_cluster`.
```
git clone https://github.com/RHarperR/SoxEvolution.git
cd SoxEvolution/
rm -r 01_hmmsearch_res
rm -r 02_gene_cluster

mkdir 01_hmmsearch_res
mkdir 02_gene_cluster
# the hmmsearch step takes ~5 min
source code_hmmsearch.sh &>/dev/null
python code_parsehmm.py 01_hmmsearch_res/ 02_gene_cluster/
```
### Step 2: identification of gene clusters based on homology search results
The sox/rDsr/sHdr gene clusters are identified using the following command. The gene clusters found in genomes will be listed in the directory `02_gene_cluster/`.
```
python code_rdsr_clusterfind.py 02_gene_cluster/00_rdsr_all.txt
python code_soxcluster_find.py 02_gene_cluster/00_sox_all.txt
```


