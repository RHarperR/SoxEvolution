# SoxEvolution
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


