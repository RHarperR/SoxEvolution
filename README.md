# SoxEvolution
## Raw data
### hmmsearch against the GTDB v214
To find Sox, rDsr, and sHdr homologs, we first search against the GTDB v214 using HMM profiles listed in 00_hmm_model with following command.  

	```
	hmmsearch --tblout <output_file> --noali <Threshold> <HMM_profile> <protein_database>
	```
