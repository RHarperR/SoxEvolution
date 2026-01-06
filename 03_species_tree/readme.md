# Species tree reconstruction
## Sequence Alignment & Trimming
Ribosomal proteins and RNA polymerase subunits were identified and aligned using markerfinder pipeline (Martinez-Gutierrez and Aylward, 2021).

```
trimal -in 03_ribo_rnap.aln -out 03_ribo_rnap.aln.trim -gt 0.1
```

## Phylogenetic Reconstruction



### LG+R10 tree
A first-pass species tree was inferred using the LG+R10 model in IQ-TREE, which accounts for among-site rate heterogeneity. Notably, the LG+R10 model is also recommended by Martinez-Gutierrez and Aylward (2021) who use the same gene set.

	```
	iqtree -s 03_ribo_rnap.aln.trim -m LG+R10 --prefix 04_species_LG_R10 -T 20
	```
### LG+C60+G tree
To better account for compositional heterogeneity across sites, we inferred a second tree under the LG+C60+G mixture model.

The LG+R10 tree was supplied as a fixed starting tree (-ft) to speed up convergence.

	```
	iqtree -s 03_ribo_rnap.aln.trim -m LG+C60+G -ft 04_species_LG_R10.treefile --prefix species_LG_R10 -T 20
	```

### partition tree

To account for gene-wise heterogeneity, the alignment was partitioned by marker genes.
Each partition was assigned its own best-fit substitution model, following Martinez-Gutierrez and Aylward (2021).

	```
	iqtree -s 03_ribo_rnap.aln.trim -p partition.nex -nt 30 -seed 123 --prefix 04_species_partition -T 20
	```

### consensus tree

The consensus topology **04_species_consensus.topology** was pruned and regraphed from the LG+R10 tree to ensure the monophyly of (1) the Archaea, (2) the DST group comprising phylum Deinococcota, Synergistota, and Thermotogota, (3) the Terrabacteria superphylum, and (4) the Gracilicutes superphylum. This topology represents the final working hypothesis of species relationships.

The consensus topology was used as a fixed tree to re-estimate branch lengths and
bootstrap support values under the LG+R10 model.

	```
	iqtree -s 03_ribo_rnap.aln.trim -t 04_species_consensus.topology -m LG+R10 --prefix 04_species_consensus -nt 20 -B 1000
	```

## topology test between the LG+R10 and consensus treefile

Topology tests were performed to evaluate whether the LG+R10 tree can be statistically rejected in favor of the consensus topology.

Rooted trees were concatenated into a single tree set and tested.

	```
	cat species_LG_R10.treefile.rooted species_consensus.treefile.rooted > 05_tree_test_all
	iqtree -s 03_ribo_rnap.aln.trim -z 05_tree_test_all -m LG+R10 -n 0 -zb 1000 -au -T 15 --prefix 06_topotest
	```



