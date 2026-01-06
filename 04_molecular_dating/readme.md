# Molecular dating of reference species phylogeny
## Overview

This repository implements a Bayesian molecular dating framework for the reference species phylogeny using PhyloBayes v4.1c, integrating:

1. Root-prior distribution (uniform vs. gamma)
2. Molecular clock model (autocorrelated and uncorrelated)
3. Alternative time constraints
4. Alternative backbone topologies (consensus, LG+C60+G, partition)

## Input data

Because PhyloBayes requires PHYLIP‐formatted alignments and taxon names shorter than 10 characters, all input files were reformatted and renamed prior to molecular dating. These preprocessing steps are purely technical and do not affect phylogenetic inference.

The input folder **input** contains:

1. ID_correspondence.csv
	- Mapping table linking original species names to PHYLIP‐compatible identifiers.
2. alignment.phylip
	- Multiple‐sequence alignment in PHYLIP format with renamed taxa.
3. consensus.topology.acc
	- Consensus species tree topology using PHYLIP‐compatible taxon IDs.
4. LG+C60+G.topology.acc
	- Alternative backbone topology inferred under the LG+C60+G model.
5. partition.topology.acc
	- Alternative backbone topology inferred from partitioned models.
6. scenario01 – scenario14
	- 14 calibration files defining the temporal constraints and root priors used in each dating scenario.
7. ITOL_label.txt
	- Tree label file.

## Dating

For each analytic scenario, two independent PhyloBayes Markov chain Monte Carlo (MCMC) chains were run in parallel, with parameters sampled every 20 cycles. The convergence of Phylobayes runs was assessed using ```tracecomp``` command, which monitors the minimum effective sampling size (ESS) and the maximum discrepancy (maxdiff) between independent MCMC chains across twelve parameters. Runs were considered converged when ESS > 50 and maxdiff < 0.3. 

Sometimes, satisfying both criteria was difficult even after extended runs (>40,000 cycles), a common issue in high-dimensional, large-scale phylogenomic analyses . In such cases, following standard practice, runs were considered acceptable when the estimated divergence times of the focal nodes were consistent between independent chains.

The convergence assessment was performed using:

```
tracecomp -x <burnin> scenario<id>_ch1 scenario<id>_ch2
```

### Primary analysis ( _scenario 1-4_ )
```
# Scenario 1
mkdir res_s01
pb -d input/alignment.phylip -t input/consensus.topology.acc -cal input/scenario01 -cat -c20 -dgam 4 -ln res_s01/scenario01_ch1
pb -d input/alignment.phylip -t input/consensus.topology.acc -cal input/scenario01 -cat -c20 -dgam 4 -ln res_s01/scenario01_ch2

# Scenario 2
mkdir res_s02
pb -d input/alignment.phylip -t input/consensus.topology.acc -cal input/scenario02 -cat -c20 -dgam 4 -ln -r 3950 230 res_s02/scenario02_ch1
pb -d input/alignment.phylip -t input/consensus.topology.acc -cal input/scenario02 -cat -c20 -dgam 4 -ln -r 3950 230 res_s02/scenario02_ch2

# Scenario 3
mkdir res_s03
pb -d input/alignment.phylip -t input/consensus.topology.acc -cal input/scenario03 -cat -c20 -dgam 4 -ugam res_s03/scenario03_ch1
pb -d input/alignment.phylip -t input/consensus.topology.acc -cal input/scenario03 -cat -c20 -dgam 4 -ugam res_s03/scenario03_ch2

# Scenario 4
mkdir res_s04
pb -d input/alignment.phylip -t input/consensus.topology.acc -cal input/scenario04 -cat -c20 -dgam 4 -ugam -r 3950 230 res_s04/scenario04_ch1
pb -d input/alignment.phylip -t input/consensus.topology.acc -cal input/scenario04 -cat -c20 -dgam 4 -ugam -r 3950 230 res_s04/scenario04_ch2
```

Based on the performance of the primary runs (**_ _scenario 1–4_ _**), the Uncorrelated Gamma (-ugam) model was selected for the remaining sensitivity analyses (**_ _scenario 5–12_ _**), because this model outperformed the autocorrelated alternative by providing more efficient MCMC sampling, higher ESS values, and more congruent age estimates between independent runs.

### Sensitivity analysis ( _scenario 5-14_ )
```
# Scenario 5
mkdir res_s05
pb -d input/alignment.phylip -t input/consensus.topology.acc -cal input/scenario05 -cat -c20 -dgam 4 -ugam res_s05/scenario05_ch1
pb -d input/alignment.phylip -t input/consensus.topology.acc -cal input/scenario05 -cat -c20 -dgam 4 -ugam res_s05/scenario05_ch2

# Scenario 6
mkdir res_s06
pb -d input/alignment.phylip -t input/consensus.topology.acc -cal input/scenario06 -cat -c20 -dgam 4 -ugam res_s06/scenario06_ch1
pb -d input/alignment.phylip -t input/consensus.topology.acc -cal input/scenario06 -cat -c20 -dgam 4 -ugam res_s06/scenario06_ch2

# Scenario 7
mkdir res_s07
pb -d input/alignment.phylip -t input/consensus.topology.acc -cal input/scenario07 -cat -c20 -dgam 4 -ugam res_s07/scenario07_ch1
pb -d input/alignment.phylip -t input/consensus.topology.acc -cal input/scenario07 -cat -c20 -dgam 4 -ugam res_s07/scenario07_ch2

# Scenario 8
mkdir res_s08
pb -d input/alignment.phylip -t input/consensus.topology.acc -cal input/scenario08 -cat -c20 -dgam 4 -ugam res_s08/scenario08_ch1
pb -d input/alignment.phylip -t input/consensus.topology.acc -cal input/scenario08 -cat -c20 -dgam 4 -ugam res_s08/scenario08_ch2

# Scenario 9
mkdir res_s09
pb -d input/alignment.phylip -t input/consensus.topology.acc -cal input/scenario09 -cat -c20 -dgam 4 -ugam res_s09/scenario09_ch1
pb -d input/alignment.phylip -t input/consensus.topology.acc -cal input/scenario09 -cat -c20 -dgam 4 -ugam res_s09/scenario09_ch2

# Scenario 10
mkdir res_s10
pb -d input/alignment.phylip -t input/consensus.topology.acc -cal input/scenario10 -cat -c20 -dgam 4 -ugam res_s10/scenario10_ch1
pb -d input/alignment.phylip -t input/consensus.topology.acc -cal input/scenario10 -cat -c20 -dgam 4 -ugam res_s10/scenario10_ch2

# Scenario 11
mkdir res_s11
pb -d input/alignment.phylip -t input/consensus.topology.acc -cal input/scenario11 -cat -c20 -dgam 4 -ugam res_s11/scenario11_ch1
pb -d input/alignment.phylip -t input/consensus.topology.acc -cal input/scenario11 -cat -c20 -dgam 4 -ugam res_s11/scenario11_ch2

# Scenario 12
mkdir res_s12
pb -d input/alignment.phylip -t input/consensus.topology.acc -cal input/scenario12 -cat -c20 -dgam 4 -ugam res_s12/scenario12_ch1
pb -d input/alignment.phylip -t input/consensus.topology.acc -cal input/scenario12 -cat -c20 -dgam 4 -ugam res_s12/scenario12_ch2

# Scenario 13
mkdir res_s13
pb -d input/alignment.phylip -t input/LG+C60+G.topology.acc -cal input/scenario13 -cat -c20 -dgam 4 -ugam res_s13/scenario13_ch1
pb -d input/alignment.phylip -t input/LG+C60+G.topology.acc -cal input/scenario13 -cat -c20 -dgam 4 -ugam res_s13/scenario13_ch2

# Scenario 14
mkdir res_s14
pb -d input/alignment.phylip -t input/partition.topology.acc -cal input/scenario14 -cat -c20 -dgam 4 -ugam res_s14/scenario14_ch1
pb -d input/alignment.phylip -t input/partition.topology.acc -cal input/scenario14 -cat -c20 -dgam 4 -ugam res_s14/scenario14_ch2
```

### Extract results
After convergence was reached, the ```readdiv``` command was used to extract divergence‐time estimates for each analytical scenario. The appropriate burn-in was chosen based on the convergence diagnostics.

**1. Generating the chronology**

This was performed using:
```
readdiv -x _<burnin> <every> <until> <chain>
```

The resulting output file, **<chainname>_sample.chronogram**, contains the chronology with the sampled divergence times across the posterior distribution. For each scenario, the is uploaded in corresponding folder. For each scenario, this file is stored in **output** folder.

**2. Summarizing posterior age distributions**

To visualize and summarize the posterior age distributions, we aggregated the divergence‐time estimates from the chronology files **at every cycle** between the burn-in and the specified end point using the script ```run_readdiv_r_analysis.py```. This procedure extracts and compiles posterior samples for the focal nodes from each chain over the defined interval, enabling quantitative comparison of divergence‐time estimates between independent runs.

```
python run_readdiv_r_analysis.py --start 250 -end 800 --ouput output/s1_ch1_out.txt res_s01/scenario01_ch1
python run_readdiv_r_analysis.py --start 800 -end 1250 --ouput output/s1_ch2_out.txt res_s01/scenario01_ch2
python run_readdiv_r_analysis.py --start 250 -end 800 --ouput output/s2_ch1_out.txt res_s02/scenario02_ch1
python run_readdiv_r_analysis.py --start 1650 -end 1950 --ouput output/s2_ch2_out.txt res_s02/scenario02_ch2
python run_readdiv_r_analysis.py --start 250 -end 800 --ouput output/s3_ch1_out.txt res_s03/scenario03_ch1
python run_readdiv_r_analysis.py --start 250 -end 800 --ouput output/s3_ch2_out.txt res_s03/scenario03_ch2
python run_readdiv_r_analysis.py --start 250 -end 800 --ouput output/s4_ch1_out.txt res_s04/scenario04_ch1
python run_readdiv_r_analysis.py --start 250 -end 800 --ouput output/s4_ch2_out.txt res_s04/scenario04_ch2
python run_readdiv_r_analysis.py --start 250 -end 800 --ouput output/s5_ch1_out.txt res_s05/scenario05_ch1
python run_readdiv_r_analysis.py --start 250 -end 800 --ouput output/s5_ch2_out.txt res_s05/scenario05_ch2
python run_readdiv_r_analysis.py --start 700 -end 1200 --ouput output/s6_ch1_out.txt res_s06/scenario06_ch1
python run_readdiv_r_analysis.py --start 700 -end 1200 --ouput output/s6_ch2_out.txt res_s06/scenario06_ch2
python run_readdiv_r_analysis.py --start 700 -end 1200 --ouput output/s7_ch1_out.txt res_s07/scenario07_ch1
python run_readdiv_r_analysis.py --start 700 -end 1200 --ouput output/s7_ch2_out.txt res_s07/scenario07_ch2
python run_readdiv_r_analysis.py --start 250 -end 800 --ouput output/s8_ch1_out.txt res_s08/scenario08_ch1
python run_readdiv_r_analysis.py --start 250 -end 800 --ouput output/s8_ch2_out.txt res_s08/scenario08_ch2
python run_readdiv_r_analysis.py --start 700 -end 1200 --ouput output/s9_ch1_out.txt res_s09/scenario09_ch1
python run_readdiv_r_analysis.py --start 700 -end 1200 --ouput output/s9_ch2_out.txt res_s09/scenario09_ch2
python run_readdiv_r_analysis.py --start 1300 -end 1700 --ouput output/s10_ch1_out.txt res_s10/scenario10_ch1
python run_readdiv_r_analysis.py --start 700 -end 1200 --ouput output/s10_ch2_out.txt res_s10/scenario10_ch2
python run_readdiv_r_analysis.py --start 250 -end 800 --ouput output/s11_ch1_out.txt res_s11/scenario11_ch1
python run_readdiv_r_analysis.py --start 250 -end 800 --ouput output/s11_ch2_out.txt res_s11/scenario11_ch2
python run_readdiv_r_analysis.py --start 250 -end 800 --ouput output/s12_ch1_out.txt res_s12/scenario12_ch1
python run_readdiv_r_analysis.py --start 250 -end 800 --ouput output/s12_ch2_out.txt res_s12/scenario12_ch2
python run_readdiv_r_analysis.py --start 250 -end 800 --ouput output/s11_ch1_out.txt res_s13/scenario13_ch1
python run_readdiv_r_analysis.py --start 250 -end 800 --ouput output/s11_ch2_out.txt res_s13/scenario13_ch2
python run_readdiv_r_analysis.py --start 550 -end 1300 --ouput output/s14_ch1_out.txt res_s14/scenario14_ch1
python run_readdiv_r_analysis.py --start 550 -end 1300 --ouput output/s14_ch2_out.txt res_s14/scenario14_ch2
```

