# Bayesian Graphical Compositional Regression --- Simulation and Application

This folder contains the code for the simulation studies and the applications in Mao et, al (2018). The examples in the paper are generated under `R` version 3.4.3 and `Rcpp` version 0.12.13 (2017-09-24). Note that the working directory are default to be each subfolder when running code inside the subfolder.

## Data process
The folder "./Data_process" contains the raw data and the file to process the raw data. The file "97_otus.tree" contains the phylogenetic tree on the 27774 OTUs in the study. The file "ag_fecal.txt" contains the covariate information on each sample. The file "ag_fecal_from_biom.txt" contains the OTU table on all the OTUs. Running the file "Data_process.R" create a dataset containing the OTU table on the top 100 OTUs and the corresponding phylogenetic tree used in the application section. This file can be simply modified to create dataset containing the top 50 and 75 OTUs used for the simulation studies. 

## Simulation
The folder "./Simulation" contains the file used to generate data for the simulation studies. `R` script "sim_s1.R", "sim_s2.R", "sim_s3.R" and "sim_s4.R" perform one round of each simulation scenario, respectively. Running the first three scripts 3000 times with seed 1-3000 and the last script 750 times with seed 1-750 provides the data used for the figures in the simulation section (Section 3) of the paper. The files "nodetest.R", "ntaxa.R", "phyloscan.R" and "phylostructure.R" are used to implement the PhyloScan test in Tang et, al (2018). We appreciate the authors for sharing their code. Typically, it should take less than 3 hours to finish a single simulation round under a certain simulation scenario (we tested the code on a machine with 16GB memory and 2.5 GHz Intel Core i7 processor).

## Application
The folder "./Application" contains the code used to generate results in the Appication section. The file "AG_application.R" applies BGCR to compare groups of samples defined by different dietary habits. Directly running this file will create a series of `R` datasets (each dataset for each dietary habit) that can be used to create the plots and tables in the application section. Computation time for different dietary habits various due to the different sample size. Typically, the code should be able to finish within 5 hours for each diet (we run the analysis on a machine with 16GB memory and 2.5 GHz Intel Core i7 processor).

## Reference

Mao J., Chen Y. and Ma L. (2018). Bayesian graphical compositional regression for microbiome data. https://arxiv.org/abs/1712.04723.

Tang, Y., Ma, L. and Nicolae, D.L., (2018). A phylogenetic scan test on a dirichlet-tree multinomial model for microbiome data. The Annals of Applied Statistics, 12(1), pp.1-26.
