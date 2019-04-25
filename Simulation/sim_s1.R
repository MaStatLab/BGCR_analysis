
library(HMP)
library(cubature)
library(R2Cuba)
library(hashmap)
library(doParallel)
library(geiger)
library(ape)

source("ntaxa.R")
source("phyloscan.R")
source("phylostructure.R")
source("nodetest.R")
source("BGCR.R")
load("AG.RData")

##################################################################################################
##################################################################################################

AG_simulate_senario_1 = function(seed, inc = c(1, 1.5, 2.5), save = TRUE){

  set.seed(seed)

  CAUCASIAN = subset(ag_fecal, ag_fecal$RACE == "Caucasian" & ag_fecal$SEX == "male")
  WEST = subset(CAUCASIAN, CAUCASIAN$CENSUS_REGION == "West")
  WEST = subset(WEST, WEST$AGE_CAT == "30s" | WEST$AGE_CAT == "40s"| WEST$AGE_CAT == "50s")
  data_id = colnames(otu_top)[colnames(otu_top) %in% WEST$SampleID]
  n = length(data_id)
  p = dim(otu_top)[1]

  num_leaf = tree$Nnode
  num_node = tree$Nnode + 1

  beta_null = 0
  sim_PJAP_null = 0
  sim_PMAP_null = 0
  sim_PJAPind_null = 0
  sim_PMAPind_null = 0
  sim_w_null = 0
  sim_w_null_ind = 0
  sim_dir_null

  beta_alt = rep(0, length(inc))
  sim_PJAP_alt = rep(0, length(inc))
  sim_PMAP_alt = rep(0, length(inc))
  sim_PJAPind_alt = rep(0, length(inc))
  sim_PMAPind_alt = rep(0, length(inc))
  sim_w_alt = rep(0, length(inc))
  sim_w_alt_ind = rep(0, length(inc))
  sim_dir_alt = rep(0, length(inc))

    id_group_1 = sample(data_id, n / 2, replace = FALSE)
    id_group_2 = setdiff(data_id, id_group_1)

    otu_group_1 = otu_top[, id_group_1]
    otu_group_2 = otu_top[, id_group_2]

    res = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                            otu_group_1, otu_group_2, X_group_1 = "default", X_group_2 = "default",
                            nu = 10 ^ (seq(-1, 4, 0.25)), sigma = sqrt(10), verbose = TRUE)

    sim_PJAP_null = res$PJAP
    sim_PMAP_null = max(res$PMAP)
    sim_PJAPind_null = res$PJAP_ind
    sim_PMAPind_null = max(res$PMAP_ind)
    beta_null = res$beta

    pstrct = phylostructure(tree)
    group.data = list(x1 = t(otu_group_1), x2 = t(otu_group_2))
    nt = nodetest(pstrct, group.data) #Generate triplet statistics
    sim_w_null = nt$w
    sim_w_null_ind = min(nt$MoMp)
    sim_dir_null = Xmcupo.sevsample(group.data)$`p value`

    x1 = rep(0, p)
    x2 = rep(0, p)

    for(i in 1:p){
      x1[i] = sum(otu_group_1[i, ] == 0)
      x2[i] = sum(otu_group_2[i, ] == 0)
    }

    otusum_1 = apply(otu_group_1, 1, mean)
    otusum_2 = apply(otu_group_2, 1, mean)

    pool = which(x1 < n/8 & x2 < n/8)
    pool = intersect(pool, which(otusum_1 >= quantile(otusum_1, 0.1) & otusum_1 <= quantile(otusum_1, 0.9)))
    pool = intersect(pool, which(otusum_2 >= quantile(otusum_2, 0.1) & otusum_2 <= quantile(otusum_2, 0.9)))
    signal = sample(pool, 1)

    for(k in 1:length(inc)){
      otu_group_2 = otu_top[, id_group_2]
      otu_group_2[signal, ] = round((1 + inc[k]) * otu_group_2[signal, ])

      res = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                              otu_group_1, otu_group_2, X_group_1 = "default", X_group_2 = "default",
                              nu = 10 ^ (seq(-1, 4, 0.25)), sigma = sqrt(10))

      sim_PJAP_alt[k] = res$PJAP
      sim_PMAP_alt[k] = max(res$PMAP)
      sim_PJAPind_alt[k] = res$PJAP_ind
      sim_PMAPind_alt[k] = max(res$PMAP_ind)
      beta_alt[k] = res$beta

      group.data = list(x1 = t(otu_group_1), x2 = t(otu_group_2))
      nt = nodetest(pstrct, group.data) #Generate triplet statistics
      sim_w_alt[k] = nt$w
      sim_w_alt_ind[k] = min(nt$MoMp)
      sim_dir_alt[k] = Xmcupo.sevsample(group.data)$`p value`
    }

  res = c(sim_PJAP_null = sim_PJAP_null,
	        sim_PJAP_alt = sim_PJAP_alt,
          sim_PJAPind_null = sim_PJAPind_null,
          sim_PJAPind_alt = sim_PJAPind_alt,
          sim_w_null = sim_w_null, sim_w_alt = sim_w_alt,
	        sim_w_null_ind = sim_w_null_ind, sim_w_alt_ind = sim_w_alt_ind,
          sim_dir_null = sim_dir_null, sim_dir_alt = sim_dir_alt,
	        beta_null = beta_null, beta_alt = beta_alt)

  if(save == TRUE){
    file = paste("./Result_sim1/sim_", seed, ".RData", sep = "")
    save(res, file = file)
  }else{
    return(res)
  }
}

##################################################################################################
##################################################################################################

seed = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
AG_simulate_senario_1(seed, inc = c(1, 1.5, 2.5), save = TRUE)

