

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

tree = reorder(tree, order = "cladewise")
tree_postorder = reorder(tree, order = "postorder")
map_for_child = tree_postorder$edge[dim(tree_postorder$edge)[1]:1, ]

####

get_child_index = function(x, map_for_child, which = 0){

  num_leaf = dim(map_for_child)[1]/2 + 1
  if(x <= num_leaf){
    print("Warning: node has no child.")
    return(NA)
  }else if(x > 2 * num_leaf - 1){
    print("Error: node index out of bound.")
    return(NA)
  }else{
    if(which == 0){
      return(map_for_child[(x - num_leaf) * 2, 2])
    }else if(which == 1){
      return(map_for_child[(x - num_leaf) * 2 - 1, 2])
    }else{
      print("Error: 'which' should be 0/1.")
      return(NULL)
    }
  }
}

####

get_descendant = function(x, map_for_child){

  num_leaf = dim(map_for_child)[1]/2 + 1
  if(x <= num_leaf){
    return(x)
  }else{
    left = get_child_index(x, map_for_child, 0)
    right = get_child_index(x, map_for_child, 1)

    des = union(get_descendant(left, map_for_child), get_descendant(right, map_for_child))
    return(des)
  }
}


####

num_descendant = function(map_for_child){

  num_leaf = dim(map_for_child)[1]/2 + 1
  num_node = num_leaf - 1
  num = rep(1, num_leaf + num_node)
  for(i in (num_leaf + 1):(num_leaf + num_node)){
    num[i] = length(get_descendant(i, map_for_child))
  }
  return(num)
}


##################################################################################################
##################################################################################################

AG_simulate_senario_3 = function(seed, inc = c(1, 1.5, 2.5), save = TRUE){

  set.seed(seed)

  label = tree$tip.label
  otu_top = otu_top[ , otu_top[label[95], ] <= 500]
  otu_top = otu_top[ , otu_top[label[96], ] <= 500]
  otu_top = otu_top[ , otu_top[label[97], ] <= 500]
  otu_top = otu_top[ , otu_top[label[98], ] <= 500]

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
  sim_dir_null = 0

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

  for(k in 1:length(inc)){
    otu_group_2 = otu_top[, id_group_2]

    otu_group_2[label[97], ] = (1 + inc[k] * 1) * otu_group_2[label[97], ]
    otu_group_2[label[96], ] = (1 + inc[k] * 0.67) * otu_group_2[label[96], ]
    otu_group_2[label[95], ] = (1 + inc[k] * 0.33) * otu_group_2[label[95], ]

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
    file = paste("./Result_sim3/sim_", seed, ".RData", sep = "")
    save(res, file = file)
  }else{
    return(res)
  }
}

##################################################################################################
##################################################################################################

seed = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
AG_simulate_senario_3(seed, inc = c(0.75, 1, 1.25), save = FALSE)
