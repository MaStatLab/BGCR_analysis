

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
load("AG_50.RData")


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

AG_simulate_senario_cov = function(seed, inc = c(1, 1.5, 2.5), inc_cov = 0.75, save = TRUE){

  set.seed(seed)

  label = tree$tip.label
  otu_top = otu_top[ , otu_top[label[47], ] <= 500]
  otu_top = otu_top[ , otu_top[label[37], ] <= 500]

  CAUCASIAN = subset(ag_fecal, ag_fecal$RACE == "Caucasian")
  WEST = subset(CAUCASIAN, CAUCASIAN$CENSUS_REGION == "West")
  WEST = subset(WEST, WEST$AGE_CAT == "30s" | WEST$AGE_CAT == "40s"| WEST$AGE_CAT == "50s")
  WEST = subset(WEST, WEST$SampleID %in% colnames(otu_top))
  WEST_male_id = sample(WEST[WEST$SEX == "male", ]$SampleID, 200, replace = FALSE)
  WEST_female_id = sample(WEST[WEST$SEX == "female", ]$SampleID, 200, replace = FALSE)
  data_id = colnames(otu_top)[colnames(otu_top) %in% c(WEST_male_id, WEST_female_id)]

  id_group_1 = c(sample(WEST_male_id, 160, replace = FALSE), sample(WEST_female_id, 40, replace = FALSE))
  id_group_2 = setdiff(data_id, id_group_1)

  id_group_1 = id_group_1[id_group_1 %in% data_id]

  p = dim(otu_top)[1]
  n = 400
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

  beta_null_cov = 0
  sim_PJAP_null_cov = 0
  sim_PJAPind_null_cov = 0


  beta_alt = rep(0, length(inc))
  sim_PJAP_alt = rep(0, length(inc))
  sim_PMAP_alt = rep(0, length(inc))
  sim_PJAPind_alt = rep(0, length(inc))
  sim_PMAPind_alt = rep(0, length(inc))
  sim_w_alt = rep(0, length(inc))
  sim_w_alt_ind = rep(0, length(inc))
  sim_dir_alt = rep(0, length(inc))

  beta_alt_cov = rep(0, length(inc))
  sim_PJAP_alt_cov = rep(0, length(inc))
  sim_PJAPind_alt_cov = rep(0, length(inc))

  ####

  otu_group_1 = otu_top[, id_group_1]
  otu_group_2 = otu_top[, id_group_2]

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

  cov1 = matrix(0, ncol = 1, nrow = length(id_group_1))
  cov2 = matrix(0, ncol = 1, nrow = length(id_group_2))

  for(i in 1:dim(cov1)[1]){
    if(WEST[ WEST$SampleID == id_group_1[i], ]$SEX == "male"){
      cov1[i, 1] = 1
    }
  }
  for(i in 1:dim(cov2)[1]){
    if(WEST[WEST$SampleID == id_group_2[i], ]$SEX == "male"){
      cov2[i, 1] = 1
    }
  }

  for(j in 1:length(cov1)){
    otu_group_1[label[47], j] = round((1 + inc_cov * cov1[j]) * otu_group_1[label[47], j])
  }

  for(j in 1:length(cov2)){
    otu_group_2[label[47], j] = round((1 + inc_cov * cov2[j]) * otu_group_2[label[47], j])
  }

  res = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                          otu_group_1, otu_group_2, X_group_1 = "default", X_group_2 = "default",
                          nu = 10 ^ (seq(-1, 4, 0.25)), sigma = sqrt(10), verbose = FALSE)

  sim_PJAP_null = res$PJAP
  sim_PJAPind_null = res$PJAP_ind
  beta_null = res$beta

  pstrct = phylostructure(tree)
  group.data = list(x1 = t(otu_group_1), x2 = t(otu_group_2))
  nt = nodetest(pstrct, group.data) #Generate triplet statistics
  sim_w_null = phyloscan(pstrct, nt$w, nthread = 2, gridInc = 0.01, reltol = 0.01)$Pu

  sim_dir_null = Xmcupo.sevsample(group.data)$`p value`

  ####

  res = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                          otu_group_1, otu_group_2, X_group_1 = cov1, X_group_2 = cov2,
                          nu = 10 ^ (seq(-1, 4, 0.25)), sigma = sqrt(10), verbose = FALSE)

  sim_PJAP_null_cov = res$PJAP
  sim_PJAPind_null_cov = res$PJAP_ind
  beta_null_cov = res$beta

  for(k in 1:length(inc)){

    otu_group_2[label[37], ] = round((1 + inc[k]) * otu_group_2[label[37], ])

    res = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                            otu_group_1, otu_group_2, X_group_1 = "default", X_group_2 = "default",
                            nu = 10 ^ (seq(-1, 4, 0.25)), sigma = sqrt(10))

    sim_PJAP_alt[k] = res$PJAP
    sim_PJAPind_alt[k] = res$PJAP_ind
    beta_alt[k] = res$beta

    group.data = list(x1 = t(otu_group_1), x2 = t(otu_group_2))
    nt = nodetest(pstrct, group.data)
    sim_w_alt[k] = phyloscan(pstrct, nt$w, nthread = 2, gridInc = 0.01, reltol = 0.01)$Pu
    sim_dir_alt[k] = Xmcupo.sevsample(group.data)$`p value`

    ####

    res = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                            otu_group_1, otu_group_2, X_group_1 = cov1, X_group_2 = cov2,
                            nu = 10 ^ (seq(-1, 4, 0.25)), sigma = sqrt(10))

    sim_PJAP_alt_cov[k] = res$PJAP
    sim_PJAPind_alt_cov[k] = res$PJAP_ind
    beta_alt_cov[k] = res$beta

  }

  res = c(sim_PJAP_null = sim_PJAP_null,
          sim_PJAP_alt = sim_PJAP_alt,
          sim_PJAPind_null = sim_PJAPind_null,
          sim_PJAPind_alt = sim_PJAPind_alt,
          sim_w_null = sim_w_null, sim_w_alt = sim_w_alt,
          sim_w_null_ind = sim_w_null_ind, sim_w_alt_ind = sim_w_alt_ind,
          sim_dir_null = sim_dir_null, sim_dir_alt = sim_dir_alt,
          beta_null = beta_null, beta_alt = beta_alt,
          sim_PJAP_null_cov = sim_PJAP_null_cov,
          sim_PJAP_alt_cov = sim_PJAP_alt_cov,
          beta_null_cov = beta_null_cov, beta_alt_cov = beta_alt_cov)

  if(save == TRUE){
    file = paste("./Result_sim4/sim_", seed, ".RData", sep = "")
    save(res, file = file)
  }else{
    return(res)
  }
}

##################################################################################################
##################################################################################################

seed = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
AG_simulate_senario_cov(seed, inc = 2, inc_cov = 2, save = FALSE)

