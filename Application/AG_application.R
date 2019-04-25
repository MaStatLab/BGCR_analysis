
rm(list = ls())

library(ape)
source("BGCR.R")
load("AG.RData")

ag_fecal = ag_fecal[ag_fecal$SUBSET_AGE == "true" & ag_fecal$SUBSET_BMI == "true", ]
otu_top = otu_top[, intersect(colnames(otu_top), ag_fecal$SampleID) ]

##################################################################################################

n = dim(ag_fecal)[1]
COV = NULL

DIA = ag_fecal$SUBSET_DIABETES
dia = rep(NA, n)

for(i in 1:n){
  if(DIA[i] == "false"){
    dia[i] = 1
  }else if(DIA[i] == "true"){
    dia[i] = 0
  }
}


IBD = ag_fecal$SUBSET_IBD
ibd = rep(NA, n)

for(i in 1:n){
  if(IBD[i] == "false"){
    ibd[i] = 1
  }else if(IBD[i] == "true"){
    ibd[i] = 0
  }
}



ANTI = ag_fecal$SUBSET_ANTIBIOTIC_HISTORY
anti = rep(NA, n)

for(i in 1:n){
  if(ANTI[i] == "false"){
    anti[i] = 1
  }else if(ANTI[i] == "true"){
    anti[i] = 0
  }
}

####

SEX = ag_fecal$SEX
SEX = factor(SEX)
lev = levels(SEX)

sex = rep(NA, n)

for(i in 1:n){
  if(SEX[i] == lev[1]){
    sex[i] = 1
  }else if(SEX[i] == lev[2]){
    sex[i] = 0
  }
}

####

ALC = ag_fecal$ALCOHOL_FREQUENCY
ALC = factor(ALC)
lev = levels(ALC)

alc = rep(NA, n)

for(i in 1:n){
  if(ALC[i] == lev[1] | ALC[i] == lev[5]){
    alc[i] = 1
  }else if(ALC[i] == lev[2] | ALC[i] == lev[3] | ALC[i] == lev[4]){
    alc[i] = 0
  }
}


####

VEGE = ag_fecal$VEGETABLE_FREQUENCY
VEGE = factor(VEGE)
lev = levels(VEGE)
n = length(VEGE)

veg = rep(NA, n)

for(i in 1:n){
  if(VEGE[i] == lev[1] | VEGE[i] == lev[5]){
    veg[i] = 1
  }else if(VEGE[i] == lev[2] | VEGE[i] == lev[3] | VEGE[i] == lev[4]){
    veg[i] = 0
  }
}

####

FRUIT = ag_fecal$FRUIT_FREQUENCY
FRUIT = factor(FRUIT)
lev = levels(FRUIT)

fruit = rep(NA, n)

for(i in 1:n){
  if(FRUIT[i] == lev[1] | FRUIT[i] == lev[5]){
    fruit[i] = 1
  }else if(FRUIT[i] == lev[2] | FRUIT[i] == lev[3] | FRUIT[i] == lev[4]){
    fruit[i] = 0
  }
}

####

SEA = ag_fecal$SEAFOOD_FREQUENCY
SEA = factor(SEA)
lev = levels(SEA)

sea = rep(NA, n)

for(i in 1:n){
  if(SEA[i] == lev[1] | SEA[i] == lev[5]){
    sea[i] = 1
  }else if(SEA[i] == lev[2] | SEA[i] == lev[3] | SEA[i] == lev[4]){
    sea[i] = 0
  }
}

####

SUGAR = ag_fecal$SUGARY_SWEETS_FREQUENCY
SUGAR = factor(SUGAR)
lev = levels(SUGAR)

sugar = rep(NA, n)

for(i in 1:n){
  if(SUGAR[i] == lev[1] | SUGAR[i] == lev[5]){
    sugar[i] = 1
  }else if(SUGAR[i] == lev[2] | SUGAR[i] == lev[3] | SUGAR[i] == lev[4]){
    sugar[i] = 0
  }
}


####

MEAT = ag_fecal$MEAT_EGGS_FREQUENCY
MEAT = factor(MEAT)
lev = levels(MEAT)

meat = rep(NA, n)

for(i in 1:n){
  if(MEAT[i] == lev[1] | MEAT[i] == lev[5]){
    meat[i] = 1
  }else if(MEAT[i] == lev[2] | MEAT[i] == lev[3] | MEAT[i] == lev[4]){
    meat[i] = 0
  }
}

####

PRO = ag_fecal$PROBIOTIC_FREQUENCY
PRO = factor(PRO)
lev = levels(PRO)

pro = rep(NA, n)

for(i in 1:n){
  if(PRO[i] == lev[1] | PRO[i] == lev[5]){
    pro[i] = 1
  }else if(PRO[i] == lev[2] | PRO[i] == lev[3] | PRO[i] == lev[4]){
    pro[i] = 0
  }
}


####

GRAIN = ag_fecal$WHOLE_GRAIN_FREQUENCY
GRAIN = factor(GRAIN)
lev = levels(GRAIN)

grain = rep(NA, n)

for(i in 1:n){
  if(GRAIN[i] == lev[1] | GRAIN[i] == lev[5]){
    grain[i] = 1
  }else if(GRAIN[i] == lev[2] | GRAIN[i] == lev[3] | GRAIN[i] == lev[4]){
    grain[i] = 0
  }
}

####

MILK = ag_fecal$MILK_CHEESE_FREQUENCY
MILK = factor(MILK)
lev = levels(MILK)

milk = rep(NA, n)

for(i in 1:n){
  if(MILK[i] == lev[1] | MILK[i] == lev[5]){
    milk[i] = 1
  }else if(MILK[i] == lev[2] | MILK[i] == lev[3] | MILK[i] == lev[4]){
    milk[i] = 0
  }
}


##################################################################################################

COV = cbind(sex, dia, ibd, anti, pro, alc, veg, fruit, sea, meat, sugar, grain, milk)
complete_case = complete.cases(COV)
complete_id = ag_fecal$SampleID[complete_case]
COV = COV[complete.cases(COV), ]

##################################################################################################
####----FRUIT

false_group_id = ag_fecal$SampleID[which(fruit == 0)]
true_group_id = ag_fecal$SampleID[which(fruit == 1)]

otu_group_1 = otu_top[, intersect(complete_id, false_group_id)]
otu_group_2 = otu_top[, intersect(complete_id, true_group_id)]

res_fruit_none = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                              otu_group_1, otu_group_2, X_group_1 = "default", X_group_2 = "default",
                              nu = 10 ^ (seq(-1, 4)), sigma = sqrt(10), verbose = TRUE)

cov1 = COV[which(fruit[complete_case] == 0), c(1, 2, 3, 4, 5)]
cov2 = COV[which(fruit[complete_case] == 1), c(1, 2, 3, 4, 5)]

res_fruit = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                                  otu_group_1, otu_group_2, X_group_1 = cov1, X_group_2 = cov2,
                                  nu = 10 ^ (seq(-1, 4)), sigma = sqrt(10), verbose = TRUE)

cov1 = COV[which(fruit[complete_case] == 0), -8]
cov2 = COV[which(fruit[complete_case] == 1), -8]

res_fruit_cov = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                              otu_group_1, otu_group_2, X_group_1 = cov1, X_group_2 = cov2,
                              nu = 10 ^ (seq(-1, 4)), sigma = sqrt(10), verbose = TRUE)


list_fruit = list(res_fruit_none = res_fruit_none, res_fruit = res_fruit, res_fruit_cov = res_fruit_cov)
save(list_fruit, file = "fruit.RData")

####---SEA


false_group_id = ag_fecal$SampleID[which(sea == 0)]
true_group_id = ag_fecal$SampleID[which(sea == 1)]

otu_group_1 = otu_top[, intersect(complete_id, false_group_id)]
otu_group_2 = otu_top[, intersect(complete_id, true_group_id)]


res_sea_none = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                            otu_group_1, otu_group_2, X_group_1 = "default", X_group_2 = "default",
                            nu = 10 ^ (seq(-1, 4)), sigma = sqrt(10), verbose = TRUE)

cov1 = COV[which(sea[complete_case] == 0), c(1, 2, 3, 4, 5)]
cov2 = COV[which(sea[complete_case] == 1), c(1, 2, 3, 4, 5)]

res_sea = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                              otu_group_1, otu_group_2, X_group_1 = cov1, X_group_2 = cov2,
                              nu = 10 ^ (seq(-1, 4)), sigma = sqrt(10), verbose = TRUE)

cov1 = COV[which(sea[complete_case] == 0), -9]
cov2 = COV[which(sea[complete_case] == 1), -9]

res_sea_cov = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                                  otu_group_1, otu_group_2, X_group_1 = cov1, X_group_2 = cov2,
                                  nu = 10 ^ (seq(-1, 4)), sigma = sqrt(10), verbose = TRUE)

list_sea = list(res_sea_none = res_sea_none, res_sea = res_sea, res_sea_cov = res_sea_cov)
save(list_sea, file = "sea.RData")


####---VEGETABLE


false_group_id = ag_fecal$SampleID[which(veg == 0)]
true_group_id = ag_fecal$SampleID[which(veg == 1)]

otu_group_1 = otu_top[, intersect(complete_id, false_group_id)]
otu_group_2 = otu_top[, intersect(complete_id, true_group_id)]


res_veg_none = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                            otu_group_1, otu_group_2, X_group_1 = "default", X_group_2 = "default",
                            nu = 10 ^ (seq(-1, 4)), sigma = sqrt(10), verbose = TRUE)

cov1 = COV[which(veg[complete_case] == 0), c(1, 2, 3, 4, 5)]
cov2 = COV[which(veg[complete_case] == 1), c(1, 2, 3, 4, 5)]

res_veg = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                            otu_group_1, otu_group_2, X_group_1 = cov1, X_group_2 = cov2,
                            nu = 10 ^ (seq(-1, 4)), sigma = sqrt(10), verbose = TRUE)

cov1 = COV[which(veg[complete_case] == 0), -7]
cov2 = COV[which(veg[complete_case] == 1), -7]

res_veg_cov = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                                otu_group_1, otu_group_2, X_group_1 = cov1, X_group_2 = cov2,
                                nu = 10 ^ (seq(-1, 4)), sigma = sqrt(10), verbose = TRUE)


list_veg = list(res_veg_none = res_veg_none, res_veg = res_veg, res_veg_cov = res_veg_cov)
save(list_veg, file = "veg.RData")

####---ALCHOHOL


false_group_id = ag_fecal$SampleID[which(alc == 0)]
true_group_id = ag_fecal$SampleID[which(alc == 1)]

otu_group_1 = otu_top[, intersect(complete_id, false_group_id)]
otu_group_2 = otu_top[, intersect(complete_id, true_group_id)]


res_alc_none = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                                 otu_group_1, otu_group_2, X_group_1 = "default", X_group_2 = "default",
                                 nu = 10 ^ (seq(-1, 4)), sigma = sqrt(10), verbose = TRUE)

cov1 = COV[which(alc[complete_case] == 0), c(1, 2, 3, 4, 5)]
cov2 = COV[which(alc[complete_case] == 1), c(1, 2, 3, 4, 5)]

res_alc = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                            otu_group_1, otu_group_2, X_group_1 = cov1, X_group_2 = cov2,
                            nu = 10 ^ (seq(-1, 4)), sigma = sqrt(10), verbose = TRUE)

cov1 = COV[which(alc[complete_case] == 0), -6]
cov2 = COV[which(alc[complete_case] == 1), -6]

res_alc_cov = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                                otu_group_1, otu_group_2, X_group_1 = cov1, X_group_2 = cov2,
                                nu = 10 ^ (seq(-1, 4)), sigma = sqrt(10), verbose = TRUE)

list_alc = list(res_alc_none = res_alc_none, res_alc = res_alc, res_alc_cov = res_alc_cov)
save(list_alc, file = "alc.RData")

####---MEAT


false_group_id = ag_fecal$SampleID[which(meat == 0)]
true_group_id = ag_fecal$SampleID[which(meat == 1)]

otu_group_1 = otu_top[, intersect(complete_id, false_group_id)]
otu_group_2 = otu_top[, intersect(complete_id, true_group_id)]


res_meat_none = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                                 otu_group_1, otu_group_2, X_group_1 = "default", X_group_2 = "default",
                                 nu = 10 ^ (seq(-1, 4)), sigma = sqrt(10), verbose = TRUE)

cov1 = COV[which(meat[complete_case] == 0), c(1, 2, 3, 4, 5)]
cov2 = COV[which(meat[complete_case] == 1), c(1, 2, 3, 4, 5)]

res_meat = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                            otu_group_1, otu_group_2, X_group_1 = cov1, X_group_2 = cov2,
                            nu = 10 ^ (seq(-1, 4)), sigma = sqrt(10), verbose = TRUE)

cov1 = COV[which(meat[complete_case] == 0), -10]
cov2 = COV[which(meat[complete_case] == 1), -10]

res_meat_cov = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                                otu_group_1, otu_group_2, X_group_1 = cov1, X_group_2 = cov2,
                                nu = 10 ^ (seq(-1, 4)), sigma = sqrt(10), verbose = TRUE)


list_meat = list(res_meat_none = res_meat_none, res_meat = res_meat, res_meat_cov = res_meat_cov)
save(list_meat, file = "meat.RData")

####---SUGAR


false_group_id = ag_fecal$SampleID[which(sugar == 0)]
true_group_id = ag_fecal$SampleID[which(sugar == 1)]

otu_group_1 = otu_top[, intersect(complete_id, false_group_id)]
otu_group_2 = otu_top[, intersect(complete_id, true_group_id)]


res_sugar_none = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                                  otu_group_1, otu_group_2, X_group_1 = "default", X_group_2 = "default",
                                  nu = 10 ^ (seq(-1, 4)), sigma = sqrt(10), verbose = TRUE)

cov1 = COV[which(sugar[complete_case] == 0), c(1, 2, 3, 4, 5)]
cov2 = COV[which(sugar[complete_case] == 1), c(1, 2, 3, 4, 5)]

res_sugar = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                             otu_group_1, otu_group_2, X_group_1 = cov1, X_group_2 = cov2,
                             nu = 10 ^ (seq(-1, 4)), sigma = sqrt(10), verbose = TRUE)

cov1 = COV[which(sugar[complete_case] == 0), -11]
cov2 = COV[which(sugar[complete_case] == 1), -11]

res_sugar_cov = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                                 otu_group_1, otu_group_2, X_group_1 = cov1, X_group_2 = cov2,
                                 nu = 10 ^ (seq(-1, 4)), sigma = sqrt(10), verbose = TRUE)

list_sugar = list(res_sugar_none = res_sugar_none, res_sugar = res_sugar, res_sugar_cov = res_sugar_cov)
save(list_sugar, file = "sugar.RData")


####---GRAIN


false_group_id = ag_fecal$SampleID[which(grain == 0)]
true_group_id = ag_fecal$SampleID[which(grain == 1)]

otu_group_1 = otu_top[, intersect(complete_id, false_group_id)]
otu_group_2 = otu_top[, intersect(complete_id, true_group_id)]


res_grain_none = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                                   otu_group_1, otu_group_2, X_group_1 = "default", X_group_2 = "default",
                                   nu = 10 ^ (seq(-1, 4)), sigma = sqrt(10), verbose = TRUE)

cov1 = COV[which(grain[complete_case] == 0), c(1, 2, 3, 4, 5)]
cov2 = COV[which(grain[complete_case] == 1), c(1, 2, 3, 4, 5)]

res_grain = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                              otu_group_1, otu_group_2, X_group_1 = cov1, X_group_2 = cov2,
                              nu = 10 ^ (seq(-1, 4)), sigma = sqrt(10), verbose = TRUE)

cov1 = COV[which(grain[complete_case] == 0), -12]
cov2 = COV[which(grain[complete_case] == 1), -12]

res_grain_cov = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                                  otu_group_1, otu_group_2, X_group_1 = cov1, X_group_2 = cov2,
                                  nu = 10 ^ (seq(-1, 4)), sigma = sqrt(10), verbose = TRUE)



list_grain = list(res_grain_none = res_grain_none, res_grain = res_grain, res_grain_cov = res_grain_cov)
save(list_grain, file = "grain.RData")

####---MILK


false_group_id = ag_fecal$SampleID[which(milk == 0)]
true_group_id = ag_fecal$SampleID[which(milk == 1)]

otu_group_1 = otu_top[, intersect(complete_id, false_group_id)]
otu_group_2 = otu_top[, intersect(complete_id, true_group_id)]


res_milk_none = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                                   otu_group_1, otu_group_2, X_group_1 = "default", X_group_2 = "default",
                                   nu = 10 ^ (seq(-1, 4)), sigma = sqrt(10), verbose = TRUE)

cov1 = COV[which(milk[complete_case] == 0), c(1, 2, 3, 4, 5)]
cov2 = COV[which(milk[complete_case] == 1), c(1, 2, 3, 4, 5)]

res_milk = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                              otu_group_1, otu_group_2, X_group_1 = cov1, X_group_2 = cov2,
                              nu = 10 ^ (seq(-1, 4)), sigma = sqrt(10), verbose = TRUE)

cov1 = COV[which(milk[complete_case] == 0), -13]
cov2 = COV[which(milk[complete_case] == 1), -13]

res_milk_cov = BGCR(PrJAP = 0.5, sum_PrMAP = "default", threshold = 0.005, gamma = 0, tree,
                                  otu_group_1, otu_group_2, X_group_1 = cov1, X_group_2 = cov2,
                                  nu = 10 ^ (seq(-1, 4)), sigma = sqrt(10), verbose = TRUE)


list_milk = list(res_milk_none = res_milk_none, res_milk = res_milk, res_milk_cov = res_milk_cov)
save(list_milk, file = "milk.RData")



##################################################################################################
####----unbalance

ub1 = table(COV[, "fruit"], COV[, "veg"])
rownames(ub1) = c("False", "True")
colnames(ub1) = c("False", "True")

ub2 = table(COV[, "fruit"], COV[, "grain"])
rownames(ub2) = c("False", "True")
colnames(ub2) = c("False", "True")


setEPS()
postscript("app_unbalance.eps", width = 8, height = 4)
par(mfrow = c(1, 2), mar = c(2,2,2,2))
mosaicplot(ub1, xlab = "Fruit", ylab = "Vegetable", main = NULL)
mosaicplot(ub2, xlab = "Fruit", ylab = "Grain", main = NULL)
dev.off()


##################################################################################################
####----PMAP plot

setEPS()
postscript("app_fruit.eps", width = 6, height = 5)
par(mar = c(1,1,1,1))
plot(list_fruit$res_fruit_cov, main = "Fruit", cex = 0.25)
dev.off()

setEPS()
postscript("app_sea.eps", width = 6, height = 5)
plot(list_sea$res_sea_cov, main = "Seafood", cex = 0.25)
dev.off()

setEPS()
postscript("app_veg.eps", width = 6, height = 5)
plot(list_veg$res_veg_cov, main = "Vegetable", cex = 0.25)
dev.off()

setEPS()
postscript("app_grain.eps", width = 6, height = 5)
plot(list_grain$res_grain_cov, main = "Grain", cex = 0.25)
dev.off()


setEPS()
postscript("app_fruit_none.eps", width = 6, height = 5)
par(mar = c(1,1,1,1))
plot(list_fruit$res_fruit_none, main = "Fruit", cex = 0.25)
dev.off()

setEPS()
postscript("app_sea_none.eps", width = 6, height = 5)
plot(list_sea$res_sea_none, main = "Seafood", cex = 0.25)
dev.off()

setEPS()
postscript("app_veg_none.eps", width = 6, height = 5)
plot(list_veg$res_veg_none, main = "Vegetable", cex = 0.25)
dev.off()

setEPS()
postscript("app_grain_none.eps", width = 6, height = 5)
plot(list_grain$res_grain_none, main = "Grain", cex = 0.25)
dev.off()


setEPS()
postscript("app_fruit_half.eps", width = 6, height = 5)
par(mar = c(1,1,1,1))
plot(list_fruit$res_fruit, main = "Fruit", cex = 0.25)
dev.off()

setEPS()
postscript("app_sea_half.eps", width = 6, height = 5)
plot(list_sea$res_sea, main = "Seafood", cex = 0.25)
dev.off()

setEPS()
postscript("app_veg_half.eps", width = 6, height = 5)
plot(list_veg$res_veg, main = "Vegetable", cex = 0.25)
dev.off()

setEPS()
postscript("app_grain_half.eps", width = 6, height = 5)
plot(list_grain$res_grain, main = "Grain", cex = 0.25)
dev.off()

