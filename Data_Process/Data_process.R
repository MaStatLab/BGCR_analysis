
rm(list = ls())

library(data.table)
library(phyloseq)
library(ape)

AG_tree = read_tree_greengenes("97_otus.tree")        #-import the phylogenetric tree
otu = fread('ag_fecal_from_biom.txt', header = TRUE)  #-import the otu table


otu_matrix = do.call(rbind, otu)
otu_matrix = apply(otu_matrix,1, as.numeric)
row.names(otu_matrix) = otu_matrix[, "OTUID"]
otu_neat = otu_matrix[, -1]
rm(otu_matrix)


#################################################################################################
####----choose subset of otus based on some criterion

otu_rowsum = apply(otu_neat,1,sum)
num_otu_top = 100  ####----number of top otus to be used in the analysis
cutoff = sort(otu_rowsum, decreasing = TRUE)[num_otu_top]
otu_top = otu_neat[otu_rowsum>=cutoff, ]

####---create a trimmed otu tree with only the otus selected above

diff = setdiff(AG_tree$tip.label, row.names(otu_top))
tree = drop.tip(AG_tree, tip=diff, trim.internal = TRUE, subtree = FALSE, root.edge = 0)

#################################################################################################
#################################################################################################
#################################################################################################
###---load the AG data

ag_fecal = fread("ag_fecal.txt")

#################################################################################################
####----remove the ununsed variables

rm(otu, otu_neat, AG_tree, cutoff, diff, num_otu_top, otu_rowsum)

#################################################################################################
####----save the data

save(otu_top, tree, ag_fecal, file="AG.RData")
