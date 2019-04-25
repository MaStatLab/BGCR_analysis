
phylostructure <- function(phylotree) {
  if(class(phylotree) != "phylo")
    stop("phylotree is not a phylo class")
  phyloparent <- numeric(phylotree$Nnode + ntaxa(phylotree))
  phylochildren <- matrix(0, phylotree$Nnode + ntaxa(phylotree), 2)
  for (i in 1:nrow(phylotree$edge)) {
    i1 <- phylotree$edge[i, 1]
    i2 <- phylotree$edge[i, 2]
    if (i1 <= ntaxa(phylotree)) 
      stop("Internal node label is not larger than ntaxa(phylotree)")
    if (i2 > ntaxa(phylotree) && i1 > i2) 
      stop("Parent node label is larger than child internal node label")
    if (all(phylochildren[i1, ] > 0)) 
      stop("Phylogenetic tree is not binary")
    
    phyloparent[i2] <- i1
    if (phylochildren[i1, 1] == 0) 
      phylochildren[i1, 1] <- i2 else phylochildren[i1, 2] <- i2
  }
  
  descendant <- matrix(FALSE, phylotree$Nnode + ntaxa(phylotree), ntaxa(phylotree))
  for (i in 1:ntaxa(phylotree)) descendant[i, i] <- TRUE
  processed <- logical(phylotree$Nnode + ntaxa(phylotree))
  processed[1:ntaxa(phylotree)] <- TRUE
  while (!all(processed)) {
    for (i in (ntaxa(phylotree) + 1):(ntaxa(phylotree) + phylotree$Nnode)) {
      if (all(processed[phylochildren[i, ]])) {
        descendant[i, descendant[phylochildren[i, 1], ]] <- TRUE
        descendant[i, descendant[phylochildren[i, 2], ]] <- TRUE
        processed[i] <- TRUE
      }
    }
  }
  list(phylotree = phylotree, phylochildren = phylochildren, phyloparent = phyloparent, 
       descendant = descendant)
}