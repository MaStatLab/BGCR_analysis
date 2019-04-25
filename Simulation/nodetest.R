

nodetest <- function(pstrct, group.data) {
  phylotree <- pstrct$phylotree
  phyloparent <- pstrct$phyloparent
  phylochildren <- pstrct$phylochildren
  
  MoMp <- matrix(0, phylotree$Nnode, 2, dimnames = list(c(), c("Node", "p-value")))
  for (gd in group.data) if (ncol(gd) != ntaxa(pstrct$phylotree)) 
    stop("Number of columns in group.data is not equal to number of leaves")
  for (i in 1:phylotree$Nnode) {
    NodeID <- i + ntaxa(phylotree)
    c1 <- phylochildren[NodeID, 1]
    c2 <- phylochildren[NodeID, 2]
    gp <- lapply(group.data, function(g) {
      k1 <- g[, pstrct$descendant[c1, ]]
      k2 <- g[, pstrct$descendant[c2, ]]
      if (is.matrix(k1)) 
        k1 <- rowSums(k1)
      if (is.matrix(k2)) 
        k2 <- rowSums(k2)
      tm <- cbind(k1, k2)
      tm[rowSums(tm) != 0, ]
    })
    MoMp[i, ] <- c(NodeID, Xmcupo.sevsample(gp)$`p value`)
  }
  
  ts <- matrix(0, phylotree$Nnode, 4, dimnames = list(
    c(), c("Node1", "Node2", "Node3", "stat")))
  for (i in 2:phylotree$Nnode) {
    NodeID <- i + ntaxa(phylotree)
    Par <- phyloparent[NodeID]
    GPar <- phyloparent[Par]
    if (GPar == 0) 
      next
    ts[i, 1:3] <- c(GPar, Par, NodeID)
    ts[i, 4] <- sum(qchisq(1 - MoMp[c(NodeID, Par, GPar) - ntaxa(phylotree), 2], 1))
  }
  ts <- ts[rowSums(ts[, 1:3]) != 0, ]
  list(MoMp = MoMp, tripletstat = ts, w = max(ts[, 4]))
}