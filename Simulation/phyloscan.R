initLabel <- function(pstrct) {
  phylotree <- pstrct$phylotree
  phylochildren <- pstrct$phylochildren
  label <- numeric(phylotree$Nnode)
  Clabel <- 0
  for (i in 1:phylotree$Nnode) {
    if (label[i] != 0) 
      next
    cl <- phylochildren[i + ntaxa(phylotree), 1]
    cr <- phylochildren[i + ntaxa(phylotree), 2]
    flag <- FALSE
    if (cl > ntaxa(phylotree)) {
      if (phylochildren[cl, 1] > ntaxa(phylotree)) {
        Clabel <- Clabel + 1
        label[c(i, cl - ntaxa(phylotree), phylochildren[cl, 1] - ntaxa(phylotree))] <- Clabel
        next
      } else if (phylochildren[cl, 2] > ntaxa(phylotree)) {
        Clabel <- Clabel + 1
        label[c(i, cl - ntaxa(phylotree), phylochildren[cl, 2] - ntaxa(phylotree))] <- Clabel
        next
      }
    }
    if (cr > ntaxa(phylotree)) {
      if (phylochildren[cr, 1] > ntaxa(phylotree)) {
        Clabel <- Clabel + 1
        label[c(i, cr - ntaxa(phylotree), phylochildren[cr, 1] - ntaxa(phylotree))] <- Clabel
        next
      } else if (phylochildren[cr, 2] > ntaxa(phylotree)) {
        Clabel <- Clabel + 1
        label[c(i, cr - ntaxa(phylotree), phylochildren[cr, 2] - ntaxa(phylotree))] <- Clabel
        next
      }
    }
  }
  for (i in 1:phylotree$Nnode) {
    if (label[i] != 0) 
      next
    cl <- phylochildren[i + ntaxa(phylotree), 1]
    cr <- phylochildren[i + ntaxa(phylotree), 2]
    flag <- FALSE
    if (cl > ntaxa(phylotree)) {
      Clabel <- Clabel + 1
      label[c(i, cl - ntaxa(phylotree))] <- Clabel
      next
    }
    if (cr > ntaxa(phylotree)) {
      Clabel <- Clabel + 1
      label[c(i, cr - ntaxa(phylotree))] <- Clabel
      next
    }
    Clabel <- Clabel + 1
    label[i] <- Clabel
  }
  label
}

initGroup3 <- function(label, pstrct) {
  phylotree <- pstrct$phylotree
  phyloparent <- pstrct$phyloparent
  phylochildren <- pstrct$phylochildren
  Group3 <- matrix(0, phylotree$Nnode - 3, 3)
  Cntr <- 0
  for (i in 3:phylotree$Nnode) {
    if (i == phylochildren[1 + ntaxa(phylotree), 2] - ntaxa(phylotree)) 
      next
    
    p <- phyloparent[i + ntaxa(phylotree)]
    G3Seq <- c(phyloparent[p] - ntaxa(phylotree), p - ntaxa(phylotree), i)
    if (label[G3Seq[1]] == label[G3Seq[2]] && label[G3Seq[2]] == label[G3Seq[3]]) 
      next
    Cntr <- Cntr + 1
    Group3[Cntr, ] <- c(phyloparent[p] - ntaxa(phylotree), p - ntaxa(phylotree), i)
  }
  Group3[1:Cntr, ]
}

initInteraction2 <- function(G3, label) {
  Interaction2 <- list()
  for (i in 1:nrow(G3)) {
    if (i > 1) 
      for (j in 1:(i - 1)) {
        if (length(intersect(G3[i, ], G3[j, ])) == 2) {
          if (length(Interaction2) < i) 
            Interaction2[[i]] <- G3[j, ] 
          else 
            Interaction2[[i]] <- c(Interaction2[[i]], G3[j, ])
        }
      }
    if (length(Interaction2) < i) {
      Interaction2[[i]] <- -1
      lseq1 <- label[G3[i, ]]
      if (length(unique(lseq1)) == 2) {
        attr(Interaction2[[i]], "type") <- "112"
        dup <- lseq1[duplicated(lseq1)]
        attr(Interaction2[[i]], "idxdist") <- c(dup, setdiff(lseq1, dup))
      } else if (length(unique(lseq1)) == 3) {
        attr(Interaction2[[i]], "type") <- "123"
        attr(Interaction2[[i]], "idxdist") <- lseq1
      }
      next
    }
    seq1 <- G3[i, ]
    seq2 <- Interaction2[[i]]
    allseq <- c(seq1, seq2)
    lseq1 <- label[seq1]
    lseq2 <- label[seq2]
    lallseq <- label[allseq]
    if (length(seq2) == 3) {
      if (length(unique(lallseq)) == 4) {
        attr(Interaction2[[i]], "type") <- "123124"
        dup <- intersect(lseq1, lseq2)
        DistSeq <- c(dup, setdiff(lseq1, dup), setdiff(lseq2, dup))
        attr(Interaction2[[i]], "idxdist") <- DistSeq
      } else if (length(unique(lallseq)) == 2) {
        attr(Interaction2[[i]], "type") <- "112122"
        dup <- (lseq1)[duplicated(lseq1)]
        DistSeq <- c(dup, setdiff(lseq1, dup))
        attr(Interaction2[[i]], "idxdist") <- DistSeq
      } else if (length(unique(lallseq)) == 3) {
        if (length(unique(lseq1)) == 2 && length(unique(lseq2)) == 2) {
          attr(Interaction2[[i]], "type") <- "112113"
          dup <- (lseq1)[duplicated(lseq1)]
          DistSeq <- c(dup, setdiff(lseq1, dup), setdiff(lseq2, dup))
          attr(Interaction2[[i]], "idxdist") <- DistSeq
        }
        if (length(unique(lseq1)) == 2 && length(unique(lseq2)) == 3) {
          attr(Interaction2[[i]], "type") <- "112123"
          dup <- (lseq1)[duplicated(lseq1)]
          DistSeq <- c(dup, setdiff(lseq1, dup), setdiff(lseq2, lseq1))
          attr(Interaction2[[i]], "idxdist") <- DistSeq
        }
        if (length(unique(lseq1)) == 3 && length(unique(lseq2)) == 2) {
          attr(Interaction2[[i]], "type") <- "123112"
          dup <- (lseq2)[duplicated(lseq2)]
          DistSeq <- c(dup, setdiff(lseq2, dup), setdiff(lseq1, lseq2))
          attr(Interaction2[[i]], "idxdist") <- DistSeq
        }
        
      }
    } else if (length(seq2) == 6) {
      if (length(unique(lseq1)) == 2 && length(unique(lseq2)) == 3) {
        attr(Interaction2[[i]], "type") <- "112122123"
        dup <- (lseq1)[duplicated(lseq1)]
        DistSeq <- c(dup, setdiff(lseq1, dup), setdiff(lseq2, lseq1))
        attr(Interaction2[[i]], "idxdist") <- DistSeq
      } else if (length(unique(lseq1)) == 2 && length(unique(lseq2)) == 4) {
        attr(Interaction2[[i]], "type") <- "112123124"
        dup <- (lseq1)[duplicated(lseq1)]
        DistSeq <- c(dup, setdiff(lseq1, dup), setdiff(lseq2, lseq1))
        attr(Interaction2[[i]], "idxdist") <- DistSeq
      } else if (length(unique(lseq1)) == 3 && length(unique(lseq2)) == 2) {
        attr(Interaction2[[i]], "type") <- "123112122"
        dup <- unique(lseq2)
        DistSeq <- c(dup, setdiff(lseq1, dup))
        attr(Interaction2[[i]], "idxdist") <- DistSeq
      } else if (length(unique(lseq1)) == 3 && length(unique(lseq2)) == 3) {
        attr(Interaction2[[i]], "type") <- "123122124"
        tt <- table(lallseq)
        d1 <- as.numeric(names(tt)[which(tt == 3)])
        d2 <- as.numeric(names(tt)[which(tt == 4)])
        DistSeq <- c(d1, d2, setdiff(lseq1, c(d1, d2)), setdiff(lseq2, lseq1))
        attr(Interaction2[[i]], "idxdist") <- DistSeq
      }
    }
  }
  Interaction2
}

initInteraction1 <- function(G3, label) {
  Interaction1 <- list()
  for (i in 1:nrow(G3)) {
    iC <- 0
    Interaction1[[i]] <- list()
    for (j in 1:(i - 1)) {
      if (length(intersect(G3[i, ], G3[j, ])) == 1) {
        iC <- iC + 1
        Interaction1[[i]][[iC]] <- G3[j, ]
      }
    }
    if (length(Interaction1[[i]]) == 0) 
      next
    for (j in 1:length(Interaction1[[i]])) {
      seq1 <- G3[i, ]
      seq2 <- Interaction1[[c(i, j)]]
      lseq1 <- label[seq1]
      lseq2 <- label[seq2]
      if (length(unique(lseq1)) == 2) {
        dup <- (lseq1)[duplicated(lseq1)]
        ndup <- setdiff(lseq1, dup)
        if (sum(lseq2 == dup) == 2) {
          attr(Interaction1[[c(i, j)]], "type") <- "e112113"
          DistSeq <- c(dup, ndup, setdiff(lseq2, dup))
          attr(Interaction1[[c(i, j)]], "idxdist") <- DistSeq
        } else if (sum(lseq2 == dup) == 1 && length(unique(lseq2)) == 2) {
          attr(Interaction1[[c(i, j)]], "type") <- "e112133"
          DistSeq <- c(dup, ndup, setdiff(lseq2, dup))
          attr(Interaction1[[c(i, j)]], "idxdist") <- DistSeq
        } else if (sum(lseq2 == dup) == 1 && length(unique(lseq2)) == 3) {
          attr(Interaction1[[c(i, j)]], "type") <- "e112134"
          DistSeq <- c(dup, ndup, setdiff(lseq2, dup))
          attr(Interaction1[[c(i, j)]], "idxdist") <- DistSeq
        } else if (sum(lseq2 == ndup) == 2) {
          attr(Interaction1[[c(i, j)]], "type") <- "e112223"
          DistSeq <- c(dup, ndup, setdiff(lseq2, ndup))
          attr(Interaction1[[c(i, j)]], "idxdist") <- DistSeq
        } else if (sum(lseq2 == ndup) == 1 && length(unique(lseq2)) == 2) {
          attr(Interaction1[[c(i, j)]], "type") <- "e112233"
          DistSeq <- c(dup, ndup, setdiff(lseq2, ndup))
          attr(Interaction1[[c(i, j)]], "idxdist") <- DistSeq
        } else if (sum(lseq2 == ndup) == 1 && length(unique(lseq2)) == 3) {
          attr(Interaction1[[c(i, j)]], "type") <- "e112234"
          DistSeq <- c(dup, ndup, setdiff(lseq2, ndup))
          attr(Interaction1[[c(i, j)]], "idxdist") <- DistSeq
        }
      } else if (length(unique(lseq1)) == 3) {
        dup <- intersect(lseq1, lseq2)
        if (sum(lseq2 == dup) == 2) {
          attr(Interaction1[[c(i, j)]], "type") <- "e123114"
          DistSeq <- c(lseq1, setdiff(lseq2, dup))
          attr(Interaction1[[c(i, j)]], "idxdist") <- DistSeq
        } else if (sum(lseq2 == dup) == 1 && length(setdiff(lseq2, dup)) == 1) {
          attr(Interaction1[[c(i, j)]], "type") <- "e123144"
          DistSeq <- c(lseq1, setdiff(lseq2, dup))
          attr(Interaction1[[c(i, j)]], "idxdist") <- DistSeq
        } else if (sum(lseq2 == dup) == 1 && length(setdiff(lseq2, dup)) == 2) {
          attr(Interaction1[[c(i, j)]], "type") <- "e123145"
          DistSeq <- c(lseq1, setdiff(lseq2, dup))
          attr(Interaction1[[c(i, j)]], "idxdist") <- DistSeq
        }
      }
    }
  }
  Interaction1
}

initInteraction0 <- function(G3, label) {
  Interaction0 <- list()
  for (i in 1:nrow(G3)) {
    iC <- 0
    Interaction0[[i]] <- list()
    if (i == 1) 
      next
    for (j in 1:(i - 1)) {
      if (length(intersect(G3[i, ], G3[j, ])) == 0 && length(
        intersect(label[G3[i, ]], label[G3[j, ]])) > 0) {
        iC <- iC + 1
        Interaction0[[i]][[iC]] <- G3[j, ]
      }
    }
  }
  Interaction0
}

initIndep <- function(G3, label) {
  Indep <- list()
  for (i in 1:nrow(G3)) {
    iC <- 0
    Indep[[i]] <- list()
    if (i == 1) 
      next
    for (j in 1:(i - 1)) {
      if (length(intersect(label[G3[i, ]], label[G3[j, ]])) == 0) {
        iC <- iC + 1
        Indep[[i]][[iC]] <- G3[j, ]
      }
    }
  }
  Indep
}


#' Bounding the scan statistic p-value
#'
#' \code{phyloscan} returns the upper and lower bound of the p-value of the 
#' maximum triplet statistic.
#' 
#' @param pstrct An object returned from function \code{phylostructure}.
#' @param w Observed maximum triplet statistic. This is included in the return 
#' object of function \code{nodetest}.
#' @param nthread Number of parallel thread.
#' @param verbose Should the progress output be printed? Default is TRUE.
#' @param gridInc Increment of sucessive grid points on which the cumulative 
#' distribution function is calculated. Default is \code{1e-3}.
#' @param reltol Relative error of numerical integration. Default is \code{1e-3}.
#' @param btol Relative error of cumulative distribution function on the grid. 
#' Default is \code{1e-6}. Recommend not to change.
#' @details This function calculates the upper and lower bound of the p-value of
#' maximum triplet statistic using the method in the reference paper. 
#' There are two stages in total, both of which are computationally intensive. 
#' Parallel processing is strongly recommended. Total computation time scales up 
#' with \code{w} and size of phylogenetic tree.
#' 
#' @return A list with the following components:
#' \item{\code{Pu}}{Upper bound of p-value}
#' \item{\code{Pl}}{Lower bound of p-value}
#' @author Yunfan Tang
#' @references Tang, Y., Ma, Li. and Nicolae, D. L. (2017). A Phylogenetic Scan Test 
#' on Dirichlet-Tree Multinomial Model for microbiome data. 


#' \href{https://arxiv.org/abs/1610.08974}{arXiv:1610.08974} [stat.AP].
#' @examples
#' library(ape)
#' set.seed(10)
#' 
#' pstrct <- phylostructure(rtree(8))
#' phyloscan(pstrct, 2, nthread = 2, gridInc = 0.02, reltol = 0.02)
#' 
#' ## This example should take about a minute
#' \dontrun{
#' set.seed(10)
#'   
#' pstrct <- phylostructure(rtree(10))
#' p1 <- c(rep(0.12, 3), rep(0.08, 3), rep(0.1, 4))
#' p2 <- p1 + 0.001 * c(c(1,-1), rep(0,8))
#' n <- 1000 #Number of sequences in each sample
#' m <- 200 #Number of samples in each group
#' group.data <- list(x1 = t(rmultinom(m, n, p1)), x2 = t(rmultinom(m, n, p2)))
#' nt <- nodetest(pstrct, group.data) #Generate triplet statistics
#' phyloscan(pstrct, nt$w, nthread = 2, gridInc = 0.01, reltol = 0.01)
#' } 
#' @importFrom R2Cuba cuhre suave
#' @import hashmap
#' @import foreach
#' @import parallel
#' @import doParallel
#' @export

phyloscan <- function(pstrct, w, nthread, verbose = TRUE, gridInc = 0.001, 
                      reltol = 0.001, btol = 1e-06) {
  psn <- 1e-09
  invGP <- round(1/gridInc)
  
  dchisqM <- function(x, df) {
    if (df > 1) 
      return(dchisq(x, df)) else {
        x[x < 0] <- -1
        x[(x >= 0) & (x < psn)] <- psn
        dchisq(x, df)
      }
  }
  
  f31 <- function(x) pchisq(w - x, 2) * dchisqM(x, 1)/pchisqC3
  f21 <- function(x) pchisq(w - x, 1) * dchisqM(x, 1)/pchisqC2
  f11 <- function(x) (x < w) * dchisqM(x, 1)/pchisqC1
  f32 <- function(x) pchisq(w - sum(x), 1) * prod(dchisqM(x, 1))/pchisqC3
  f32s <- function(x) pchisq(w - x, 1) * dchisqM(x, 2)/pchisqC3
  f22 <- function(x) (sum(x) < w) * prod(dchisqM(x, 1))/pchisqC2
  f22s <- function(x) (x < w) * dchisqM(x, 2)/pchisqC2
  f33 <- function(x) (sum(x) < w) * prod(dchisqM(x, 1))/pchisqC3
  f3132s <- function(x) (sum(x) < w) * dchisqM(x[1], 1) * dchisqM(x[2], 2)/pchisqC3
  
  f_1 <- function(x, d) switch(d, f11(x), f21(x), f31(x))
  f_2 <- function(x, d) switch(d - 1, f22(x), f32(x))
  f_2s <- function(x, d) switch(d - 1, f22s(x), f32s(x))
  
  CDFappx <- function(val) {
    tm <- val * invGP
    s1 <- floor(tm)
    if (s1 == tm) 
      return(c(s1, s1, gridInc, 0))
    s2 <- s1 + 1
    c(s1, s2, IntGrids2[s2] - val, val - IntGrids2[s1])
  }
  
  F_1_1s <- function(x, d1, d2) {
    if (x >= 2 * w - psn) 
      return(1)
    if (x < psn) 
      return(0)
    dd <- paste(d1, d2, sep = "")
    Fg <- switch(dd, `33` = F3131s_g, `32` = F3121s_g, `23` = F3121s_g, 
                 `22` = F2121s_g, `21` = F2111s_g, `12` = F2111s_g, 
                 `31` = F3111s_g, `13` = F3111s_g, `11` = F1111s_g)
    
    if (x < gridInc) 
      return(Fg[1] * x * invGP)
    tc <- CDFappx(x)
    sum(Fg[tc[1:2]] * tc[3:4]) * invGP
  }
  F_2s <- function(x, d) {
    if (x >= w - psn) 
      return(1)
    if (x < psn) 
      return(0)
    
    if (x < gridInc) 
      return(switch(d - 1, pchisq(x, 2)/pchisq(w, 2), F32s_g[1] * x * invGP))
    tc <- CDFappx(x)
    switch(d - 1, pchisq(x, 2)/pchisq(w, 2), sum(F32s_g[tc[1:2]] * tc[3:4]) * invGP)
  }
  F_1 <- function(x, d) {
    if (x >= w - psn) 
      return(1)
    if (x < psn) 
      return(0)
    if (x < gridInc) 
      return(switch(d, pchisq(x, 1)/pchisq(w, 1), F21_g[1] * x * invGP, 
                    F31_g[1] * x * invGP))
    tc <- CDFappx(x)
    switch(d, pchisq(x, 1)/pchisq(w, 1), sum(F21_g[tc[1:2]] * tc[3:4]) * invGP, 
           sum(F31_g[tc[1:2]] * tc[3:4]) * invGP)
  }
  
  
  f11_a <- function(x, a) (x < a) * dchisqM(x, 1)/pchisq(a, 1)
  f21_a <- function(x, a) pchisq(a - x, 1) * dchisqM(x, 1)/pchisq(a, 2)
  f31_a <- function(x, a) pchisq(a - x, 2) * dchisqM(x, 1)/pchisq(a, 3)
  f_1_a <- function(x, a, df) switch(df, f11_a(x, a), f21_a(x, a), f31_a(x, a))
  F_1_a <- function(x, a, df) {
    if (x <= 0) 
      return(0) else if (x > a) 
        return(1)
    ff <- switch(df, f11_a, f21_a, f31_a)
    integrate(function(y) ff(y, a), 0, x, subdivisions = 1e+05, rel.tol = reltol, 
              abs.tol = 0)$value
  }
  
  
  
  P_112 <- function(df) {
    IntegrandF <- function(x) f_1(x, df[2]) * (1 - F_2s(w - x, df[1]))
    cuhre(ndim = 1, ncomp = 1, integrand = IntegrandF, lower = 0, upper = w, 
          max.eval = 5e+07, flags = list(verbose = 0), rel.tol = reltol)$value
  }
  P_123 <- function(df) {
    IntegrandF <- function(x) f_1(x, df[1]) * (1 - F_1_1s(w - x, df[2], df[3]))
    cuhre(ndim = 1, ncomp = 1, integrand = IntegrandF, lower = 0, upper = w, 
          max.eval = 5e+07, flags = list(verbose = 0), rel.tol = reltol)$value
  }
  P_112122 <- function(df) {
    IntegrandFa <- function(x) f_1(x[1], df[1]) * f_1(x[2], df[2]) * 
      (1 - F_1_a(w - sum(x), w - x[1], df[1] - 1)) * F_1_a(w - sum(x), w - x[2], df[2] - 1)
    cuhre(ndim = 2, ncomp = 1, integrand = IntegrandFa, lower = rep(0, 2), 
          upper = rep(w, 2), max.eval = 5e+07, flags = list(verbose = 0), 
          rel.tol = reltol)$value
  }
  P_112113 <- function(df) {
    IntegrandF <- function(x) f_2(x, df[1]) * (1 - F_1(w - sum(x), df[2])) * 
      F_1(w - sum(x), df[3])
    cuhre(ndim = 2, ncomp = 1, integrand = IntegrandF, lower = rep(0, 2), 
          upper = rep(w, 2), max.eval = 5e+07, flags = list(verbose = 0), 
          rel.tol = reltol)$value
  }
  P_112123 <- function(df) {
    IntegrandFa <- function(x) f_1(x[1], df[1]) * f_1(x[2], df[2]) * 
      (1 - F_1_a(w - sum(x), w - x[1], df[1] - 1)) * F_1(w - sum(x), df[3])
    cuhre(ndim = 2, ncomp = 1, integrand = IntegrandFa, lower = rep(0, 2), 
          upper = rep(w, 2), max.eval = 5e+07, flags = list(verbose = 0), 
          rel.tol = reltol)$value
  }
  P_123112 <- function(df) {
    IntegrandFa <- function(x) f_1(x[1], df[1]) * f_1(x[2], df[2]) * 
      F_1_a(w - sum(x), w - x[1], df[1] - 1) * (1 - F_1(w - sum(x), df[3]))
    cuhre(ndim = 2, ncomp = 1, integrand = IntegrandFa, lower = rep(0, 2), 
          upper = rep(w, 2), max.eval = 5e+07, flags = list(verbose = 0), 
          rel.tol = reltol)$value
  }
  P_123124 <- function(df) {
    IntegrandF <- function(x) f_1(x[1], df[1]) * f_1(x[2], df[2]) * 
      (1 - F_1(w - sum(x), df[3])) * F_1(w - sum(x), df[4])
    cuhre(ndim = 2, ncomp = 1, integrand = IntegrandF, lower = rep(0, 2), 
          upper = rep(w, 2), max.eval = 5e+07, flags = list(verbose = 0), 
          rel.tol = reltol)$value
  }
  P_112122123 <- function(df) {
    IntegrandFa <- function(x) f_1(x[1], df[1]) * f_1(x[2], df[2]) * 
      (1 - F_1_a(w - sum(x), w - x[1], df[1] - 1)) * 
      F_1_a(w - sum(x), w - x[2], df[2] - 1) * F_1(w - sum(x), df[3])
    cuhre(ndim = 2, ncomp = 1, integrand = IntegrandFa, lower = rep(0, 2), 
          upper = rep(w, 2), max.eval = 5e+07, flags = list(verbose = 0), 
          rel.tol = reltol)$value
  }
  P_112123124 <- function(df) {
    IntegrandFa <- function(x) f_1(x[1], df[1]) * f_1(x[2], df[2]) * 
      (1 - F_1_a(w - sum(x), w - x[1], df[1] - 1)) * F_1(w - sum(x), df[3]) * 
      F_1(w - sum(x), df[4])
    cuhre(ndim = 2, ncomp = 1, integrand = IntegrandFa, lower = rep(0, 2), 
          upper = rep(w, 2), max.eval = 5e+07, flags = list(verbose = 0), 
          rel.tol = reltol)$value
  }
  P_123112122 <- function(df) {
    IntegrandFa <- function(x) f_1(x[1], df[1]) * f_1(x[2], df[2]) * 
      (1 - F_1(w - sum(x), df[3])) * F_1_a(w - sum(x), w - x[1], df[1] - 1) *
      F_1_a(w - sum(x), w - x[2], df[2] - 1)
    cuhre(ndim = 2, ncomp = 1, integrand = IntegrandFa, lower = rep(0, 2), 
          upper = rep(w, 2), max.eval = 5e+07, flags = list(verbose = 0), 
          rel.tol = reltol)$value
  }
  P_123122124 <- function(df) {
    IntegrandFa <- function(x) f_1(x[1], df[1]) * f_1(x[2], df[2]) * 
      (1 - F_1(w - sum(x), df[3])) * F_1_a(w - sum(x), w - x[2], df[2] - 1) * 
      F_1(w - sum(x), df[4])
    cuhre(ndim = 2, ncomp = 1, integrand = IntegrandFa, lower = rep(0, 2), 
          upper = rep(w, 2), max.eval = 5e+07, flags = list(verbose = 0), 
          rel.tol = reltol)$value
  }
  Pe_112113 <- function(df) {
    IntegrandF <- function(x) f33(x) * (1 - F_1(w - x[1] - x[2], df[2])) * 
      (1 - F_1(w - x[2] - x[3], df[3]))
    suave(ndim = 3, ncomp = 1, integrand = IntegrandF, lower = rep(0, 3), 
          upper = rep(w, 3), max.eval = 5e+07, flags = list(verbose = 0), 
          rel.tol = reltol)$value
  }
  Pe_112133 <- function(df) {
    IntegrandF <- function(x) f_2(x, df[1]) * (1 - F_1(w - sum(x), df[2])) * 
      (1 - F_2s(w - x[2], df[3]))
    cuhre(ndim = 2, ncomp = 1, integrand = IntegrandF, lower = rep(0, 2), 
          upper = rep(w, 2), max.eval = 5e+07, flags = list(verbose = 0), 
          rel.tol = reltol)$value
  }
  Pe_112134 <- function(df) {
    IntegrandF <- function(x) f_2(x, df[1]) * (1 - F_1(w - sum(x), df[2])) * 
      (1 - F_1_1s(w - x[2], df[3], df[4]))
    cuhre(ndim = 2, ncomp = 1, integrand = IntegrandF, lower = rep(0, 2), 
          upper = rep(w, 2), max.eval = 5e+07, flags = list(verbose = 0), 
          rel.tol = reltol)$value
  }
  Pe_112223 <- function(df) {
    IntegrandF <- function(x) f_2(x, df[2]) * (1 - F_2s(w - x[1], df[1])) * 
      (1 - F_1(w - sum(x), df[3]))
    cuhre(ndim = 2, ncomp = 1, integrand = IntegrandF, lower = rep(0, 2), 
          upper = rep(w, 2), max.eval = 5e+07, flags = list(verbose = 0), 
          rel.tol = reltol)$value
  }
  Pe_112233 <- function(df) {
    IntegrandF <- function(x) f_1(x, df[2]) * (1 - F_2s(w - x, df[1])) * 
      (1 - F_2s(w - x, df[3]))
    cuhre(ndim = 1, ncomp = 1, integrand = IntegrandF, lower = 0, upper = w, 
          max.eval = 5e+07, flags = list(verbose = 0), rel.tol = reltol)$value
  }
  Pe_112234 <- function(df) {
    IntegrandF <- function(x) f_1(x, df[2]) * (1 - F_2s(w - x, df[1])) * 
      (1 - F_1_1s(w - x, df[3], df[4]))
    cuhre(ndim = 1, ncomp = 1, integrand = IntegrandF, lower = 0, upper = w,
          max.eval = 5e+07, flags = list(verbose = 0), rel.tol = reltol)$value
  }
  Pe_123114 <- function(df) {
    IntegrandF <- function(x) f_2(x, df[1]) * (1 - F_1_1s(w - x[1], df[2], df[3])) * 
      (1 - F_1(w - sum(x), df[4]))
    cuhre(ndim = 2, ncomp = 1, integrand = IntegrandF, lower = rep(0, 2), 
          upper = rep(w, 2), max.eval = 5e+07, flags = list(verbose = 0), 
          rel.tol = reltol)$value
  }
  Pe_123144 <- function(df) {
    IntegrandF <- function(x) f_1(x, df[1]) * (1 - F_1_1s(w - x, df[2], df[3])) * 
      (1 - F_2s(w - x, df[4]))
    cuhre(ndim = 1, ncomp = 1, integrand = IntegrandF, lower = 0, upper = w, 
          max.eval = 5e+07, flags = list(verbose = 0), rel.tol = reltol)$value
  }
  Pe_123145 <- function(df) {
    IntegrandF <- function(x) f_1(x, df[1]) * (1 - F_1_1s(w - x, df[2], df[3])) * 
      (1 - F_1_1s(w - x, df[4], df[5]))
    cuhre(ndim = 1, ncomp = 1, integrand = IntegrandF, lower = 0, upper = w, 
          max.eval = 5e+07, flags = list(verbose = 0), rel.tol = reltol)$value
  }
  
  
  GenericP <- function(type, ldf) {
    HashKey <- paste(type, paste(ldf, collapse = "-"))
    if (is.na(GenericP_hm[[HashKey]])) {
      PFunc <- switch(type, `112` = P_112, `123` = P_123, `112122` = P_112122, 
                      `112113` = P_112113, `112123` = P_112123, 
                      `123112` = P_123112, `123124` = P_123124, 
                      `112122123` = P_112122123, `112123124` = P_112123124, 
                      `123112122` = P_123112122, `123122124` = P_123122124)
      GenericP_hm[[HashKey]] <<- PFunc(ldf)
    }
    GenericP_hm[[HashKey]]
  }
  GenericPe <- function(type, ldf) {
    HashKey <- paste(type, paste(ldf, collapse = "-"))
    if (is.na(GenericPe_hm[[HashKey]])) {
      PeFunc <- switch(type, e112113 = Pe_112113, e112133 = Pe_112133, 
                       e112134 = Pe_112134, e112223 = Pe_112223, 
                       e112233 = Pe_112233, e112234 = Pe_112234, 
                       e123114 = Pe_123114, e123144 = Pe_123144, 
                       e123145 = Pe_123145)
      GenericPe_hm[[HashKey]] <<- PeFunc(ldf)
    }
    GenericPe_hm[[HashKey]]
  }
  Ps <- function(l1, l2, labeldf) {
    cl <- intersect(l1, l2)
    ll1 <- l1[-which(l1 == cl)]
    ll2 <- l2[-which(l2 == cl)]
    
    componentf <- function(x, ll) {
      if (length(ll) == 1) 
        1 - F_1(w - x, labeldf[ll]) 
      else if (length(ll) == 2 && length(unique(ll)) == 2) 
        1 - F_1_1s(w - x, labeldf[ll[1]], labeldf[ll[2]]) 
      else if (length(ll) == 2 && length(unique(ll)) == 1) 
        1 - F_2s(w - x, labeldf[ll[1]])
    }
    
    if (sum(l1 == cl) + sum(l2 == cl) == 3) 
      IntegrandF <- function(x) f3132s(x) * componentf(x[1], ll1) * 
      componentf(x[2], ll2)
    else if (sum(l1 == cl) + sum(l2 == cl) == 2) 
      IntegrandF <- function(x) f_2(x, labeldf[cl]) * componentf(x[1], ll1) * 
      componentf(x[2], ll2)
    
    cuhre(ndim = 2, ncomp = 1, integrand = IntegrandF, lower = rep(0, 2), 
          upper = rep(w, 2), max.eval = 5e+07, flags = list(verbose = 0), 
          rel.tol = reltol)$value
  }
  IndepProb <- function(l) {
    if (length(unique(l)) == 3) {
      HashKey <- paste("123", paste(sort(labeldf[l]), collapse = "-"))
      if (is.na(IndepProb_hm[[HashKey]])) 
        IndepProb_hm[[HashKey]] <- P_123(labeldf[l])
      return(IndepProb_hm[[HashKey]])
    }
    dup <- l[duplicated(l)]
    cl <- c(dup, setdiff(l, dup))
    HashKey <- paste("112", paste(labeldf[cl], collapse = "-"))
    if (is.na(IndepProb_hm[[HashKey]])) 
      IndepProb_hm[[HashKey]] <<- P_112(labeldf[cl])
    IndepProb_hm[[HashKey]]
  }
  
  GenericP_hm <- hashmap("", 0)
  GenericPe_hm <- hashmap("", 0)
  IndepProb_hm <- hashmap("", 0)
  
  label <- initLabel(pstrct)
  labeldf <- sapply(1:length(unique(label)), function(x) sum(label == x))
  G3 <- initGroup3(label, pstrct)
  if (is.vector(G3)) 
    G3 <- matrix(G3, nrow = 1)
  Interaction2 <- initInteraction2(G3, label)
  Interaction1 <- initInteraction1(G3, label)
  Interaction0 <- initInteraction0(G3, label)
  Indep <- initIndep(G3, label)
  
  int2Verify <- function(n) ifelse(length(Interaction2[[n]]) == 1, 0, 
                                   length(Interaction2[[n]])/3)
  psum <- sapply(1:nrow(G3), function(n) length(Interaction0[[n]]) + 
                   length(Interaction1[[n]]) + int2Verify(n) + length(Indep[[n]]))
  if (!all(psum == 0:(nrow(G3) - 1))) 
    stop("Triplet initialization error")
  
  pchisqC3 <- pchisq(w, 3)
  pchisqC2 <- pchisq(w, 2)
  pchisqC1 <- pchisq(w, 1)
  
  ts1 <- proc.time()
  IntGrid <- seq(gridInc, ceiling(w * invGP) * gridInc, by = gridInc)
  IntGrids2 <- seq(gridInc, ceiling(2 * w * invGP) * gridInc, by = gridInc)
  
  F31_g <- vapply(IntGrid, function(upper) {
    integrate(f31, 0, upper, subdivisions = 1e+06, rel.tol = btol, abs.tol = 0)$value
  }, numeric(1))
  F21_g <- vapply(IntGrid, function(upper) {
    integrate(f21, 0, upper, subdivisions = 1e+06, rel.tol = btol, abs.tol = 0)$value
  }, numeric(1))
  F32s_g <- vapply(IntGrid, function(upper) {
    integrate(f32s, 0, upper, subdivisions = 1e+06, rel.tol = btol, abs.tol = 0)$value
  }, numeric(1))
  if (verbose) message("Stage I: 1/7 complete")
  
  cl <- makeCluster(nthread)
  registerDoParallel(cl)
  
  ig2 <- round(seq(1, length(IntGrids2) + 1, len = nthread + 1))
  F3131s_g <- unlist(foreach(k = 1:nthread) %dopar% {
    lseq <- IntGrids2[ig2[k]:(ig2[k + 1] - 1)]
    unlist(lapply(lseq, function(u) integrate(function(x) vapply(x, function(xx) 
      f31(xx) * F_1(u - xx, 3), numeric(1)), 0, u, subdivisions = 1e+06, abs.tol = 0)$value))
  })
  if (verbose) message("Stage I: 2/7 complete")
  
  F3121s_g <- unlist(foreach(k = 1:nthread) %dopar% {
    lseq <- IntGrids2[ig2[k]:(ig2[k + 1] - 1)]
    unlist(lapply(lseq, function(u) integrate(function(x) vapply(x, function(xx) 
      f31(xx) * F_1(u - xx, 2), numeric(1)), 0, u, subdivisions = 1e+06, abs.tol = 0)$value))
  })
  if (verbose) message("Stage I: 3/7 complete")
  
  F3111s_g <- unlist(foreach(k = 1:nthread) %dopar% {
    lseq <- IntGrids2[ig2[k]:(ig2[k + 1] - 1)]
    unlist(lapply(lseq, function(u) integrate(function(x) vapply(x, function(xx) 
      f31(xx) * F_1(u - xx, 1), numeric(1)), 0, u, subdivisions = 1e+06, abs.tol = 0)$value))
  })
  if (verbose) message("Stage I: 4/7 complete")
  
  F2121s_g <- unlist(foreach(k = 1:nthread) %dopar% {
    lseq <- IntGrids2[ig2[k]:(ig2[k + 1] - 1)]
    unlist(lapply(lseq, function(u) integrate(function(x) vapply(x, function(xx) 
      f21(xx) * F_1(u - xx, 2), numeric(1)), 0, u, subdivisions = 1e+06, abs.tol = 0)$value))
  })
  if (verbose) message("Stage I: 5/7 complete")
  
  F2111s_g <- unlist(foreach(k = 1:nthread) %dopar% {
    lseq <- IntGrids2[ig2[k]:(ig2[k + 1] - 1)]
    unlist(lapply(lseq, function(u) integrate(function(x) vapply(x, function(xx) 
      f21(xx) * F_1(u - xx, 1), numeric(1)), 0, u, subdivisions = 1e+06, abs.tol = 0)$value))
  })
  if (verbose) message("Stage I: 6/7 complete")
  
  F1111s_g <- unlist(foreach(k = 1:nthread) %dopar% {
    lseq <- IntGrids2[ig2[k]:(ig2[k + 1] - 1)]
    unlist(lapply(lseq, function(u) integrate(function(x) vapply(x, function(xx) 
      f11(xx) * F_1(u - xx, 1), numeric(1)), 0, u, subdivisions = 1e+06, abs.tol = 0)$value))
  })
  if (verbose) message("Stage I: 7/7 complete")
  
  stopCluster(cl)
  if (verbose) message("Stage I running time: ", round((proc.time() - ts1)[3]), "s")
  
  
  ts2 <- proc.time()
  GridCProb <- prod(unlist(lapply(labeldf, function(x) pchisq(w, x))))
  Pu <- 1 - GridCProb
  for (i in 1:length(Interaction2)) {
    CP <- GenericP(attributes(Interaction2[[i]])$type, 
                   labeldf[attributes(Interaction2[[i]])$idxdist])
    Pu <- Pu + GridCProb * CP
  }
  if (verbose) 
    message("Stage II: 1/4 complete")
  
  Eu_int1 <- 0
  for (i in 1:nrow(G3)) {
    if (length(Interaction1[[i]]) == 0) 
      next
    for (j in Interaction1[[i]]) {
      CP <- GenericPe(attributes(j)$type, labeldf[attributes(j)$idxdist])
      Eu_int1 <- Eu_int1 + GridCProb * CP
    }
  }
  if (verbose) 
    message("Stage II: 2/4 complete")
  
  cl <- makeCluster(nthread)
  registerDoParallel(cl)
  Pu3 <- unlist(foreach(i = 1:nrow(G3), .export = c("cuhre")) %dopar% {
    if (length(Interaction0[[i]]) == 0) 
      return(0)
    temp <- 0
    for (j in Interaction0[[i]]) temp <- temp + Ps(label[G3[i, ]], label[j], labeldf)
    GridCProb * temp
  })
  stopCluster(cl)
  Eu_int0 <- sum(Pu3)
  if (verbose) 
    message("Stage II: 3/4 complete")
  
  Eu_indep <- 0
  for (i in 1:nrow(G3)) {
    if (length(Indep[[i]]) == 0) 
      next
    for (j in Indep[[i]]) {
      Eu_indep <- Eu_indep + GridCProb * IndepProb(label[G3[i, ]]) * IndepProb(label[j])
    }
  }
  if (verbose) 
    message("Stage II: 4/4 complete")
  if (verbose) 
    message("Stage II running time: ", round((proc.time() - ts2)[3]), "s")
  
  list(Pu = min(1, Pu), Pl = Pu - (Eu_int1 + Eu_int0 + Eu_indep))
}