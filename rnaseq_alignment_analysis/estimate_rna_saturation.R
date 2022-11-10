# ---------------------------------------------------------------
# Author.... Merritt Khaipho-Burch
# Contact... mbb262@cornell.edu
# Date...... 2022-01-05 
#
# Description 
#   - estimate rna saturation and extract counts version 1
#   - taken from: 
# ---------------------------------------------------------------


estimate_saturation <-
  function(counts, max_reads=Inf,
           method="sampling",
           ndepths=6, nreps=5,
           min_counts=1, min_cpm=NULL,
           verbose=FALSE) {
    if (sum(!is.null(min_counts), !is.null(min_cpm)) != 1)
      stop("One of min_counts or min_cpm must be specified, but not both.")
    method <- match.arg(method, choices=c("division", "sampling"))
    
    counts <-
      extract_counts(counts, return_class="matrix") # extract counts and/or convert to matrix
    
    readsums <- colSums(counts)
    max_reads <- min(max(readsums), max_reads)
    depths <- round(seq(from=0, to=max_reads, length.out=ndepths+1))
    saturation <-
      data.frame(sample=as.vector(sapply(colnames(counts), FUN=rep, times=ndepths+1)),
                 depth=rep(depths, time=ncol(counts)))
    sat.estimates <- as.numeric(rep(NA, ncol(counts) * length(depths)))
    counter <- 0
    if (method=="sampling")
      sat.var.estimates <- as.numeric(rep(NA, ncol(counts) * length(depths)))
    for (i in 1:ncol(counts)) {
      if (verbose) cat("Working on library", i, "of", ncol(counts), "\n")
      
      # adjust to min_cpm if specified
      if (!is.null(min_cpm)) {
        min_counts.lib <- readsums[i] / 1000000
      } else {
        min_counts.lib <- min_counts
      }
      
      probs <- counts[,i, drop=TRUE] / readsums[i] # calculate gene probabilities for the library
      probs <- probs[probs > 0] # zero counts add nothing but computational time!
      ngenes <- length(probs)
      for (j in depths) {
        counter <- counter + 1
        if (j == 0) {
          sat.estimates[counter] <- 0
          if (method=="sampling")
            sat.var.estimates[counter] <- 0
        } else if (j > readsums[i]) {
          sat.estimates[counter] <- NA
          if (method=="sampling")
            sat.var.estimates[counter] <- NA
        } else if (method=="division") {
          sat.estimates[counter] <-
            sum((probs * j) >= min_counts.lib)
        } else if (method=="sampling") {
          est <- as.numeric(rep(NA, nreps))
          for (k in 1:nreps) {
            reads <- as.matrix(sample.int(n=ngenes, size=j, replace=TRUE, prob=probs))
            est[k] <- sum(bigtabulate::bigtable(reads, ccol=1) >= min_counts.lib)
          }
          sat.estimates[counter] <- mean(est)
          sat.var.estimates[counter] <- var(est)
        }
      }
    }
    saturation$sat <- sat.estimates
    if (method=="sampling")
      saturation$sat.var <- sat.var.estimates
    
    return(saturation)
  }


extract_counts <-
  function(counts, return_class=NULL) {
    if (inherits(counts, "EList")) {
      counts <- counts[["E"]]
    } else if (inherits(counts, "DGEList")) {
      counts <- counts[["counts"]]
    } else if (inherits(counts, "ExpressionSet") | inherits(counts, "eSet")) {
      if (is.environment(counts@assayData)) {
        counts <- counts@assayData[["exprs"]]
      } else if (is.matrix(counts@assayData)) {
        counts <- counts@assayData
      }
    } else if (!(is.matrix(counts) | is.data.frame(counts)))
      stop("Class of input \"counts\" object not recognized. Please check object class for compatability with this function.")
    
    if (!is.null(return_class)) counts <- as(counts, return_class)
    
    counts
  }
