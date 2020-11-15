#' Batch Direction Analysis in 2-dimentional space
#'
#' Rotate to the direction of interest in polar coordinates by degree (e.g. pi/4).
#' @usage directExplorer2d(Tc, annotation=NULL, gene.method="OSP", 
#' path.method="Stouffer", top=10, nd=8, ...)
#' @param Tc a numeric matrix with 2 columns. The rows are genes or phosphorylation sites and the columns 
#' are treatments vs control statistics.
#' @param annotation a list with names correspond to pathways or kinases and elements correspond to 
#' genes or substrates belong to each pathway or kinase, respectively.
#' @param gene.method the method to be used for integrating statistics across treatments for each gene
#' or phosphorylation site.  Available methods are Stouffer, OSP, Fisher, and maxP. 
#' Default method is OSP.
#' @param path.method the method to be used for integrating statistics of all genes or phosphorylation
#' sites that belongs to a pathway or kinase. Available methods are Stouffer, OSP, Fisher, and maxP. 
#' Default method is Stouffer. 
#' @param top the number of entries to be highlighted in the plot.
#' @param nd the number of directions to plot (4 or 8)
#' @param ... parameters for controlling the plot.
#' @return The the list of enrichment analysis in tables. 
#' @export
#' @examples
#' # load the phosphoproteomics dataset
#' data(HEK)
#' 
#' # load the kinase-substrate annoations
#' data(PhosphoSite)
#' 
#' # test enrichment on 8 directions in polar coordinate system.
#' bda <- directExplorer2d(Tc=HEK, annotation=PhosphoSite.mouse)
#' 
#' # the direction are denoted as follow for the two treatments vs control:
#' # ++: up-regulated in both treatments
#' # +*: up-regulated in the first treatment and unchanged in the second treatment
#' # +-: up-regulated in the first treatment and down-regulated in the second treatment
#' # *-: unchanged in the first treatment and down-regulated in the second treatment
#' # --: down-regulated in both treatments
#' # -*: down-regulated in the first treatment and unchanged in the second treatment
#' # -+: down-regulated in the first treatment and up-regulated in the second treatment
#' # *+: unchanged in the first treatment and up-regulated in the second treatment
#' 
#' # sort the most enriched phosphorylation sites and kinases on down-regulaiton from both 
#' # treatments (i.e. "--") and displa the top-10 entries
#' bda$gene.tab[order(bda$gene.tab[,"--"]),][1:10,]
#' bda$path.tab[order(bda$path.tab[,"--"]),][1:10,]
#'
directExplorer2d <- function(Tc, annotation=NULL, gene.method="OSP", path.method="Stouffer", top=10, nd=8, ...) {
  
  directionCode <- c('++','+*','+-','*-','--','-*','-+','*+')
  gene.tab <- matrix(NA, nrow(Tc), nd)
  rownames(gene.tab) <- rownames(Tc)
  if(nd != 8) {
    colnames(gene.tab) <- directionCode[c(1,3,5,7)]
  } else {
    colnames(gene.tab) <- directionCode
  }
  
  path.tab <- NULL
  if (!is.null(annotation)) {
    path.tab <- matrix(NA, length(annotation), nd+1)
    rownames(path.tab) <- names(annotation)
    if(nd != 8) {
      colnames(path.tab) <- c("size", directionCode[c(1,3,5,7)])
    } else {
      colnames(path.tab) <- c("size", directionCode)
    }
  }
  
  ds <- c()
  if (nd != 8) {
    ds <- c(0, 2, 4, 6)
  } else {
    ds <- 0:7
  }
  count <- 0
  
  for(i in ds){
    count <- count + 1
    
    # gene or phosphorylation site level
    Tc.zscores <- apply(Tc, 2, function(x){qnorm(rank(x)/(nrow(Tc)+1))})
    Tc.rotated <- rotate2d(Tc.zscores, pi/4*i)
    gene.pvalues <- apply(Tc.rotated, 1, geneStats, gene.method)
    gene.tab[,count] <- gene.pvalues
	  
	  # pathway or kinase level
    if (!is.null(annotation)) {
      gene.zscores <- qnorm(gene.pvalues, lower.tail = FALSE) 
      gst <- t(sapply(annotation, pathwayStats, gene.zscores, minSize=5, path.method))
	    if (i == 0) {
		    path.tab[,1] <- unlist(gst[,"size"])
	    }
	    path.tab[,(count+1)] <- unlist(gst[,"pvalue"])
    }
  }
  
  # highlight candidates
  plot(Tc, col="gray", pch=16, ...)
  abline(h=0, v=0, col="gold", lty=2)
  abline(a=0, b=1, lty = 2, col="darkgreen")
  color <- c("red", "orange3", "green4", "#169B48", "blue4", "#0F8AB5", "purple", "gray50")
  count <- 0
  for (i in ds) {
    count <- count + 1
     ids <- names(sort(gene.tab[,count])[1:top])
     points(Tc[ids,], col=color[i+1], pch=16)
     textxy(Tc[ids,1], Tc[ids,2], ids, col = color[i+1])
  }

  results <- list()
  results$gene.tab <- gene.tab
  results$path.tab <- path.tab
  return(results)
}

