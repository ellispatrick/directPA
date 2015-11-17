#' Direction Analysis for Pathways and Kinases
#'
#' The main function of direction Analysis. This function takes in a matrix of test statistics with 
#' two (2-dimensional space) or three (3-dimensional space) columns, the direction of interests, and 
#' the annotation list such as pathway annotation or kinase-substrate annotation, and test for 
#' enrichment of pathways or kinases on the specified direction.
#' 
#' @usage directPA(Tc, direction, annotation, minSize=5, gene.method="OSP", 
#' path.method="Stouffer", visualize=TRUE)
#' 
#' @param Tc a numeric matrix. The columns are genes or phosphorylation sites and the columns are 
#' treatments vs control statistics.
#' @param direction the direction to be tested for enrichment. Either specified as a degree for 
#' two-dimensional analysis or as contrast (in a triplet) for three-dimensional analysis.
#' @param annotation a list with names correspond to pathways or kinases and elements correspond to 
#' genes or substrates belong to each pathway or kinase, respectively.
#' @param minSize the size of annotation groups to be considered for calculating enrichment. Groups 
#' that are smaller than the minSize will be removed from the analysis.
#' @param gene.method the method to be used for integrating statistics across treatments for each gene
#' or phosphorylation site.  Available methods are Stouffer, OSP, Fisher, and maxP. 
#' Default method is OSP.
#' @param path.method the method to be used for integrating statistics of all genes or phosphorylation
#' sites that belongs to a pathway or kinase. Available methods are Stouffer, OSP, Fisher, and maxP. 
#' Default method is Stouffer. 
#' @param visualize whether to visualize the 
#' @return a list of enrichment for pathways or kinases
#' @export
#' @examples
#' 
#' ## I. for kinase anlaysis on a phosphoproteomics dataset (Humphrey et al. Cell Metab., 2013)
#' # load the phosphoproteomics dataset
#' data(HEK)
#' 
#' # load the kinase-substrate annoations
#' data(PhosphoSite)
#' 
#' # direction pathway analysis in 2-dimensional space. Implemented as rotating by degree 
#' # (1) test combined effect of Torin1 and Rapamycin vs insul both on "down-regulation"
#' # (180 degree to original direction)
#' kst1 <- directPA(Tc=HEK, direction=pi, annotation=PhosphoSite.mouse)
#' kst1[order(unlist(kst1[,1])),][1:20,]
#' 
#' # (2) test combined effect of Torin1 and Rapamycin vs insul on "no change and down-regulation"
#' # (135 degree to the original direction) 
#' kst2 <- directPA(Tc=HEK, direction=pi*3/4, annotation=PhosphoSite.mouse)
#' kst2[order(unlist(kst2[,1])),][1:20,]
#' 
#' # (3) test combined effect of Torin1 and Rapamycin vs insul on "down-regulation and no change"
#' # (225 degree to the original direction) 
#' kst3 <- directPA(Tc=HEK, direction=pi*5/4, annotation=PhosphoSite.mouse)
#' kst3[order(unlist(kst3[,1])),][1:20,]
#' 
#' 
#' ## II. for pathway analysis on a proteomics dataset (Yang et al. Bioinformatics, 2014)
#' # load the proteomics dataset
#' data(PM)
#' 
#' # load pathway annotations
#' data(Pathways)
#' 
#' # display reactome pathways. Could be replaced by any other pathway databases
#' Pathways.reactome[1:5]
#' 
#' # direction pathway analysis in 3-dimensional space. Implemnted as rotating by contrast
#' # (1) test combined effect of all 3 treatments (stimulation and inhibitions) vs control (basal) 
#' # on the original direction.
#' gst1 <- directPA(Tc=PM, direction=c(1,1,1), annotation=Pathways.reactome)
#' gst1[order(unlist(gst1[,1])),][1:20,]
#' 
#' # (2) test combined effect of all 3 treatments vs controls on direction c(1,-1, 0)
#' # this rotates Ins by 0 degree, Wmn by 90 degree, and MK by 45 degree.
#' gst2 <- directPA(Tc=PM, direction=c(1,-1,0), annotation=Pathways.reactome)
#' gst2[order(unlist(gst2[,1])),][1:20,]
#' 
#' # (3) test combined effect of all 3 perturbations vs controls on direction c(1,-1, 1)
#' # this rotates Ins by 0 degree, Wmn by 90 degree, and MK by 0 degree.
#' gst3 <- directPA(Tc=PM, direction=c(1,-1,1), annotation=Pathways.reactome)
#' gst3[order(unlist(gst3[,1])),][1:20,]
#'
directPA <- function(Tc, direction, annotation, minSize=5, gene.method="OSP", path.method="Stouffer", visualize=TRUE, ...){
  ## spherical coordinates for three-dimensional rotation
  if (length(direction) == 3) {
    
    # step 1. convert statistics into z-scores
    Tc.zscores <- apply(Tc, 2, function(x){qnorm(rank(x)/(nrow(Tc)+1))})
    
    # step 2. rotate z-scores
    Tc.rotated <- rotate3d(Tc.zscores, direction)
    
    # step 3. integrate statistics across treatments
    gene.pvalues <- apply(Tc.rotated, 1, geneStats, gene.method)
    
    if (visualize == TRUE) {
      HC = rainbow(length(gene.pvalues)*1.2)
      plot3d(Tc, col=HC[rank(gene.pvalues)], size=5, ...)
      abclines3d(x=0, y=0, z=0, a=diag(3), col="black", lwd=3)
      abclines3d(x=0, a=direction, col="pink", lwd=5)
    }
    
    # step 4. integrate statistics for pathways or kinases
    gene.zscores <- qnorm(gene.pvalues, lower.tail = FALSE)
    gst <- t(sapply(annotation, pathwayStats, gene.zscores, minSize=5, path.method))
    
    result <- list()
    result$gene.pvalues <- gene.pvalues
    result$gst <- gst
    return(result)
  }
  
  ## polar coordinates for two-dimensional rotation
  if (length(direction) == 1) {
    
    # step 1. convert statistics into z-scores
    Tc.zscores <- apply(Tc, 2, function(x){qnorm(rank(x)/(nrow(Tc)+1))})
    
    # step 2. rotate z-scores
    Tc.rotated <- rotate2d(Tc.zscores, direction)
    
    # step 3. integrate statistics across treatments
    gene.pvalues <- apply(Tc.rotated, 1, geneStats, gene.method)
    
    if (visualize == TRUE) {
      HC = rainbow(length(gene.pvalues)*1.2)
      plot(Tc, col=HC[rank(gene.pvalues)], pch=16, ...)
      abline(v = 0,h = 0, lty=2, col="gold")
      abline(a=0, b=1, col="darkgreen", lty=2)
      abline(a=0, b=-1, col="darkgreen", lty=2)
    }
    
    # step 4. integrate statistics for pathways or kinases
    gene.zscores <- qnorm(gene.pvalues, lower.tail = FALSE)
    gst <- t(sapply(annotation, pathwayStats, gene.zscores, minSize=5, path.method))

    result <- list()
    result$gene.pvalues <- gene.pvalues
    result$gst <- gst
    return(result)
  }
}
