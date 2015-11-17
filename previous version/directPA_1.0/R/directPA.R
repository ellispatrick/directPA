##### 

directPA <- function(Tc, direction, pathway.list, minSize=5, gene.method="OSP", path.method="Stouffer", visualize=FALSE){
   # spherical coordinates for three-dimensional rotation
   if (length(direction) == 3) {
      Tc.rotated <- rotate3Sphere(Tc, direction)
	  gene.pvalues <- apply(Tc.rotated, 1, geneStats, gene.method)
	  
	  if (visualize == TRUE) {
	         library(rgl)
                 HC = rainbow(length(gene.pvalues)*1.2)
		 plot3d(Tc, col=HC[rank(gene.pvalues)], size=5)
		 abclines3d(x=0, y=0, z=0, a=diag(3), col="black", lwd=3)
		 abclines3d(x=0, a=direction, col="pink", lwd=5)
	  }
	  
	  gene.zscores <- qnorm(gene.pvalues, lower.tail = FALSE)
	  gst <- t(sapply(pathway.list, pathwayStats, gene.zscores, minSize=5, path.method))
	  return(gst)
   }
   
   # polar coordinates for two-dimensional rotation
   if (length(direction) == 1) {
      Tc.rotated <- rotate2Sphere(Tc, direction)
	  gene.pvalues <- apply(Tc.rotated, 1, geneStats, gene.method)
	  
	  if (visualize == TRUE) {
         HC = rainbow(length(gene.pvalues)*1.2)
         plot(Tc, col=HC[rank(gene.pvalues)], pch=16)
		 abline(v = 0,h = 0, lty=2, col="gold")  
	  }
	  
	  gene.zscores <- qnorm(gene.pvalues, lower.tail = FALSE)
	  gst <- t(sapply(pathway.list, pathwayStats, gene.zscores, minSize=5, path.method))
	  return(gst)
   }
}
