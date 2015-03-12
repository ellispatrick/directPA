# directional kinase enrichment visualization
directEnrichPlot <- function(set.list, Tc, ...) {

   DKE = lapply(set.list, function(x){
      if(sum(rownames(Tc)%in%x) >= 5) {
         X <- Tc[rownames(Tc)%in%x,]
		 n = nrow(X)
		 Cbar = sum(X[,1])/sqrt(n)
		 Sbar = sum(X[,2])/sqrt(n)
		 list(Cbar=Cbar, Sbar=Sbar)
      }
    })

	# filter DKE that has 0 element
   DKE <- DKE[which(sapply(DKE, length) != 0)]
   Cbar <- unlist(sapply(DKE, function(x){x[1]}))
   Sbar <- unlist(sapply(DKE, function(x){x[2]}))
   
   # visualization
   plot(Cbar, Sbar, col="darkblue", pch=16, ...)
   textxy(Cbar,Sbar, names(DKE), col="black", cex=1)
   abline(v=0, h=0, col="gold", lty=2)
   abline(a=0, b=1, col="darkgreen", lty=2)
   abline(a=0, b=-1, col="darkgreen", lty=2)

   r <- ceiling(max(sqrt(Cbar^2 + Sbar^2)))
   for(i in seq(0, r, r/5)) {
      theta = seq(-3.14,3.14,0.05)
      lines(i*cos(theta),i*sin(theta),col = 'gray', type="l")
   }

   points(Cbar, Sbar, col="darkblue", pch=16, cex=2)

   ## return the results
   result <- list()
   result$DKE <- DKE
   result$Cbar <- Cbar
   result$Sbar <- Sbar
   return(result)
}
