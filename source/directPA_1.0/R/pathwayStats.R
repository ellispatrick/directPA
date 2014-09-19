# Pathway level statistics is calculated from here
# The default method for integrating pathway level information is Stouffer

pathwayStats = function(PGs, T, minSize, method="Stouffer"){
  pvalue <- 0
  Z <- T[names(T) %in% PGs]

  if (length(Z) >= minSize) {
    if (method == "Stouffer") {
	   pvalue <- pnorm(sum(Z), 0, sqrt(length(Z)), lower.tail=FALSE)
	} else if (method == "OSP") {
       p <- pnorm(Z, lower.tail = TRUE)
	   pvalue <- pchisq(-2*sum(log(p)), 2*length(p), lower.tail = TRUE)
	} else if (method == "Fisher") {
	   p <- pnorm(Z, lower.tail = FALSE)
	   pvalue <- pchisq(-2*sum(log(p)), 2*length(p), lower.tail = FALSE)
    } else if (method == "maxP") {
	   pvalue <- pnorm(max(Z), lower.tail = FALSE)
    }
  } else { 
    pvalue <- NA
  }
  
  result <- list()
  result$pvalue <- pvalue
  result$size <- length(Z)
  return (result)
}
