# gene level statistics is calculated from here
# the default method for integration gene level information is OSP

geneStats <- function(T, method="OSP") {
   pvalue <- 0

   if (method == "Stouffer") {
      pvalue <- pnorm(sum(T), 0, sqrt(length(T)), lower.tail =FALSE)
   } else if (method == "OSP") {
      p <- pnorm(T, lower.tail = TRUE)
	  pvalue <- pchisq(-2*sum(log(p)), 2*length(p), lower.tail = TRUE)
   } else if (method == "Fisher") {
	  p <- pnorm(T, lower.tail = FALSE)
	  pvalue <- pchisq(-2*sum(log(p)), 2*length(p), lower.tail = FALSE)
   } else if (method == "maxP") {
	  pvalue <- pnorm(max(T), lower.tail = FALSE)
   }

   return (pvalue)
}
