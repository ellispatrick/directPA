#' Gene Level Statistics
#' 
#' Takes a vector of statistics with each element corresponds to a treatment vs control comparison, 
#' and calculates a combined statistics accross multiple treatments.
#' 
#' @param T a vector of statistics (z-scores) with each element correspond to a treatment vs control comparison.
#' @param method the p-value integration method for combining accross multiple treatments. Available methods are 
#' Stouffer, OSP, Fisher, and maxP. The default method is OSP.
#' @return a p-value after integration across treatments.
#' @export
#' 
#' @examples
#' # Load the example data.
#' data(PM)
#' 
#' # (1) For three perturbatins vesus controls, use three dimentional rotation.
#' Tc <- cbind(Ins, Wmn, MK)
#' 
#' 
#' 
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
