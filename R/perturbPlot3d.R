#' Perturbation Plot 3D
#'
#' This function takes in a matrix of test statistics with two columns (3-dimensional space) and the 
#' annotation list such as pathway annotation or kinase-substrate annotation, and visualize the enrichment
#' of pathways or kinases in direction specific manner.
#' 
#' @usage perturbPlot3d(Tc, annotation, minSize=5, ...)
#' @param Tc a numeric matrix. The columns are genes or phosphorylation sites and the columns are treatments 
#' vs control statistics.
#' @param annotation a list with names correspond to pathways or kinases and elements correspond to genes or
#' substrates belong to each pathway or kinase, respectively.
#' @param minSize the size of annotation groups to be considered for calculating enrichment. Groups 
#' that are smaller than the minSize will be removed from the analysis.
#' @param ... parameters for controling the plot.
#' @return a list of coordinates for pathways or kinases
#' @export
#' 
perturbPlot3d <- function(Tc, annotation, minSize=5, ...) {
  
  # step 1. convert statistics into z-scores
  Tc.zscores <- apply(Tc, 2, function(x){qnorm(rank(x)/(nrow(Tc)+1))})
  
  # step 2. filter the groups that are smaller than the minimun cutoff
  DE = lapply(annotation, function(x){
    if(sum(rownames(Tc.zscores) %in% x) >= minSize) {
      X <- Tc.zscores[rownames(Tc.zscores)%in%x,]
      n = nrow(X)
      Z1 = sum(X[,1])/sqrt(n)
      Z2 = sum(X[,2])/sqrt(n)
      Z3 = sum(X[,3])/sqrt(n)
      list(Z1=Z1, Z2=Z2, Z3=Z3)
    }
  })
  
  # step3. filter DE that has 0 element
  DE <- DE[which(sapply(DE, length) != 0)]

  # step4. visualise
  df <- data.frame(do.call(rbind, lapply(DE, function(x) unlist(unname(x)))))
  colnames(df) <- c("Z1", "Z2", "Z3")
  df$pathway <- rownames(df)
  
  my_col = colorRampPalette(rainbow(12))(100)
  t <- list(family = "sans serif", size = 14, color = "black")
  
  p <- plotly::plot_ly(df, x=~Z1, y=~Z2, z=~Z3, size=5, text=df$pathway)
  p <- plotly::add_markers(p)
  p <- plotly::add_text(p, textfont = t, textposition = "top right")
  p <- plotly::layout(p, scene = list(xaxis = list(title = colnames(df)[[2]]),
                                yaxis = list(title = colnames(df)[[3]]),
                                zaxis = list(title = colnames(df)[[4]])))
  print(p)
  
  
  ## return the results
  result <- list()
  result$Z1 <- df$Z1
  result$Z2 <- df$Z2
  result$Z3 <- df$Z3
  return(result)
}