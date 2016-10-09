\name{directPA}
\alias{directPA}
\title{
Direction Pathway Analysis
}
\description{
The main function of direction pathway analysis package.
}
\usage{
directPA(Tc, direction, pathway.list, minSize = 5, gene.method = "OSP", 
path.method = "Stouffer", visualize = FALSE)
}
\arguments{
  \item{Tc}{
The matrix of test statistics where rows correspond to proteins/genes and columns correspond to perturbations compared to controls.
}
  \item{direction}{
The direction to be tested. For two-dimensional rotation (experiments with two perturbations), directions are specified
as degree. For three-dimensional rotation (experiments with three perturbations), directions are specified as contrast.
See examples below for more details.
}
  \item{pathway.list}{
A pathway database in list format. See example below for details.
}
  \item{minSize}{
The minimum size of the pathway to be included in the report.
}
  \item{gene.method}{
A p-value integration method for protein/gene level combination. That is integrating information for each protein/gene accross
all perturbations (i.e. for each row of matrix Tc, integrating accross the columns). Available methods are Stouffer, OSP, Fisher, and maxP. 
Default method is OSP.
}
  \item{path.method}{
A p-value integration method for pathway level combination. That is integrating information of protein/gene included in a pathway, given
a pathway database. Available methods are Stouffer, OSP, Fisher, and maxP. Default method is Stouffer. 
}
  \item{visualize}{
A boolean value indicating whether to visualize the protein/gene level integration statistics in a scatter plot.
}
}
\value{
Return a matrix of test statistics for all given pathway that passed the minimum size cutoff on a specified test direction.
}
\author{
Pengyi Yang & Ellis Patrick
}
\references{
Pengyi Yang, Ellis Patrick, Shi-Xiong Tan, Daniel J. Fazakerley, James Burchfield, Christopher Gribben, Matthew J. Prior, David E. James, Yee Hwa Yang, Direction pathway analysis of large-scale proteomics data reveals novel features of the insulin action pathway, submitted.
}

\examples{

# Load the example dataset.
data(PM)

# Load the reactome pathway.
data(reactome)

# Display reactome pathways. Could be replaced by any other pathway databases.
reactome.list

# Combine test statistics of each perturbation vs control into a matrix of test statistics.
Tc = cbind(Ins,Wmn,MK)

## Direction pathway analysis in two-dimensional space. Implemented as rotating by degree. 
# (1) Testing combined effect of Ins and Wmn vs controls on original direction.
gst1 <- directPA(Tc[,c(1,2)], direction=0, pathway.list=reactome.list)
# Display the top 20 pathways ranked by integrated p-value.
gst1[order(unlist(gst1[,1])),][1:20,]

# (2) Testing combined effect of Ins and MK vs controls on the pi/2 direction.
gst2 <- directPA(Tc[,c(1,3)], direction=pi/2, pathway.list=reactome.list)
gst2[order(unlist(gst2[,1])),][1:20,]

## Direction pathway analysis in three-dimensional space. Implemnted as rotating by contrast.
# (1) Testing combined effect of all 3 perturbations vs controls on the original direction.
gst3 <- directPA(Tc, direction=c(1,1,1), pathway.list=reactome.list)
gst3[order(unlist(gst3[,1])),][1:20,]

# (2) Testing combined effect of all 3 perturbations vs controls on direction c(1,-1, 0). 
# This rotates Ins by 0 degree, Wmn by 90 degree, and MK by 45 degree
gst4 <- directPA(Tc, direction=c(1,-1,0), pathway.list=reactome.list)
gst4[order(unlist(gst4[,1])),][1:20,]

# (3) Testing combined effect of all 3 perturbations vs controls on direction c(1,-1, 1). 
# This rotates Ins by 0 degree, Wmn by 90 degree, and MK by 0 degree
gst5 <- directPA(Tc, direction=c(1,-1,1), pathway.list=reactome.list)
gst5[order(unlist(gst5[,1])),][1:20,]

}
