\name{geneStats}
\alias{geneStats}
\title{
Gene Level Statistics
}
\description{
Calculating protein/gene level statistics accross multiple perturbations.
}
\usage{
geneStats(T, method = "OSP")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{T}{
A matrix of test statistics (normally after rotation), where rows correspond to proteins/genes and columns correspond to perturbations.
}
  \item{method}{
The p-value integration method for combining accross multiple perturbations. 
Available methods are Stouffer, OSP, Fisher, and maxP. The default method is OSP.
}
}
\value{
Return the integration p-value for a protein/gene.
}
\author{
Pengyi Yang & Ellis Patrick
}


\examples{

# Load the example data.
data(PM)


## (1) For three perturbatins vesus controls, use three dimentional rotation.
Tc <- cbind(Ins, Wmn, MK)

# Rotate the test statistics matrix to the direction of interest.
Tc.rotated <- rotate3Sphere(Tc, direction=c(1,-1,0))

# Integrate protein/gene statistics accross perturmations.
gene.pvalues <- apply(Tc.rotated, 1, geneStats)

# Visualize the rotation statistics.
library(rgl)
HC <- rainbow(length(gene.pvalues)*1.2)
plot3d(Tc, col=HC[rank(gene.pvalues)], size=5)
abclines3d(x=0, y=0, z=0, a=diag(3), col="black", lwd=3)
abclines3d(x=0, a=c(1,-1,0), col="pink", lwd=5)


## (2) For two perturbations vesus controls, use two dimentional rotation.
Tc <- cbind(Ins, Wmn)
Tc.rotated <- rotate2Sphere(Tc, degree=pi/4)
gene.pvalues <- apply(Tc.rotated, 1, geneStats)

# Visualize the rotation statistics.
HC = rainbow(length(gene.pvalues)*1.2)
plot(Tc, col=HC[rank(gene.pvalues)], pch=16)
abline(v = 0,h = 0, lty=2, col="gold")  
}


