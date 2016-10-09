\name{pathwayStats}
\alias{pathwayStats}
\title{
Pathway Level Statistics
}
\description{
Calculating pathway level statistics by combining protein/gene belong to this pathway.
}
\usage{
pathwayStats(PGs, T, minSize, method = "Stouffer")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PGs}{
An array of names indicating proteins/genes belong to a given pathway.
}
  \item{T}{
An array of test statistics of all proteins/genes. This should be the zsore transformed output of geneStats.
}
  \item{minSize}{
The minimum size of a pathway to be reported.
}
  \item{method}{
A p-value integration method for combining proteins/genes in a given pathway.
}
}
\value{
  \item{pvalue }{Integrated pvalue of the given pathway}
  \item{size }{size of the given pathway}
}
\author{
Pengyi Yang & Ellis Patrick
}

\examples{

# Load the example data.
data(PM)
Tc = cbind(Ins, Wmn, MK)

# Load reactome pathway.
data(reactome)

# Rotate the test statistics matrix to the direction of interest.
Tc.rotated <- rotate3Sphere(Tc, direction=c(1, -1, 0))

# Integrate protein/gene statistics accross perturmations.
gene.pvalues <- apply(Tc.rotated, 1, geneStats)

# Transform the pvlaue computed from geneStats to zsocres.
gene.zscores <- qnorm(gene.pvalues, lower.tail = FALSE)

# Compute statistics for all reactome pathway
gst <- t(sapply(reactome.list, pathwayStats, gene.zscores, minSize=5))

}
