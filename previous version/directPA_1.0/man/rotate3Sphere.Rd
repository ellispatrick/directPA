\name{rotate3Sphere}
\alias{rotate3Sphere}
\title{
Spherical Coordinates Rotation
}
\description{
Rotate to the direction of interest in spherical coordinates by contrast (i.e 1, -1, -1).
}
\usage{
rotate3Sphere(T, direction = c(1, 1, 1))
}
\arguments{
  \item{T}{
A matrix of test statistics. The rows correspond to protein/genes and the columns correspond to perturbations. The number of columns should be 3. 
}
  \item{direction}{
The angle to be rotated in contrast format (see example below).
}
}
\value{
A transformed matrix with respect to the direction of interest.
}
\author{
Pengyi Yang & Ellis Patrick
}

\examples{
# Load the example data.
data(PM)

# Combine three perturbations into a single matrix. 
# The matrix contains test statistics of each protein w.r.t each perturbations.
Tc <- cbind(Ins, Wmn, MK)

# Rotate the matrix by contrast 1, -1, -1 in spherical coordinates.
Tc.rotated <- rotate3Sphere(Tc, direction = c(1, -1, -1))
}
