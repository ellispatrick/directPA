\name{rotate2Sphere}
\alias{rotate2Sphere}
\title{
Polar Coordinates Rotation
}
\description{
Rotate to the direction of interest in polar coordinates by degree (i.e pi/4). 
}
\usage{
rotate2Sphere(T, degree = 0)
}
\arguments{
  \item{T}{
A matrix of test statistics. The rows correspond to protein/genes and the columns correspond to perturbations. The number of colmuns should be 2.
}
  \item{degree}{
The degree to be rotated. 
}
}
\value{
A transformed matrix with respect to the direction of interest.
}
\author{
Pengyi Yang & Ellis Patrick
}

\examples{
# Load the example dataset.
data(PM)

# Combine two perturbations into a single matrix. 
# The matrix contains test statistics of each protein w.r.t each perturbation.
Tc <- cbind(Ins, Wmn)

# Rotate the matrix by pi/2 in polar coordinates.
Tc.rotated <- rotate2Sphere(Tc, degree = pi/2)
}


