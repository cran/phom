\name{phom-package}
\alias{phom-package}
\alias{phom}
\docType{package}
\title{
Persistent Homology in R
}
\description{
This package computes persistent homology of filtered simplicial complexes, 
and provides facilities for constructing complexes from geometric data. It has
been developed by members of the computational topology workgroup at Stanford University.
More information may be found at \url{http://comptop.stanford.edu}.
}
\details{
\tabular{ll}{
Package: \tab phom\cr
Type: \tab Package\cr
Version: \tab 1.0.3\cr
Date: \tab 11/5/11\cr
License: \tab BSD\cr
LazyLoad: \tab yes\cr
}

The pipeline for computing the persistent homology consists of two stages: the construction
of a filtered simplicial complex, and the computation of the homology of the filtered 
complex. This package has two methods for creating filtered simplicial complexes from
geometric data:

\itemize{
\item Vietoris-Rips mode: This uses the Vietoris-Rips construction to build a sequence of 
simplicial complexes on the underlying point-set. Note that this construction
may produce complexes that have a very large number of simplices in them.
\item Lazy-Witness mode: This uses the lazy-witness construction to build a filtered
simplicial complex on a subset of the points. This construction has the advantage of 
reducing the noise in the resulting set of intervals, and is considerably faster than
the Vietoris-Rips filtration. 
}

The algorithms operate on two types of data:

\itemize{
\item Euclidean Data: This assumes that the point set is represented as a matrix, where
each point is a row. The number of columns is the dimensionality of the dataset. A wide
variety of metrics are available on this sort of data.
\item Distance Matrix: For more general metric spaces, one may specify the interpoint
distances between all of the points. Note that the algorithms do not check that this
indeed specifies a metric.
}

}
\author{
Andrew Tausz \email{atausz@stanford.edu}
}
\references{

Gunnar Carlsson. Topology and Data. Bull. Amer. Math. Soc. 46 (2009), 255-308.

}
%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
%%~~ the R documentation directory ~~
\keyword{ package }

\examples{

### Example 1

library(phom)

x <- runif(100)
y <- runif(100)
points <- t(as.matrix(rbind(x, y)))

max_dim <- 2
max_f <- 0.2

intervals <- pHom(points, max_dim, max_f, metric="manhattan")

plotPersistenceDiagram(intervals, max_dim, max_f, 
	title="Random Points in Cube with l_1 Norm")


### Example 2

library(phom)

t <- 2 * pi * runif(100)
x <- cos(t); y <- sin(t)
X <- t(as.matrix(rbind(x, y)))

max_dim <- 1
max_f <- 0.6

intervals <- pHom(X, max_dim, max_f)

plotPersistenceDiagram(intervals, max_dim, max_f,
	title="Random Points on S^1 with Euclidean Norm")


}


