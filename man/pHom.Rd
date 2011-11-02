
\name{pHom}
\alias{pHom}
\title{Compute Persistent Homology}
\description{

This function computes persistent homology on a given dataset. This is a two-step process 
which involves 1) creating a filtered simplicial complex on the dataset, and 2) computing
the persistent homology of the filtered complex.

It outputs a matrix with three columns. Each row in the output 
matrix corresponds to a persistence interval. The first column stores the dimension of the 
interval, the second stores the starting point, and the third stores the ending point.

The method provides two ways to construct a filtered simplicial complex on the data points:
the Vietoris-Rips construction, and the lazy-witness construction.

The Vietoris-Rips construction on a metric space, \eqn{(X, d)}, is performed as follows. 
\itemize{
\item For each \eqn{u \in X}, \eqn{u \in VR(X, 0)}
\item For each \eqn{u, v \in X}, \eqn{[u, v] \in VR(X, r)} if and only if \eqn{d(u, v) \leq r}
\item \eqn{VR(X, r)} is maximal among all simpicial complexes containing its \eqn{1}-skeleton. 
In other words, \eqn{[v_0, \ldots v_n] \in VR(X, r)} if and only if for all pairs, \eqn{[v_i, v_j] \in VR(X, r)}.
}

Note that the Vietoris-Rips construction can produce complexes that are very large. To circumvent this,
the witness construction may be used. This method produces complexes that are significantly smaller.


The lazy-witness construction is outlined as follows. Suppose that we have a metric space, \eqn{(X, d)}, and
a subset of points, \eqn{L \subset X}, which we refer to as the landmark set. The lazy-witness complex, 
parameterized by \eqn{\nu} is defined as follows:

\itemize{
\item The vertex set of \eqn{LW_{\nu}(X, L, r)} is \eqn{L}.
\item The edge \eqn{[u, v]} is in \eqn{LW_{\nu}(X, L, r)} if there exists a witness, \eqn{z}, such that
\eqn{\max(d(u, z), d(v, z)) \leq r + m(z)}. Here \eqn{m(z)} is the distance from \eqn{z} to the \eqn{\nu}-th
closest point in the landmark set.
\item A higher dimensional simplex is in \eqn{LW_{\nu}(X, L, r)} if all of its edges are.
}

The landmark set selection process is performed as follows. Initially a random point is added to \eqn{L}. 
Inductively, points are added which maximize the minimum distance to the existing landmark set. 
\deqn{
L_{i+1} = L_{i} \cup \arg \max_{x \in X} (\min_{y \in L_{i}} d(x, y))
}

The argument \code{mode} is used to select between the two different constructions, by either specififying 
\code{"vr"} or \code{"lw"}. 

}
\usage{
pHom(X, dimension, max_filtration_value, 
mode = "vr", metric = "euclidean", p = 2, 
landmark_set_size = 2 * ceiling(sqrt(length(X))), 
maxmin_samples = min(1000, length(X)))
}
\arguments{
\item{X}{A matrix which has one of the two following interpretations. In the case where \code{metric = "distance_matrix"}, \code{X} is required to be a
 square matrix whose entries indicate the distance between two points. The number of rows (and columns) equals the number of points in the dataset.
In the case where \code{metric} is something other than \code{"distance_matrix"}, \code{X} is required to be a matrix in which the rows are points in Euclidean
space. The number of columns is the dimensionality of the dataset. Note that for different choices of the \code{metric} argument, the points will
be endowed with different metrics.
}
\item{dimension}{The maximum dimension to compute persistent homology to.}
\item{max_filtration_value}{The maximum filtration value to use in constructing the filtered complex.}
\item{mode}{This indicates the type of filtration to use. The two possible choices are \code{"vr"} (default) and \code{"lw"}.
The choice \code{"vr"} indicates that the Vietoris-Rips filtration will be used, and the choice \code{"lw"} indicates
that the lazy-witness construction will be used. For Vietoris-Rips filtrations, the parameters \code{landmark_set_size} and
\code{maxmin_samples} are ignored. }
\item{metric}{This indicates the type of metric that will be used. Valid choices include the following: 
\code{"distance_matrix", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"}.
When \code{"distance_matrix"} is specified, the paramter \code{X} is interpreted as a distance matrix. Otherwise,
\code{X} is regarded as a set of points, where each row is one point. }
\item{p}{This is the value of the power to use in the minkowski metric.}
\item{landmark_set_size}{The number of points to include in the landmark set. A sensible value for this is in between 20 and 100. The default value is taken to be \eqn{2 \sqrt{|X|}}. 
This parameter is only relevant for the lazy-witness filtration.}
\item{maxmin_samples}{The number of samples to use when performing the maxmin selection. The default value is taken to be \eqn{\min(|X|, 1000)}.
This parameter is only relevant for the lazy-witness filtration.}
}


