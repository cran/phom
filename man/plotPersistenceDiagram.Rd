\name{plotPersistenceDiagram}
\alias{plotPersistenceDiagram}
\title{Create a persistence diagram plot from a set of intervals}
\description{
This function plots a persistence diagram from a given set of intervals.
A persistence diagram is a set of unordered points \eqn{\{(x_i, y_i) | y_i \geq x_i \} } for
each dimension, where \eqn{x_i} indicates the starting point of interval \eqn{i}
and \eqn{y_i} indicates the ending point of interval \eqn{i}.
This function plots these points as coordinates in the plane, using difference colors and
markers for each dimension.

We warn the user that this function may result in two identical points occluding each other. 
In other words, the multiplicity of each interval is not clearly visible. The function
\code{\link{plotBarcodeDiagram}} does not suffer from this.
}
\usage{
plotPersistenceDiagram(intervals, max_dim, max_f, title="Persistence Diagram")
}
\arguments{
\item{intervals}{A matrix with three columns that specifies the persistence intervals.
Entries in the first column indicate the dimension of an interval. The entries in the
second and third columns indicate the start and end points of the intervals, respectively.
The function \code{\link{pHom}} produces outputs that are in this form.
}
\item{max_dim}{The maximum dimension to plot. Note that unlike \code{\link{plotBarcodeDiagram}},
this function plots points (intervals) for all dimensions up to the maximum dimension.}
\item{max_f}{The maximum filtration value to use in the persistence diagram.}
\item{title}{The title on the persistence diagram.}
}


