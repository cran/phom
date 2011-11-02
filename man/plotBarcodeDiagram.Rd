\name{plotBarcodeDiagram}
\alias{plotBarcodeDiagram}
\title{Create a barcode diagram plot from a set of intervals}
\description{
This function plots a barcode diagram from a given set of intervals. It creates a slightly different
plot from \code{\link{plotPersistenceDiagram}}. As opposed to plotting the starting and
ending points as the horizontal and vertical coordinates of points in the plane, this function displays
a list of barcodes. 
}
\usage{
plotBarcodeDiagram(intervals, dimension, max_f, title="Barcode Diagram")
}
\arguments{
\item{intervals}{A matrix with three columns that specifies the persistence intervals.
Entries in the first column indicate the dimension of an interval. The entries in the
second and third columns indicate the start and end points of the intervals, respectively.
The function \code{\link{pHom}} produces outputs that are in this form.
}
\item{dimension}{The dimension to plot intervals for. Unlike \code{\link{plotPersistenceDiagram}}, this function
only plots intervals for one dimension, and not all of them.}
\item{max_f}{The maximum filtration value to use in the persistence diagram.}
\item{title}{The title on the barcode diagram.}
}

