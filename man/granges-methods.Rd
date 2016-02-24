\name{GRanges methods}

\alias{trim,GInteractions-method}
\alias{trim,InteractionSet-method}
\alias{trim,ContactMatrix-method}

\alias{width,GInteractions-method}
\alias{width,InteractionSet-method}
\alias{width,ContactMatrix-method}

\title{Equivalents of GRanges methods}

\description{Methods for GInteractions, InteractionSet and ContactMatrix that operate on the internal GenomicRanges.}

\usage{
\S4method{trim}{GInteractions}(x, use.names=FALSE)
\S4method{width}{GInteractions}(x)

# Equivalent calls for InteractionSet and ContactMatrix.
}

\arguments{
\item{x}{A GInteractions, InteractionSet or ContactMatrix object.}
\item{use.names}{See \code{?"\link{trim,GenomicRanges-method}"}.}
}

\details{
The \code{trim} method will return an object of the same class as \code{x} where out-of-bound ranges in \code{regions(x)} are truncated to fit.
This is followed by resorting of the regions and removal of redundant regions.
As such, that the length of \code{regions(x)} and the anchor indices for each interaction may change before and after trimming.
However, the number and order of the interactions will not change.

For GInteractions and InteractionSet objects, the \code{width} method will return a DataFrame with integer fields \code{anchor1} and \code{anchor2}.
Each row contains the width of the two anchor regions corresponding to each interaction.
For ContactMatrix objects, the method will return a list with the same fields (a different format is necessary as the lengths of \code{anchor1} and \code{anchor2} may be different).
}

\value{
Depending on the method, an object of the same class as \code{x}, or information regarding the genomic regions -- see Details.
}

\seealso{
\code{\link{trim}},
\code{\link{width}}
}

\examples{
example(GInteractions, echo=FALSE)

trim(gi)
width(gi)

example(ContactMatrix, echo=FALSE)
trim(x)
width(x)
}

\author{
Aaron Lun, based on suggestions from Liz Ing-Simmons.
}

