\name{Convert classes}

\alias{inflate}
\alias{inflate,InteractionSet-method}
\alias{deflate}
\alias{deflate,ContactMatrix-method}

\title{Convert between classes}

\description{
Inflate an InteractionSet into a ContactMatrix, or deflate a ContactMatrix to an InteractionSet.
}

\usage{
\S4method{inflate}{InteractionSet}(x, rows, columns, assay=1L, sample=1L, fill=NULL, ...)

\S4method{deflate}{ContactMatrix}(x, unique=TRUE, ...)
}

\arguments{
\item{x}{
    An InteractionSet object for \code{inflate}, or a ContactMatrix object for \code{deflate}.
}

\item{rows, columns}{
    An integer vector, character vector or GRanges object indicating the regions of interest to be used as the rows or columns of the ContactMatrix.

    If it is an integer vector, it is assumed to refer to intervals in the \code{regions} slot of \code{x}.
    Values of the vector need not be sorted or unique, but must lie within \code{[1, regions(x)]}.

    If it is a character vector, it is assumed to contain the names of the reference sequences of interest (i.e., chromosome names).

    If it is a GRanges object, \code{\link[GenomicRanges]{overlapsAny}} will be called to identify the overlapping intervals of the \code{regions} slot.
}

\item{assay}{
    A numeric scalar indicating the assay of the InteractionSet object, from which values are extracted to fill the ContactMatrix.
}

\item{sample}{
    A numeric scalar indicating the sample (i.e., column) of the assay to extract values to fill the ContactMatrix.
}

\item{fill}{
    A numeric vector of length equal to \code{nrow(x)}, with alternative values to fill the ContactMatrix.
    If specified, this overrides extraction of assay values.
}

\item{unique}{
    A  logical scalar indicating whether duplicated interactions should be removed from \code{x} prior to deflation.
}

\item{...}{
    For \code{inflate,InteractionSet-method}, additional arguments to pass to \code{\link[GenomicRanges]{overlapsAny}} when \code{rows} or \code{columns} is a GRanges.

    For \code{deflate,ContactMatrix-method}, additional arguments to pass to the InteractionSet constructor.
}
}

\details{
The \code{inflate} method will return a ContactMatrix where the rows and columns correspond to specified regions of interest in \code{rows} and \code{columns}.
Values in the matrix are filled in based on the assay values in the InteractionSet (or \code{fill}, if specified).
Some entries of the matrix will correspond to pairwise interactions that are not present in the InteractionSet - these are filled with \code{NA} values.
If duplicated interactions are present in \code{x}, only the first interaction will be used.

The \code{deflate} method will return an InteractionSet where each non-\code{NA} entry in the ContactMatrix is converted into a pairwise interaction.
The values of all matrix entries are stored as a one-sample assay.
Duplicate interactions are not included by default, but can be returned by setting \code{unique=FALSE}.
Additional arguments can be used to specify the \code{colData} and \code{metadata}, which are stored in the ContactMatrix itself.
}

\examples{
example(InteractionSet, echo=FALSE)

inflate(x, 1:10, 1:10)
inflate(x, 1:10, 1:5+10)
inflate(x, "chrA", 1:5+10)
inflate(x, "chrA", "chrB")
inflate(x, "chrA", GRanges("chrB", IRanges(1, 10)))

y <- inflate(x, 1:10, 1:10)
x2 <- deflate(y)
x2
assay(x2)

# Testing with different fillings:
y <- inflate(x, 1:10, 1:10, sample=2)
x2 <- deflate(y)
assay(x2)

y <- inflate(x, 1:10, 1:10, fill=rowSums(assay(x)))
x2 <- deflate(y)
assay(x2)

# Effect of 'unique'
y <- inflate(x, c(8, 1:10), 1:10)
deflate(y)
deflate(y, unique=FALSE)
}

\author{
Aaron Lun
}

\seealso{
\code{\link{InteractionSet-class}},
\code{\link{ContactMatrix-class}}
}
