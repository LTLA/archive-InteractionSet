\name{ContactMatrix class}
\docType{class}

\alias{ContactMatrix-class}
\alias{ContactMatrix}
\alias{ContactMatrix,ANY,numeric,numeric-method}
\alias{ContactMatrix,ANY,GRanges,GRanges-method}

\title{ContactMatrix class and constructors}

\description{
The ContactMatrix class contains a matrix where rows and columns represent genomic loci.
Each entry of the matrix contains information about the interaction between the loci represented by the corresponding row/column, e.g., contact frequencies.
Coordinates of the loci are also contained within this class.
}

\usage{
\S4method{ContactMatrix}{ANY,numeric,numeric}(matrix, anchor1, anchor2, regions)

\S4method{ContactMatrix}{ANY,GRanges,GRanges}(matrix, anchor1, anchor2, regions)
}

\arguments{
\item{matrix}{
    A matrix containing interaction data.
}

\item{anchor1,anchor2}{
    For \code{ContactMatrix,ANY,numeric,numeric-method}, an integer vector containing indices to \code{regions}.
    
    For \code{ContactMatrix,ANY,GRanges,GRanges-method}, a GRanges object specifying the interacting loci.
    
    In both cases, \code{anchor1} and \code{anchor2} should have lengths equal to the number of rows and columns in \code{matrix}, respectively.
}

\item{regions}{
    A GRanges object containing the coordinates of the interacting regions.
    This argument is optional for \code{InteractionSet,ANY,GRanges,GRanges-method}.
} 
}

\details{
The ContactMatrix class contains several slots:
\describe{
\item{\code{matrix}:}{A matrix.}
\item{\code{anchor1}:}{An integer vector specifying the index of the first interacting region.}
\item{\code{anchor2}:}{An integer vector specifying the index of the second interacting region.}
\item{\code{regions}:}{A sorted GRanges object containing the coordinates of all interacting regions.}
}
Each entry of \code{anchor1} corresponds to a row in \code{matrix}, while each entry of \code{anchor2} corresponds to a column.
Each entry of \code{matrix} represents an interaction between the corresponding entries in \code{anchor1} and \code{anchor2}, i
    which point to the relevant coordinates in \code{regions} for each locus.

ContactMatrix objects can be constructed by specifying integer vectors in \code{ContactMatrix,ANY,numeric,numeric-method} to define pairwise interactions.
Coordinates of all of the interacting loci must be specified in the \code{regions} argument.
This means that the range of entries must lie within \code{[1, length(regions)]}.:w

Alternatively, ContactMatrix objects can be constructed by directly supplying the GRanges of the interacting loci to \code{ContactMatrix,ANY,GRanges,GRanges-method}.
If \code{regions} is not specified, this will be constructed automatically from the two sets of supplied GRanges.
If \code{regions} is supplied, exact matching will be performed to identify the indices corresponding to \code{anchor1} and \code{anchor2}.
Missing values are not tolerated and will cause an error.

Both methods will return an ContactMatrix object containing all of the specified information.
Sorting of \code{regions} is also performed automatically, with re-indexing of all anchor indices to preserve the correct pairings between regions.
}

\examples{
set.seed(1000)
N <- 30
all.starts <- round(runif(N, 1, 100))
all.ends <- all.starts + round(runif(N, 5, 20))
all.regions <- GRanges(rep(c("chrA", "chrB"), c(N-10, 10)),
    IRanges(all.starts, all.ends))

Nr <- 10
Nc <- 20
all.anchor1 <- sample(N, Nr)
all.anchor2 <- sample(N, Nc)
counts <- matrix(rpois(Nr*Nc, lambda=10), Nr, Nc)
x <- ContactMatrix(counts, all.anchor1, all.anchor2, all.regions)

# Equivalent construction:
x <- ContactMatrix(counts, all.regions[all.anchor1], 
        all.regions[all.anchor2])
x <- ContactMatrix(counts, all.regions[all.anchor1], 
        all.regions[all.anchor2], all.regions)
}

\author{
Aaron Lun
}

\seealso{
\code{\link{ContactMatrix-access}},
\code{\link{ContactMatrix-subset}},
\code{\link{ContactMatrix-sort}}
}

