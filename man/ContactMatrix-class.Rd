\name{ContactMatrix class}
\docType{class}

\alias{ContactMatrix-class}
\alias{ContactMatrix}
\alias{ContactMatrix,ANY,numeric,numeric,GRanges-method}
\alias{ContactMatrix,ANY,GRanges,GRanges,GenomicRangesORmissing-method}
\alias{ContactMatrix,missing,missing,missing,GenomicRangesORmissing-method}

\title{ContactMatrix class and constructors}

\description{
The ContactMatrix class contains a matrix where rows and columns represent genomic loci.
Each entry of the matrix contains information about the interaction between the loci represented by the corresponding row/column, e.g., contact frequencies.
Coordinates of the loci are also contained within this class.
}

\usage{
\S4method{ContactMatrix}{ANY,numeric,numeric,GRanges}(matrix, anchor1, anchor2, regions, metadata=list())

\S4method{ContactMatrix}{ANY,GRanges,GRanges,GenomicRangesORmissing}(matrix, anchor1, anchor2, regions, metadata=list())

\S4method{ContactMatrix}{missing,missing,missing,GenomicRangesORmissing}(matrix, anchor1, anchor2, regions, metadata=list())
}

\arguments{
\item{matrix}{
    Any object containing interaction data that can be coerced into a Matrix object.
}

\item{anchor1, anchor2}{
    Either a pair of numeric vectors containing indices to \code{regions} or a pair of GRanges objects.
    In both cases, \code{anchor1} and \code{anchor2} should have lengths equal to the number of rows and columns in \code{matrix}, respectively.
}

\item{regions}{
    A GRanges object containing the coordinates of the interacting regions.
    This argument is optional for \code{InteractionSet,ANY,GRanges,GRanges-method}.
} 

\item{metadata}{
    A list containing experiment-wide metadata -  see \code{?\link[S4Vectors]{Annotated}} for more details.
}
}

\value{
For the constructors, a ContactMatrix object is returned.
}

\details{
The ContactMatrix class inherits from the \code{\link[S4Vectors]{Annotated}} class, with several additional slots:
\describe{
\item{\code{matrix}:}{A Matrix object.}
\item{\code{anchor1}:}{An integer vector specifying the index of the first interacting region.}
\item{\code{anchor2}:}{An integer vector specifying the index of the second interacting region.}
\item{\code{regions}:}{A sorted GRanges object containing the coordinates of all interacting regions.}
}
Each entry of \code{anchor1} corresponds to a row in \code{matrix}, while each entry of \code{anchor2} corresponds to a column.
Each entry of \code{matrix} represents an interaction between the corresponding entries in \code{anchor1} and \code{anchor2}, i
    which point to the relevant coordinates in \code{regions} for each locus.

ContactMatrix objects can be constructed by specifying numeric vectors as \code{anchor1} and \code{anchor2} in the \code{ContactMatrix} function.
These vectors will define the regions corresponding to the rows and columns of the matrix.
Specifically, each value of the vector acts as an index to specify the relevant coordinates from \code{regions}.
This means that the range of entries must lie within \code{[1, length(regions)]}.

Alternatively, ContactMatrix objects can be constructed by directly supplying the GRanges of the interacting loci in \code{ContactMatrix}.
If \code{regions} is not specified, this will be constructed automatically from the two sets of supplied GRanges.
If \code{regions} is supplied, exact matching will be performed to identify the indices in \code{regions} corresponding to the regions in the supplied GRanges.
Missing values are not tolerated and will cause an error to be raised.

Both methods will return an ContactMatrix object containing all of the specified information.
Sorting of \code{regions} is also performed automatically, with re-indexing of all anchor indices to preserve the correct pairings between regions.
}

\section{Choosing between Matrix classes}{
The ContactMatrix class provides support for Matrix objects from the Matrix package.
This provides several classes for storing sparse and symmetric matrices in a memory-efficient manner, e.g., via the sparseMatrix class.
These classes may be preferable for handling very large matrices spanning substantial areas in the interaction space.

There are, however, some caveats to the use of the various Matrix subclasses.
All numeric values are coerced to double-precision types, which may take up more memory than a direct integer representation.
Character data is not supported, though this is unlikely to be relevant in most cases.
Another issue is how missing values should be interpreted in the sparseMatrix -- see \code{?\link{inflate}} for more details.

In general, the choice of class depends on the type of data and the intended application.
A dense Matrix is often sufficient for most use cases where small or dense regions of the interaction space are being examined.
Switching to sparseMatrix may be desirable in emptier parts of the space, where there can be substantial memory savings.
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

# Also works directly with Matrix objects.
counts2 <- Matrix::Matrix(counts)
x2 <- ContactMatrix(counts2, all.anchor1, all.anchor2, all.regions)
counts2 <- as(counts2, "dgCMatrix")
x2 <- ContactMatrix(counts2, all.anchor1, all.anchor2, all.regions)
}

\author{
Aaron Lun
}

\seealso{
\code{\link{ContactMatrix-access}},
\code{\link{ContactMatrix-subset}},
\code{\link{ContactMatrix-sort}},
\code{\link{Annotated-class}},
\code{\link{Matrix-class}}
}

