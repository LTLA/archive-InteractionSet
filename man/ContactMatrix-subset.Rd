\name{ContactMatrix subsetting}
\alias{ContactMatrix-subset}

\alias{[,ContactMatrix,ANY-method}
\alias{[,ContactMatrix,ANY,ANY-method}
\alias{[,ContactMatrix,ANY,ANY,ANY-method}
\alias{[<-,ContactMatrix,ANY,ANY,ContactMatrix-method}
\alias{subset,ContactMatrix-method}

\alias{cbind,ContactMatrix-method}
\alias{rbind,ContactMatrix-method}

\title{ContactMatrix subsetting and combining}

\description{
Methods to subset or combine ContactMatrix objects.
}

\usage{
### Subsetting

\S4method{[}{ContactMatrix,ANY,ANY}(x, i, j, ..., drop=TRUE)
\S4method{[}{ContactMatrix,ANY,ANY,ContactMatrix}(x, i, j) <- value
\S4method{subset}{ContactMatrix}(x, i, j)

### Combining

\S4method{cbind}{ContactMatrix}(..., deparse.level=1)
\S4method{rbind}{ContactMatrix}(..., deparse.level=1)
}

\arguments{
\item{x}{A ContactMatrix object.}

\item{i, j}{A vector of subscripts, indicating the rows and columns to be subsetted for \code{i} and \code{j}, respectively.}

\item{...}{    
    For \code{cbind}, \code{rbind} and \code{c}, \code{...} contains ContactMatrix objects to be combined.
    Otherwise, this argument is ignored during subsetting.
}

\item{drop}{
    A logical scalar, ignored by \code{[,ContactMatrix,ANY,ANY-method}.
}

\item{value}{
    A ContactMatrix object with dimensions equal to the length of the two subscripts (or the corresponding dimensions of \code{x}, if either subscript is missing.
}

\item{deparse.level}{An integer scalar; see \code{?base::\link[base]{cbind}} for a description of this argument.}
}

\value{
A subsetted or combined ContactMatrix object.
}

\details{
Subsetting of ContactMatrix objects behaves like that for standard matrices.
Rows and columns can be selected and rearranged, with concomitant changes to the \code{anchor1} and \code{anchor2} slots.
All subsetting operations will return an ContactMatrix with the specified rows and columns.
However, note that the value of \code{regions} will not be modified by subsetting.

\code{cbind} will combines objects with the same rows but different columns.
Errors will occur if the row regions are not identical between objects (i.e., must have same values in the slots for \code{regions} and \code{anchor1}).
Conversely, \code{rbind} will combines objects with the same columns but different rows.
Again, errors will occur if the columns are not identical (i.e., must have same values in the slots for \code{regions} and \code{anchor2}).

If subsetted anchors are required, see \code{?"\link{interaction-subset}"} for why subsetting should be done before calling the \code{\link{anchors}} method.
}

\examples{
example(ContactMatrix, echo=FALSE)

# Subsetting:
x[1:5,]
x[,10:15]
x[1:5,10:15]

x2 <- x
x2[1:5,] <- x[6:10,]
as.matrix(x2[,1]) <- 20

# Combining
cbind(x, x)
rbind(x, x)
}

\author{
Aaron Lun
}

\seealso{
\code{\link{ContactMatrix-class}}
}


