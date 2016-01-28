\name{ContactMatrix sorting}

\alias{ContactMatrix-sort}
\alias{order,ContactMatrix-method}
\alias{sort,ContactMatrix-method}
\alias{duplicated,ContactMatrix-method}
\alias{unique,ContactMatrix-method}

\title{ContactMatrix sorting and ordering}

\description{
Methods to sort and order ContactMatrix objects, based on the anchor indices.
Also, methods to remove duplicate rows or columns in each ContactMatrix.
}

\usage{
\S4method{order}{ContactMatrix}(..., na.last=TRUE, decreasing=FALSE)
\S4method{sort}{ContactMatrix}(x, decreasing=FALSE, ...)
\S4method{duplicated}{ContactMatrix}(x, incomparables=FALSE, fromLast=FALSE, ...)
\S4method{unique}{ContactMatrix}(x, incomparables=FALSE, fromLast=FALSE, ...)
}

\arguments{
\item{...}{
    For \code{sort,ContactMatrix-method}, one or more ContactMatrix objects with the same dimensions.
    Otherwise, ignored in all other methods.
}

\item{x}{A ContactMatrix object.}

\item{na.last}{
    A logical scalar indicating whether \code{NA} values should be ordered last.
    This should not be relevant as anchor indices should be finite.
}

\item{decreasing}{
    A logical scalar indicating whether sorting should be performed in decreasing order.
}

\item{incomparables}{
    A logical scalar, ignored.
}

\item{fromLast}{
    A logical scalar indicating whether the last entry of a repeated set in \code{x} should be considered as a non-duplicate.
}
}

\value{
For \code{sort} and \code{unique}, a ContactMatrix object is returned with sorted or unique rows/columns, respectively.

For \code{order}, a list of two integer vectors is returned containing row/column permutations.

For \code{duplicated}, a list of logical vectors is returned specifying which rows/columns are duplicated.
}

\details{
Sorting is performed based on the anchor indices of the ContactMatrix object.
Rows are ordered for increasing values of the \code{anchor1} slot, while columns are ordered for increasing values of the \code{anchor2} slot.
This equates to ordering by the coordinates directly, as the GRanges in the \code{regions} slot is always sorted.
Based on this, \code{sort,ContactMatrix-method} will return a sorted copy of \code{x} with permuted rows/columns in increasing order.
This can be set to decreasing order by specifying \code{decreasing=TRUE}.

\code{order,ContactMatrix-method} returns a list of 2 integer vectors.
The first vector contains the permutation to rearrange the rows of \code{x} in increasing order, while the second vector does the same for the columns of \code{x}.
If multiple objects are supplied to \code{order}, ordering will be computed using anchor indices from successive objects.
In other words, ordering will be performed using anchor indices from the first object; 
    any rows with the same \code{anchor1} or columns with the same \code{anchor2} will be split using the corresponding indices in the second object; and so on.

\code{duplicated,ContactMatrix-method} returns a list of two logical vectors.
The first vector indicates whether rows are duplicated, based on identical values in the \code{anchor1} slot.
The second does the same for columns based on the \code{anchor2} slot.
For a set of duplicated rows or columns, the first occurrence of that row/column is marked as the non-duplicate if \code{fromLast=FALSE}, and the last entry otherwise.

\code{unique,ContactMatrix-method} returns an ContactMatrix object where all duplicate rows and columns have been removed from \code{x}.
This is equivalent to subsetting based on the non-duplicate rows and columns defined in \code{duplicated,ContactMatrix-method}.

Note that sorting and duplicate identification only use the anchor indices.
The values of the \code{matrix} slot are \emph{not} used in distinguishing between rows or columns with the same index.
}

\examples{
example(ContactMatrix, echo=FALSE)

anchors(x)
x2 <- sort(x)
x2
anchors(x2)

# Can also order them.
o <- order(x)
o
stopifnot(all.equal(x[o$row,o$column], x2))

# Checking duplication.
duplicated(x)
duplicated(rbind(x, x))
stopifnot(all.equal(x, unique(rbind(x, x))))
}

\author{
Aaron Lun
}

\seealso{
\code{\link{ContactMatrix-class}}
}

