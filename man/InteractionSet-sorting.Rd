\name{InteractionSet sorting}

\alias{InteractionSet-sort}
\alias{order,InteractionSet-method}
\alias{sort,InteractionSet-method}
\alias{duplicated,InteractionSet-method}
\alias{unique,InteractionSet-method}

\title{InteractionSet sorting and ordering}

\description{
Methods to sort and order InteractionSet objects, based on the anchor indices.
}

\usage{
\S4method{order}{InteractionSet}(..., na.last=TRUE, decreasing=FALSE)
\S4method{sort}{InteractionSet}(x, decreasing=FALSE, ...)
\S4method{duplicated}{InteractionSet}(x, incomparables=FALSE, fromLast=FALSE, ...)
\S4method{unique}{InteractionSet}(x, incomparables=FALSE, fromLast=FALSE, ...)
}

\arguments{
\item{...}{
    For \code{sort,InteractionSet-method}, one or more InteractionSet objects with the same number of rows.
    Otherwise, ignored in all other methods.
}

\item{x}{An InteractionSet object.}

\item{na.last}{
    A logical scalar indicating whether \code{NA} values should be ordered last.
    This should not be relevant as anchor indices should be finite.
}

\item{decreasing}{
    A logical scalar indicating whether rows should be sorted in decreasing order.
}

\item{incomparables}{
    A logical scalar, ignored.
}

\item{fromLast}{
    A logical scalar indicating whether the last entry of a repeated set in \code{x} should be considered as a non-duplicate.
}
}

\details{
Sorting is performed based on the anchor indices of the InteractionSet object.
Rows are ordered for increasing values of the \code{anchor1} slot - if these are equal, ordering is performed with values of the \code{anchor2} slot.
This equates to ordering by the coordinates directly, as the GRanges in the \code{regions} slot is always sorted.
Based on this, \code{sort,InteractionSet-method} will return a sorted copy of \code{x} with permuted rows in increasing order.

\code{order,InteractionSet-method} returns an integer vector indicating the permutation to rearrange \code{x} in increasing order.
If multiple objects are supplied to \code{order}, ordering will be computed using anchor indices from successive objects.
For example, ordering is first performed using anchor indices from the first object; 
    any rows with the same \code{anchor1} and \code{anchor2} will be split using the second object; and so on.

For both of these methods, the sorting can be reversed by setting \code{decreasing=TRUE}.
This will sort or order for decreasing values of \code{anchor1} and \code{anchor2}, rather than for increasing values.

\code{duplicated,InteractionSet-method} returns a logical vector indicating whether the rows of \code{x} are duplicated.
Duplicated rows are identified on the basis of identical entries in the \code{anchor1} and \code{anchor2} slots.
It makes little sense to compare other slots (e.g., in \code{assays}) which may be numeric and inherently incomparable.
For a set of duplicated rows, the first occurrence in \code{x} is marked as the non-duplicate if \code{fromLast=FALSE}, and the last entry otherwise.

\code{unique,InteractionSet-method} returns an InteractionSet object where all duplicate rows have been removed from \code{x}.
This is equivalent to \code{x[!duplicated(x),]}, with any additional arguments passed to \code{duplicated} as specified.
}

\examples{
example(InteractionSet, echo=FALSE)

anchors(x)
x2 <- sort(x)
anchors(x2)

# Can also order them.
o <- order(x)
stopifnot(all.equal(x[o], x2))

# Checking duplication.
summary(duplicated(x))
summary(duplicated(rbind(x, x)))
stopifnot(all.equal(x, unique(rbind(x, x))))
}

\author{
Aaron Lun
}

\seealso{
\code{\link{InteractionSet-class}}
}

