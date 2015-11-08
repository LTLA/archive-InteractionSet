\name{InteractionSet sorting}

\alias{InteractionSet-sort}
\alias{order,InteractionSet-method}
\alias{sort,InteractionSet-method}

\title{InteractionSet sorting and ordering}

\description{
Methods to sort and order InteractionSet objects, based on the anchor indices.
}

\usage{
\S4method{order}{InteractionSet}(..., na.last=TRUE, decreasing=FALSE)
\S4method{sort}{InteractionSet}(x, decreasing=FALSE, ...)
}

\arguments{
\item{...}{
    For \code{sort,InteractionSet-method}, one or more InteractionSet objects with the same number of rows.
    Otherwise, ignored in \code{order,InteractionSet-method}.
}

\item{x}{An InteractionSet object.}

\item{na.last}{
    A logical scalar indicating whether \code{NA} values should be ordered last.
    This should not be relevant as anchor indices should be finite.
}

\item{decreasing}{
    A logical scalar indicating whether rows should be sorted in decreasing order.
}
}

\details{
Sorting is performed based on the anchor indices of the InteractionSet object.
Rows are ordered for increasing values of the \code{anchor1} slot - if these are equal, ordering is performed with values of the \code{anchor2} slot.
This equates to ordering by the coordinates directly, as the GRanges in the \code{regions} slot is always sorted.

If multiple objects are supplied to \code{order,InteractionSet-method}, ordering will be computed using anchor indices from successive objects.
Ordering will be performed using anchor indices from the first object; any rows with the same \code{anchor1} and \code{anchor2} will be split using the second object; and so on.
}

\examples{
example(InteractionSet, echo=FALSE)

anchors(x)
x2 <- sort(x)
anchors(x2)

# Can also order them.
o <- order(x)
stopifnot(all.equal(x[o], x2))
}

\author{
Aaron Lun
}

\seealso{
\code{\link{InteractionSet-class}}
}

