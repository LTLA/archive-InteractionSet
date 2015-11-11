\name{InteractionSet subsetting}
\alias{InteractionSet-subset}

\alias{[,InteractionSet,ANY-method}
\alias{[,InteractionSet,ANY,ANY-method}
\alias{subset,InteractionSet-method}

\alias{cbind,InteractionSet-method}
\alias{rbind,InteractionSet-method}
\alias{c,InteractionSet-method}
\alias{split,InteractionSet-method}

\title{InteractionSet subsetting and combining}

\description{
Methods to subset or combine InteractionSet objects.
}

\usage{
### Subsetting

\S4method{[}{InteractionSet,ANY,ANY}(x, i, j, ..., drop=TRUE)
\S4method{subset}{InteractionSet}(x, i, j)

### Combining

\S4method{cbind}{InteractionSet}(..., deparse.level=1)
\S4method{rbind}{InteractionSet}(..., deparse.level=1)
\S4method{c}{InteractionSet}(x, ..., recursive=FALSE)

### Splitting

\S4method{split}{InteractionSet}(x, f, drop=FALSE)
}

\arguments{
\item{x}{An InteractionSet object.}

\item{i, j}{A vector of subscripts, indicating the rows and columns to be subsetted for \code{i} and \code{j}, respectively.}

\item{...}{    
    For \code{cbind}, \code{rbind} and \code{c}, \code{...} contains InteractionSet objects to be combined.
    For \code{c}, these are additional to that already specified in \code{x}.
    Otherwise, this argument is ignored during subsetting.
}

\item{f}{
    A vector-like object of length equal to \code{nrow(x)}, indicating how rows of \code{x} are to be split into separate objects.
    This is usually a factor.
}

\item{drop}{
    A logical scalar, indicating if levels of \code{f} that do not occur should be dropped in \code{split,InteractionSet-method}.
    This is ignored by \code{[,InteractionSet,ANY,ANY-method}.
}

\item{deparse.level}{An integer scalar; see \code{?base::\link[base]{cbind}} for a description of this argument.}

\item{recursive}{An integer scalar, ignored.}
}

\details{
Subsetting behaves in much the same way as that for the SummarizedExperiment0 class.
Interactions are treated as rows and will be subsetted as such.
All subsetting operations will return an InteractionSet with the specified interactions (rows) or samples (columns).
However, note that the value of \code{regions} will not be modified by subsetting.

\code{cbind} will combines objects with the same interactions but different samples.
Errors will occur if the interactions are not identical between objects (i.e., must have same values in the slots for \code{regions}, \code{anchor1} and \code{anchor2}).
Additional restrictions apply on the column and assay names - see \code{\link{cbind,SummarizedExperiment0-method}} for details.

\code{rbind} will combines objects with the same samples but different interactions.
In this case, the \code{anchor1} and \code{anchor2} slots need not be identical, but the \code{regions} slot must be the same across objects.
Again, additional restrictions apply - see \code{\link{rbind,SummarizedExperiment0-method}} for details.

\code{c} does the same as \code{rbind}, but does not require that the \code{regions} slot is the same across objects.
Instead, a new GRanges is constructed containing the (sorted) union of all \code{regions} across the input objects.
This is more flexible but is more computationally intensive, hence the distinction from \code{rbind}.

\code{split} produces a list of InteractionSet objects, where each entry of the list corresponds to a level of \code{f}.
Each InteractionSet object contains all rows associated with the corresponding level of \code{f}. 
}

\examples{
example(InteractionSet, echo=FALSE)

# Subsetting:    
x[1,]
x[1:2,]
x[,1]
x[,1:2]
x[3,3]
x[3:4,3:4]

# Combining:
cbind(x, x)
rbind(x, x)

old.x <- x
example(InteractionSet, echo=FALSE)
try(rbind(x, old.x))
c(x, old.x)

# Splitting:
f <- sample(4, nrow(x), replace=TRUE)
out <- split(x, f)
out[[1]]
}

\author{
Aaron Lun
}

\seealso{
\code{\link{InteractionSet-class}}
}


