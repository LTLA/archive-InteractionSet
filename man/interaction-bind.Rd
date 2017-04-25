\name{Interaction binding}
\alias{interaction-bind}

\alias{cbind,InteractionSet-method}
\alias{rbind,InteractionSet-method}
\alias{c,GInteractions-method}

\title{Interaction combining}

\description{
Methods to combine InteractionSet or GInteractions objects.
}

\usage{
\S4method{c}{GInteractions}(x, ..., recursive=FALSE)

\S4method{rbind}{InteractionSet}(..., deparse.level=1)
\S4method{cbind}{InteractionSet}(..., deparse.level=1)
}

\arguments{
\item{x}{A GInteractions or InteractionSet object.}

\item{...}{  
    For \code{rbind} and \code{c}, \code{...} contains GInteractions or InteractionSet objects to be combined row-wise.
    All objects must be of the same class.
    For \code{c}, any objects are additional to that already specified in \code{x}.
    For \code{cbind}, \code{...} contains InteractionSet objects to be combined column-wise.
}

\item{deparse.level}{An integer scalar; see \code{?base::\link[base]{cbind}} for a description of this argument.}

\item{recursive}{An integer scalar, ignored.}
}

\value{
A combined object of the same class as \code{x}.
}

\section{Details for GInteractions}{
\code{c} will concatenate GInteractions objects.
It will check whether the \code{regions} slot of all supplied objects are the same, in which case the regions and anchor indices are used directly.
Otherwise, the \code{regions} slot is set to a new GRanges object consisting of the (sorted) union of all \code{regions} across the input objects.
Anchor indices in each object are refactored appropriately to refer to the relevant entries in the new GRanges. 

Note that the column names in \code{mcols} must be identical across all supplied objects.
The column names of \code{mcols} for the \code{regions} slot must also be identical across objects.
If GInteractions objects of different strictness (i.e., StrictGInteractions and ReverseGInteractions) are concatenated,
    the returned object will be of the same class as the first supplied object.
}

\section{Details for InteractionSet}{
\code{cbind} will combine objects with the same interactions but different samples.
Errors will occur if the interactions are not identical between objects (i.e., must have same values in the \code{interactions} slots).
Additional restrictions apply on the column and assay names - see \code{\link{cbind,SummarizedExperiment-method}} for details.

\code{rbind} will combine objects with the same samples but different interactions.
In this case, the interactions need not be identical, and will be concatenated using the methods described above for GInteractions objects.
Again, additional restrictions apply - see \code{\link{rbind,SummarizedExperiment-method}} for details.
}

\examples{
example(GInteractions, echo=FALSE)
c(gi, gi)

new.gi <- gi
regions(new.gi) <- resize(regions(new.gi), width=20, fix="start")
c(gi, new.gi)

#################
# Same can be done for an InteractionSet object:

example(InteractionSet, echo=FALSE)
cbind(iset, iset)
rbind(iset, iset)

new.iset <- iset
regions(new.iset) <- resize(regions(new.iset), width=20, fix="start")
rbind(iset, new.iset)
}

\author{
Aaron Lun
}

\seealso{
\code{\link{InteractionSet-class}}
\code{\link{GInteractions-class}}
}

