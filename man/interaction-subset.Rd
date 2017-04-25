\name{Interaction subsetting}
\alias{interaction-subset}

\alias{[,InteractionSet,ANY-method}
\alias{[,InteractionSet,ANY,ANY-method}
\alias{[,InteractionSet,ANY,ANY,ANY-method}
\alias{[<-,InteractionSet,ANY,ANY,InteractionSet-method}
\alias{subset,InteractionSet-method}

\alias{[<-,GInteractions,ANY,GInteractions-method}

\title{Interaction subsetting and combining}

\description{
Methods to subset or combine InteractionSet or GInteractions objects.
}

\usage{
\S4method{[}{InteractionSet,ANY,ANY}(x, i, j, ..., drop=TRUE)

\S4method{[}{InteractionSet,ANY,ANY,InteractionSet}(x, i, j) <- value

\S4method{subset}{InteractionSet}(x, i, j)

\S4method{[}{GInteractions,ANY,GInteractions}(x, i) <- value
}

\arguments{
\item{x}{A GInteractions or InteractionSet object.}

\item{i, j}{
    A vector of logical or integer subscripts.
    For InteractionSet objects, these indicate the rows and columns to be subsetted for \code{i} and \code{j}, respectively.
    Rows correspond to pairwise interactions while columns correspond to samples.
    For GInteractions objects, \code{i} indicates the genomic interactions to be retained.
    \code{j} is ignored as there is no concept of samples in this class.
}

\item{..., drop}{Additional arguments that are ignored.}

\item{value}{
    A GInteractions or InteractionSet object with length or number of rows equal to length of \code{i} (or that of \code{x}, if \code{i} is not specified).
    For InteractionSet objects, the number of columns must be equal to the length of \code{j} (or number of columns in \code{x}, if \code{j} is not specified).
}
}

\value{
A subsetted object of the same class as \code{x}.
}

\section{Details for GInteractions}{
Subsetting operations are not explicitly listed above as they inherit from the \code{\link{Vector}} class.
They will return a GInteractions object containing the specified interactions.
Values of the \code{anchor1} and \code{anchor2} slots will be appropriately subsetted in the returned object, along with any metadata in \code{mcols}.
However, note that the value of \code{regions} will not be modified by subsetting.

For short index vectors, subsetting a GInteractions object prior to calling \code{\link{anchors}} may be much faster than the reverse procedure.
This is because the \code{\link{anchors}} getter will construct a GRanges(List) containing the genomic loci for all pairwise interactions.
Subsetting beforehand ensures that only loci for the desired interactions are included.
This avoids constructing the entire object just to subset it later.

Subset assignment will check if the \code{regions} are identical between \code{x} and \code{value}.
If not, the \code{regions} slot in the output object will be set to a sorted union of all \code{regions} from \code{x} and \code{value}.
Indices are refactored appropriately to point to the entries in the new \code{regions}.
}

\section{Details for InteractionSet}{
Subsetting behaves in much the same way as that for the SummarizedExperiment class.
Interactions are treated as rows and will be subsetted as such.
All subsetting operations will return an InteractionSet with the specified interactions (rows) or samples (columns).
Again, note that the value of \code{regions} will not be modified by subsetting.
}

\examples{
example(GInteractions, echo=FALSE)

# Subsetting:    
gi[1,]
gi[1:2,]
gi[3]
gi[3:4]

temp.gi <- gi
temp.gi[3:4] <- gi[1:2]

# Splitting:
f <- sample(4, length(gi), replace=TRUE)
out <- split(gi, f)
out[[1]]

#################
# Same can be done for an InteractionSet object:

example(InteractionSet, echo=FALSE)

# Subsetting:    
iset[1,]
iset[1:2,]
iset[,1]
iset[,1:2]
iset[3,3]
iset[3:4,3:4]

# Splitting:
out <- split(iset, f)
out[[1]]
}

\author{
Aaron Lun
}

\seealso{
\code{\link{InteractionSet-class}}
\code{\link{GInteractions-class}}
}


