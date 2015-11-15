\name{Interaction subsetting}
\alias{interaction-subset}

\alias{[,InteractionSet,ANY-method}
\alias{[,InteractionSet,ANY,ANY-method}
\alias{[<-,InteractionSet,ANY,ANY,InteractionSet-method}
\alias{subset,InteractionSet-method}

\alias{[,GInteractions,ANY-method}
\alias{[,GInteractions,ANY,ANY-method}
\alias{[<-,GInteractions,ANY,ANY,GInteractions-method}
\alias{subset,GInteractions-method}

\alias{cbind,InteractionSet-method}
\alias{rbind,InteractionSet-method}
\alias{c,InteractionSet-method}
\alias{combine,InteractionSet,InteractionSet-method}
\alias{split,InteractionSet-method}

\alias{rbind,GInteractions-method}
\alias{c,GInteractions-method}
\alias{combine,GInteractions,GInteractions-method}
\alias{split,GInteractions-method}

\title{Interaction subsetting and combining}

\description{
Methods to subset or combine InteractionSet or GInteractions objects.
}

\usage{
\S4method{[}{GInteractions,ANY,ANY}(x, i, j, ..., drop=TRUE)
\S4method{[}{GInteractions,ANY,ANY,GInteractions}(x, i, j) <- value

\S4method{rbind}{GInteractions}(..., deparse.level=1)
\S4method{c}{GInteractions}(x, ..., recursive=FALSE)
\S4method{combine}{GInteractions,GInteractions}(x, y, ...)

\S4method{split}{GInteractions}(x, f, drop=FALSE)

## Each of the above methods has an identical equivalent for
## InteractionSet objects (not shown for brevity).

\S4method{subset}{InteractionSet}(x, i, j)
\S4method{cbind}{InteractionSet}(..., deparse.level=1)
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

\item{value}{
    A GInteractions or InteractionSet object with length or number of rows equal to length of \code{i} (or that of \code{x}, if \code{i} is not specified).
    For InteractionSet objects, the number of columns must be equal to the length of \code{j} (or number of columns in \code{x}, if \code{j} is not specified).
}

\item{...}{  
    For \code{rbind}, \code{c} and \code{combine}, \code{...} contains GInteractions or InteractionSet objects to be combined row-wise, depending on the method. 
    All objects must be of the same class.
    For all \code{c} and \code{combine} methods, any objects are additional to that already specified in \code{x} (and \code{y}).
    For \code{cbind}, \code{...} contains InteractionSet objects to be combined column-wise.
    For all other methods, this argument is ignored.
}

\item{f}{
    A vector-like object of length equal to \code{nrow(x)}, indicating how rows of \code{x} are to be split into separate objects.
    This is usually a factor.
}

\item{drop}{
    A logical scalar, indicating if levels of \code{f} that do not occur should be dropped in \code{split,InteractionSet-method}.
    This is ignored by the other methods.
}

\item{deparse.level}{An integer scalar; see \code{?base::\link[base]{cbind}} for a description of this argument.}

\item{recursive}{An integer scalar, ignored.}

\item{y}{An object of the same class as \code{x}.}
}

\section{Details for GInteractions}{
Subsetting operations will return a GInteractions object containing the specified interactions.
Values of the \code{anchor1} and \code{anchor2} slots will be appropriately subsetted in the returned object, along with any metadata in \code{mcols}.
However, note that the value of \code{regions} will not be modified by subsetting.

\code{rbind}, \code{c} and \code{combine} will concatenate GInteractions objects.
\code{rbind} and \code{c} are synonyms, and require that the \code{regions} slot of all supplied objects are the same.
\code{combine} does not have this requirement - instead, a new GRanges is constructed containing the (sorted) union of all \code{regions} across the input objects.
This is more flexible but is more computationally intensive.
Note that in both methods, the column names in \code{mcols} must be identical across all supplied objects.

\code{split} produces a list of GInteractions objects, where each entry of the list corresponds to a level of \code{f}.
Each GInteractions object contains all rows associated with the corresponding level of \code{f}. 
}

\section{Details for InteractionSet}{
Subsetting behaves in much the same way as that for the SummarizedExperiment0 class.
Interactions are treated as rows and will be subsetted as such.
All subsetting operations will return an InteractionSet with the specified interactions (rows) or samples (columns).
Again, note that the value of \code{regions} will not be modified by subsetting.

\code{cbind} will combines objects with the same interactions but different samples.
Errors will occur if the interactions are not identical between objects (i.e., must have same values in the \code{interactions} slots).
Additional restrictions apply on the column and assay names - see \code{\link{cbind,SummarizedExperiment0-method}} for details.

\code{rbind} (or equivalently, \code{c}) will combine objects with the same samples but different interactions.
In this case, the interactions need not be identical, but the value of \code{regions(x)} must be the same across objects.
\code{combine} does the same as \code{rbind}, but does not require identical \code{regions}, as described for the GInteractions class.
Again, additional restrictions apply - see \code{\link{rbind,SummarizedExperiment0-method}} for details.

\code{split} produces a list of InteractionSet objects, where each entry of the list corresponds to a level of \code{f}.
Each InteractionSet object contains all rows associated with the corresponding level of \code{f}. 
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

# Combining:
rbind(gi, gi)
c(gi, gi)

new.gi <- gi
regions(new.gi) <- resize(regions(new.gi), width=20, fix="start")
try(rbind(gi, new.gi))
combine(gi, new.gi)

# Splitting:
f <- sample(4, nrow(gi), replace=TRUE)
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

# Combining:
cbind(iset, iset)
rbind(iset, iset)
c(iset, iset)

new.iset <- iset
regions(new.iset) <- resize(regions(new.iset), width=20, fix="start")
try(rbind(iset, new.iset))
combine(iset, new.iset)

# Splitting:
out <- split(iset, f)
out[[1]]
}

\author{
Aaron Lun
}

\seealso{
\code{\link{InteractionSet-class}}
}


