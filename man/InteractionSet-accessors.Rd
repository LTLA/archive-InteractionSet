\name{InteractionSet accessors}
\alias{InteractionSet-access}

\alias{anchors}
\alias{anchors,InteractionSet-method}
\alias{anchors<-}
\alias{anchors<-,InteractionSet-method}

\alias{regions}
\alias{regions,InteractionSet-method}
\alias{regions<-}
\alias{regions<-,InteractionSet-method}

\alias{newRegions<-}
\alias{newRegions<-,InteractionSet-method}

\alias{show}
\alias{show,InteractionSet-method}

\title{InteractionSet accessors}

\description{
Methods to get and set fields in an InteractionSet object.
} 

\usage{
\S4method{anchors}{InteractionSet}(x, type="both", id=FALSE)
\S4method{anchors}{InteractionSet}(x) <- value

\S4method{regions}{InteractionSet}(x)
\S4method{regions}{InteractionSet}(x) <- value
\S4method{newRegions}{InteractionSet}(x) <- value

\S4method{show}{InteractionSet}(x)
}

\arguments{
\item{x}{An InteractionSet object.}

\item{type}{a string specifying which anchors are to be extracted or replaced.}

\item{id}{a scalar logical indicating whether indices or GRanges should be returned.}

\item{value}{
    For \code{anchors<-,InteractionSet-method}, a list of two integer vectors, each of length equal to the number of rows in \code{x}.

    For \code{regions<-,InteractionSet-method}, a GRanges object of length equal to that of \code{regions(x)}.
    
    For \code{newRegions<-,InteractionSet-method}, a GRanges object that is a superset of all entries in \code{regions(x)} involved in interactions.
}
}

\details{
The return value of \code{anchors,InteractionSet-method} varies depending on \code{type} and \code{id}.
If \code{id=FALSE}, a GRanges object will be returned containing the regions specified by the \code{anchor1} or \code{anchor2} slots in \code{x}, 
    for \code{type=="first"} or \code{"second"}, respectively.
If \code{type="both"}, a GRangesList will be returned with two entries - \code{first} and \code{second}, containing regions specified by \code{anchor1} and \code{anchor2} respectively.
    
However, if \code{id=TRUE}, the integer vectors in the \code{anchor1} or \code{anchor2} slots of \code{x} are returned directly 
    (a list of length two is returned for \code{type="both"}, containing both of the vectors).
Replacement in \code{anchors<-,InteractionSet-method} can only be performed using a list of two integer vectors in \code{value}.
Swapping is performed automatically to ensure that the entries in the \code{anchor1} slot are not less than those in \code{anchor2} in the modified object.
   
For \code{regions,InteractionSet-method}, a GRanges is returned equal to the \code{regions} slot in \code{x}.
For \code{regions<-,InteractionSet-method}, resorting of the input GRanges in \code{value} is performed automatically, along with re-indexing of the anchors.    
In addition, the input GRanges must be of the same length as the existing object in the \code{regions} slot.
The \code{newRegions} replacement method can take variable length GRanges, but requires that the replacement contain (at least) all ranges contained in \code{anchors(x)}.

The \code{show} method will print information equivalent to that done for a SummarizedExperiment object.
An additional line is added indicating the number of regions in the \code{regions} slot of the object.

For access and setting of all other slots, see \code{\link{SummarizedExperiment0-class}} for details on the appropriate methods.
This includes \code{assays}, \code{colData}, \code{mcols} or \code{metadata}, which can all be applied to InteractionSet objects.

Note that there are several options for metadata - experiment-wide metadata, which goes into \code{metadata(x)<-};
    region-specific metadata (e.g., adjacent genes, promoter/enhancer identity, GC content), which goes into \code{elementMetadat(regions(x))<-};
    and interaction-specific metadata (e.g., interaction distance, interaction type), which goes directly into \code{mcols(x)<-}.
}

\examples{
example(InteractionSet, echo=FALSE) # Generate a nice object.
show(x)

# Accessing anchor ranges or indices:
anchors(x)
anchors(x, type="first")
anchors(x, id=TRUE)

anchors(x, id=TRUE, type="first")
anchors(x, id=TRUE, type="second")

# Modifying anchor indices:
nregs <- length(regions(x))
anchors(x) <- list(sample(nregs, nrow(x), replace=TRUE),
                   sample(nregs, nrow(x), replace=TRUE))
anchors(x, id=TRUE, type="first")
anchors(x, id=TRUE, type="second")

# Alternative with Hits:
hits <- Hits(sample(nregs, nrow(x), replace=TRUE),
             sample(nregs, nrow(x), replace=TRUE),
             nregs, nregs)
anchors(x) <- as.data.frame(hits)[,1:2]

# Accessing or modifying regions:
regions(x)
regions(x)$score <- runif(length(regions(x)))

new.ranges <- c(regions(x), resize(regions(x), fix="center", width=50))
try(regions(x) <- new.ranges) # Fails
newRegions(x) <- new.ranges # Succeeds

# Standard SE methods:
colData(x)
metadata(x)
mcols(x)
}

\author{
Aaron Lun
}

\seealso{
\code{\link{InteractionSet-class}},
\code{\link{SummarizedExperiment0-class}}
}
