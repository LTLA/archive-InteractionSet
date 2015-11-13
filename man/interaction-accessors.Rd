\name{Interaction accessors}
\alias{interaction-access}

\alias{anchors}
\alias{anchors,InteractionSet-method}
\alias{anchors,GInteractions-method}
\alias{anchors<-}
\alias{anchors<-,InteractionSet-method}
\alias{anchors<-,GInteractions-method}

\alias{regions}
\alias{regions,InteractionSet-method}
\alias{regions,GInteractions-method}
\alias{regions<-}
\alias{regions<-,InteractionSet-method}
\alias{regions<-,GInteractions-method}

\alias{newRegions<-}
\alias{newRegions<-,InteractionSet-method}
\alias{newRegions<-,GInteractions-method}

\alias{show}
\alias{show,InteractionSet-method}
\alias{show,GInteractions-method}

\alias{interactions}
\alias{interactions,InteractionSet-method}
\alias{interactions<-}
\alias{interactions<-,InteractionSet-method}

\title{Interaction accessors}

\description{
Methods to get and set fields in an InteractionSet or GInteractions object.
} 

\usage{
\S4method{anchors}{GInteractions}(x, type="both", id=FALSE)
\S4method{anchors}{InteractionSet}(x, type="both", id=FALSE)
\S4method{anchors}{GInteractions}(x) <- value
\S4method{anchors}{InteractionSet}(x) <- value

\S4method{regions}{GInteractions}(x)
\S4method{regions}{InteractionSet}(x)
\S4method{regions}{GInteractions}(x) <- value
\S4method{regions}{InteractionSet}(x) <- value

\S4method{newRegions}{GInteractions}(x) <- value
\S4method{newRegions}{InteractionSet}(x) <- value

\S4method{show}{GInteractions}(x)
\S4method{show}{InteractionSet}(x)

\S4method{interactions}{InteractionSet}(x)
\S4method{interactions}{InteractionSet}(x) <- value
}

\arguments{
\item{x}{An InteractionSet or GInteractions object.}

\item{type}{a string specifying which anchors are to be extracted or replaced.}

\item{id}{a scalar logical indicating whether indices or GRanges should be returned.}

\item{value}{
    For \code{anchors<-}, a list of two integer vectors, each of length equal to the number of rows in \code{x}.

    For \code{regions<-}, a GRanges object of length equal to that of \code{regions(x)}.
    
    For \code{newRegions<-}, a GRanges object that is a superset of all entries in \code{regions(x)} involved in interactions.

    For \code{interactions<-}, a GInteractions object of length equal to that of \code{interactions(x)}.
}
}

\section{Details for GInteractions}{
The return value of \code{anchors} varies depending on \code{type} and \code{id}.
If \code{id=FALSE}, a GRanges object will be returned containing the regions specified by the \code{anchor1} or \code{anchor2} slots in \code{x}, 
    for \code{type=="first"} or \code{"second"}, respectively.
If \code{type="both"}, a GRangesList will be returned with two entries - \code{first} and \code{second}, containing regions specified by \code{anchor1} and \code{anchor2} respectively.
    
However, if \code{id=TRUE}, the integer vectors in the \code{anchor1} or \code{anchor2} slots of \code{x} are returned directly 
    (a list of length two is returned for \code{type="both"}, containing both of the vectors).
Replacement in \code{anchors<-} can only be performed using a list of two integer vectors in \code{value}.
Swapping is performed automatically to ensure that the entries in the \code{anchor1} slot are not less than those in \code{anchor2} in the modified object.
   
For \code{regions}, a GRanges is returned equal to the \code{regions} slot in \code{x}.
For \code{regions<-}, resorting of the input GRanges in \code{value} is performed automatically, along with re-indexing of the anchors.    
In addition, the input GRanges must be of the same length as the existing object in the \code{regions} slot.
The \code{newRegions} replacement method can take variable length GRanges, but requires that the replacement contain (at least) all ranges contained in \code{anchors(x)}.

For access and setting of all other slots, see \code{\link{Vector-class}} for details on the appropriate methods.
This includes \code{mcols} or \code{metadata}, to store interactions-specific or experiment-wide metadata, respectively.

The \code{show} method will print out the class, the number of pairwise interactions, and the total number of regions in the GInteractions object.
The number and names of the various metadata fields will also be printed.
}

\section{Details for InteractionSet}{
Almost all InteractionSet methods operate by calling the equivalent method for the GInteractions object, and returning the resulting value.
The only exception is \code{interactions}, which returns the GInteractions object in the \code{interactions} slot of the InteractionSet object.
This slot can also be set by supplying a valid GInteractions object in \code{interactions<-}.

The \code{show} method will print information equivalent to that done for a SummarizedExperiment object.
An additional line is added indicating the number of regions in the \code{regions} slot of the object.

For access and setting of all other slots, see \code{\link{SummarizedExperiment0-class}} for details on the appropriate methods.
This includes \code{assays}, \code{colData}, \code{mcols} or \code{metadata}, which can all be applied to InteractionSet objects.

Note that getting or setting of \code{mcols} or \code{metadata} will use the fields in the SummarizedExperiment0 class, not in the GInteractions class.
If the latter is desired, these functions should be applied after extracting the GInteractions object with \code{interactions(x)}.
}

\section{Handling different metadata}{
Note that there are several options for metadata - experiment-wide metadata, which goes into \code{metadata(x)<-};
    region-specific metadata (e.g., adjacent genes, promoter/enhancer identity, GC content), which goes into \code{elementMetadata(regions(x))<-};
    and interaction-specific metadata (e.g., interaction distance, interaction type), which goes directly into \code{mcols(x)<-}.
This is applicable to both GInteractions and InteractionSet objects.
}

\examples{
example(GInteractions, echo=FALSE) # Generate a nice object.
show(gi)

# Accessing anchor ranges or indices:
anchors(gi)
anchors(gi, type="first")
anchors(gi, id=TRUE)

anchors(gi, id=TRUE, type="first")
anchors(gi, id=TRUE, type="second")

# Modifying anchor indices:
nregs <- length(regions(gi))
anchors(gi) <- list(sample(nregs, length(gi), replace=TRUE),
                   sample(nregs, length(gi), replace=TRUE))
anchors(gi, id=TRUE, type="first")
anchors(gi, id=TRUE, type="second")

# Alternative with Hits:
hits <- Hits(sample(nregs, length(gi), replace=TRUE),
             sample(nregs, length(gi), replace=TRUE),
             nregs, nregs)
anchors(gi) <- as.data.frame(hits)[,1:2]

# Accessing or modifying regions:
regions(gi)
regions(gi)$score <- runif(length(regions(gi)))

new.ranges <- c(regions(gi), resize(regions(gi), fix="center", width=50))
try(regions(gi) <- new.ranges) # Fails
newRegions(gi) <- new.ranges # Succeeds

#################
# Same can be done for an InteractionSet object:

example(InteractionSet, echo=FALSE)

anchors(iset)
regions(iset)
regions(iset)$score <- regions(gi)$score
newRegions(iset) <- new.ranges

# Standard SE methods also available:
colData(iset)
metadata(iset)
mcols(iset)

# Note the differences in metadata storage:
metadata(iset)$name <- "metadata for SE0"
metadata(interactions(iset))$name <- "metadata for GI"
}

\author{
Aaron Lun
}

\seealso{
\code{\link{GInteractions-class}},
\code{\link{InteractionSet-class}},
\code{\link{Vector-class}}
\code{\link{SummarizedExperiment0-class}}
}
