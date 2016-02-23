\name{Interaction accessors}
\alias{interaction-access}

\alias{anchors}
\alias{anchors,InteractionSet-method}
\alias{anchors,GInteractions-method}
\alias{anchors<-}
\alias{anchors<-,InteractionSet-method}
\alias{anchors<-,GInteractions-method}
\alias{anchors<-,StrictGInteractions-method}
\alias{anchors<-,ReverseStrictGInteractions-method}

\alias{regions}
\alias{regions,InteractionSet-method}
\alias{regions,GInteractions-method}
\alias{regions<-}
\alias{regions<-,InteractionSet-method}
\alias{regions<-,GInteractions-method}

\alias{show,GInteractions-method}
\alias{show,InteractionSet-method}

\alias{replaceRegions<-}
\alias{replaceRegions<-,InteractionSet-method}
\alias{replaceRegions<-,GInteractions-method}

\alias{appendRegions<-}
\alias{appendRegions<-,InteractionSet-method}
\alias{appendRegions<-,GInteractions-method}

\alias{reduceRegions}
\alias{reduceRegions,InteractionSet-method}
\alias{reduceRegions,GInteractions-method}

\alias{seqinfo,InteractionSet-method}
\alias{seqinfo,GInteractions-method}
\alias{seqinfo<-,InteractionSet-method}
\alias{seqinfo<-,GInteractions-method}

\alias{interactions}
\alias{interactions,InteractionSet-method}
\alias{interactions<-}
\alias{interactions<-,InteractionSet-method}

\alias{length,GInteractions-method}
\alias{as.data.frame,GInteractions-method}
\alias{$,GInteractions-method}
\alias{$<-,GInteractions-method}

\alias{mcols,InteractionSet-method}
\alias{mcols<-,InteractionSet-method}

\alias{names,GInteractions-method}
\alias{names<-,GInteractions-method}
\alias{names,InteractionSet-method}
\alias{names<-,InteractionSet-method}

\title{Interaction accessors}

\description{
Methods to get and set fields in an InteractionSet or GInteractions object.
} 

\usage{
\S4method{anchors}{GInteractions}(x, type="both", id=FALSE)
\S4method{regions}{GInteractions}(x)
\S4method{seqinfo}{GInteractions}(x)
\S4method{show}{GInteractions}(x)

\S4method{anchors}{GInteractions}(x, type="both") <- value
\S4method{regions}{GInteractions}(x) <- value
\S4method{seqinfo}{GInteractions}(x) <- value

\S4method{replaceRegions}{GInteractions}(x) <- value
\S4method{appendRegions}{GInteractions}(x) <- value
\S4method{reduceRegions}{GInteractions}(x)

\S4method{names}{GInteractions}(x)
\S4method{names}{GInteractions}(x) <- value

## Each above method has an identical equivalent for 
## InteractionSet objects (not shown here for brevity).

\S4method{anchors}{StrictGInteractions}(x, type="both") <- value
\S4method{anchors}{ReverseStrictGInteractions}(x, type="both") <- value

\S4method{length}{GInteractions}(x)
\S4method{as.data.frame}{GInteractions}(x, row.names=NULL, optional=FALSE, ...)
\S4method{$}{GInteractions}(x, name)
\S4method{$}{GInteractions}(x, name) <- value

\S4method{interactions}{InteractionSet}(x)
\S4method{interactions}{InteractionSet}(x) <- value
\S4method{mcols}{InteractionSet}(x, use.names=FALSE)
\S4method{mcols}{InteractionSet}(x) <- value
}

\arguments{
\item{x}{An InteractionSet or GInteractions object.}

\item{type}{a string specifying which anchors are to be extracted or replaced.}

\item{id}{a scalar logical indicating whether indices or GRanges should be returned.}

\item{name}{a string indicating the field of \code{mcols} to be accessed or modified for a GInteractions object.}

\item{value}{
    For \code{anchors<-} and \code{type="first"} or \code{"second"}, an integer vector of length equal to the number of rows in \code{x}.
    For \code{type="both"}, a list of two such vectors must be supplied.

    For \code{regions<-}, a GRanges object of length equal to that of \code{regions(x)}.
    For \code{replaceRegions<-}, a GRanges object that is a superset of all entries in \code{regions(x)} involved in interactions.
    For \code{appendRegions<-}, a GRanges of any length containing additional regions.

    For \code{seqinfo<-}, a SeqInfo object containing chromosome length data for all regions.
    For \code{interactions<-}, a GInteractions object of length equal to that of \code{interactions(x)}.  
    For \code{mcols<-}, a DataFrame with number of rows equal to the length of \code{x}.
    For \code{names<-}, a character vector of length equal to that of \code{x}.
}

\item{row.names, optional, ...}{
    Additional arguments, see \code{?\link[base]{as.data.frame}} for more details.
}

\item{use.names}{
    A logical scalar, see \code{?\link[S4Vectors]{mcols}} for more details.
}
}

\section{Anchor manipulations for GInteractions}{
The return value of \code{anchors} varies depending on \code{type} and \code{id}.
If \code{id=FALSE}, a GRanges object will be returned containing the regions specified by the \code{anchor1} or \code{anchor2} slots in \code{x}, 
    for \code{type=="first"} or \code{"second"}, respectively.
If \code{type="both"}, a GRangesList will be returned with two entries - \code{first} and \code{second}, containing regions specified by \code{anchor1} and \code{anchor2} respectively.

If \code{id=TRUE}, the integer vectors in the \code{anchor1} or \code{anchor2} slots of \code{x} are returned directly.
A list of length two is returned for \code{type="both"}, containing both of the vectors.
Otherwise, each vector is returned corresponding to the requested value of \code{type}.

Replacement in \code{anchors<-} can only be performed using anchor indices.
If \code{type="both"}, a list of two integer vectors is required in \code{value}, specifying the indices of the interacting regions in \code{regions(x)}.
If \code{type="first"} or \code{"second"}, an integer vector is required to replace the existing values in the \code{anchor1} or \code{anchor2} slot, respectively.
If the object is a StrictGInteractions, indices are automatically swapped so that \code{anchor1 >= anchor2} for each interaction.
The opposite applies if the object is a ReverseStrictGInteractions.
}

\section{Region manipulations for GInteractions}{
For \code{regions}, a GRanges is returned equal to the \code{regions} slot in \code{x}.
For \code{regions<-}, the GRanges in \code{value} is used to replace the \code{regions} slot.
Resorting of the replacement GRanges is performed automatically, along with re-indexing of the anchors.    
In addition, the replacement must be of the same length as the existing object in the \code{regions} slot.

The \code{replaceRegions} replacement method can take variable length GRanges, but requires that the replacement contain (at least) all ranges contained in \code{anchors(x)}.
The \code{appendRegions} replacement method appends extra intervals to the existing \code{regions} slot of \code{x}.
This is more efficient than calling \code{replaceRegions} on a concatenated object with \code{regions(x)} and the extra intervals.

The \code{reduceRegions} method will return a GInteractions object where the \code{regions} slot is reduced to only those entries used in \code{anchors(x)}.
This may save some memory in each object by removing unused regions.
However, this is not recommended for large workflows with many GInteractions objects.
R uses a copy-on-write memory management scheme, so only one copy of the GRanges in \code{regions} should be stored so long as it is not modified in different objects.
}

\section{Other methods for GInteractions}{
For access and setting of all other slots, see \code{\link{Vector-class}} for details on the appropriate methods.
This includes \code{mcols} or \code{metadata}, to store interactions-specific or experiment-wide metadata, respectively.
The \code{length} method will return the number of interactions stored in \code{x}.

The \code{show} method will print out the class, the number of pairwise interactions, and the total number of regions in the GInteractions object.
The number and names of the various metadata fields will also be printed.
The \code{as.data.frame} method will return a data.frame object containing the coordinates for the two anchor regions as well as any element-wise metadata.

The \code{seqinfo} method will return the sequence information of the GRanges in the \code{regions} slot.
This can be replaced with the \code{seqinfo<-} method -- see \code{?\link{Seqinfo}} for more details.
}

\section{Details for InteractionSet}{
Almost all InteractionSet methods operate by calling the equivalent method for the GInteractions object, and returning the resulting value.
The only exception is \code{interactions}, which returns the GInteractions object in the \code{interactions} slot of the InteractionSet object.
This slot can also be set by supplying a valid GInteractions object in \code{interactions<-}.

The \code{show} method will print information equivalent to that done for a SummarizedExperiment object.
An additional line is added indicating the number of regions in the \code{regions} slot of the object.

For access and setting of all other slots, see \code{\link{SummarizedExperiment-class}} for details on the appropriate methods.
This includes \code{assays}, \code{colData}, \code{mcols} or \code{metadata}, which can all be applied to InteractionSet objects.
}

\section{Handling different metadata}{
Note that there are several options for metadata - experiment-wide metadata, which goes into \code{metadata(x)<-};
    region-specific metadata (e.g., adjacent genes, promoter/enhancer identity, GC content), which goes into \code{mcols(regions(x))<-};
    and interaction-specific metadata (e.g., interaction distance, interaction type), which goes directly into \code{mcols(x)<-}.
This is applicable to both GInteractions and InteractionSet objects.
In addition, library-specific data (e.g., library size) should be placed into \code{colData(x)<-} for InteractionSet objects.

Users should take care with the differences in the \code{$} and \code{$<-} operators between these two classes.
For GInteractions objects, this will access or modify fields in the \code{mcols} slot, i.e., for interaction-specific metadata.
For InteractionSet objects, this will access or modify fields in the \code{colData} slot, i.e., for library-specific data.
The difference in behaviour is due to the concept of libraries in the InteractionSet, which is lacking in the GInteractions class.

For InteractionSet objects, the \code{mcols} getter and setter functions operate on the GInteractions object stored in \code{interactions} slot,
    rather than accessing the \code{elementMetadata} slot of the SummarizedExperiment base class.
This makes no difference for practical usage in the vast majority of cases, 
    except that any metadata stored in this manner is carried over when the GInteractions object is extracted with \code{interactions(x)}.
Similarly, the \code{names} getter and setter will operate the names of the internal GInteractions object.
However, the \code{metadata} getter and setter will operate on the SummarizedExperiment base class, \emph{not} on the internal GInteractions object.
}

\value{
For the getters, values in various slots of \code{x} are returned, while for the setters, the slots of \code{x} are modified accordingly -- see below for more details.
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
mod <- list(sample(nregs, length(gi), replace=TRUE),
            sample(nregs, length(gi), replace=TRUE))
anchors(gi) <- mod
anchors(gi, id=TRUE, type="first")
anchors(gi, id=TRUE, type="second")

anchors(gi, type="both") <- mod
anchors(gi, type="first") <- mod[[1]]
anchors(gi, type="first") <- mod[[2]]

# Accessing or modifying regions:
regions(gi)
reduceRegions(gi)
regions(gi)$score <- runif(length(regions(gi)))

new.ranges <- c(regions(gi), resize(regions(gi), fix="center", width=50))
try(regions(gi) <- new.ranges) # Fails
replaceRegions(gi) <- new.ranges # Succeeds

length(regions(gi))
appendRegions(gi) <- GRanges("chrA", IRanges(5:10+1000, 1100+5:10), score=runif(6))
length(regions(gi))

seqinfo(gi)
seqinfo(gi) <- Seqinfo(seqnames=c("chrA", "chrB"), seqlengths=c(1000, 2000))

# Accessing or modifying metadata.
gi$score <- runif(length(gi))
mcols(gi)
as.data.frame(gi)

#################
# Same can be done for an InteractionSet object:

example(InteractionSet, echo=FALSE)

anchors(iset)
regions(iset)
reduceRegions(iset)
regions(iset)$score <- regions(gi)$score
replaceRegions(iset) <- new.ranges

seqinfo(iset)
seqinfo(iset) <- Seqinfo(seqnames=c("chrA", "chrB"), seqlengths=c(1000, 2000))

# Standard SE methods also available:
colData(iset)
metadata(iset)
mcols(iset)

# Note the differences in metadata storage:
metadata(iset)$name <- "metadata for SE0"
metadata(interactions(iset))$name <- "metadata for GI"

iset$lib.size <- runif(ncol(iset))*1e6
colData(iset)
mcols(iset) # untouched by "$" operator

mcols(iset)$whee <- runif(nrow(iset))
mcols(iset)
mcols(interactions(iset)) # preserved

names(iset) <- paste0("X", seq_along(iset))
names(iset)
names(interactions(iset))
}

\author{
Aaron Lun
}

\seealso{
\code{\link{GInteractions-class}},
\code{\link{InteractionSet-class}},
\code{\link{Vector-class}},
\code{\link{SummarizedExperiment-class}}
}
