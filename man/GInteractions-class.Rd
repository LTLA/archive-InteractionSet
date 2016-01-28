\name{GInteractions class}
\docType{class}

\alias{GInteractions-class}
\alias{StrictGInteractions-class}
\alias{ReverseStrictGInteractions-class}
\alias{GInteractions}
\alias{GInteractions,numeric,numeric,GRanges-method}
\alias{GInteractions,GRanges,GRanges,GenomicRangesORmissing-method}
\alias{GInteractions,missing,missing,GenomicRangesORmissing-method}

\title{GInteractions class and constructors}

\description{
The GInteractions class stores pairwise genomic interactions, and is intended for use in data analysis from Hi-C or ChIA-PET experiments.
Each row of the GInteractions corresponds to a pairwise interaction between two loci, with indexing to improve computational efficiency.
}

\usage{
\S4method{GInteractions}{numeric,numeric,GRanges}(anchor1, anchor2, regions, metadata=list(), mode="normal", ...)

\S4method{GInteractions}{GRanges,GRanges,GenomicRangesORmissing}(anchor1, anchor2, regions, metadata=list(), mode="normal", ...)

\S4method{GInteractions}{missing,missing,GenomicRangesORmissing}(anchor1, anchor2, regions, metadata=list(), mode="normal", ...)
}

\arguments{
\item{anchor1, anchor2}{
    Either a pair of numeric vectors containing indices to \code{regions}, or a pair of GRanges objects specifying the interacting loci.
    Lengths of both \code{anchor1} and \code{anchor2} must be equal.
}

\item{regions}{
    A GRanges object containing the coordinates of the interacting regions.
    This is only mandatory if \code{anchor1} and \code{anchor2} are numeric vectors.
} 

\item{metadata}{
    An optional \code{list} of arbitrary content describing the overall experiment.
}

\item{mode}{
    A string indicating what type of GInteractions object should be constructed.
}

\item{...}{
    Optional metadata columns.
}
}

\value{
For the constructors, a GInteractions (or StrictGInteractions, or ReverseStrictGInteractions) object is returned.
}

\section{Description of the class}{
The GInteractions class inherits from the Vector class and has access to all of its data members and methods 
    (e.g, \code{metadata} and \code{elementMetadata} - see \code{\link{Vector-class}} for more details). 
It also contains several additional slots:
\describe{
\item{\code{anchor1}:}{An integer vector specifying the index of the first interacting region.}
\item{\code{anchor2}:}{An integer vector specifying the index of the second interacting region.}
\item{\code{regions}:}{A sorted GRanges object containing the coordinates of all interacting regions.}
}
Each interaction is defined by the corresponding entries in the \code{anchor1} and \code{anchor2} slots, which point to the relevant coordinates in \code{regions} for each locus.

The StrictGInteractions class inherits from the GInteractions class and has the same slots.
The only difference is that, for each interaction, \code{anchor1} must be less than or equal to \code{anchor2}.
This means that the first interacting region has a start position that is "lower" than the second interacting region.
This condition is useful for comparing interactions within and between objects, as it ensures that redundant permutations of the same interaction are not being overlooked.
However, it is not used by default as there may conceivably be instances where the order of interactions is informative.
The ReverseStrictGInteractions class has the opposite behaviour, where \code{anchor1} must be greater than or equal to \code{anchor2}.
}

\section{Class construction}{
GInteractions objects can be constructed by specifying integer vectors to define the pairwise interactions in the \code{GInteractions} call.
For entry \code{x}, the corresponding interaction is defined between \code{regions[anchor1[x]]} and \code{regions[anchor2[x]]}.
Obviously, coordinates of all of the interacting loci must be specified in the \code{regions} argument.
Any metadata in \code{regions} will be preserved.
Note that \code{regions} will be resorted in the returned object, so the anchor indices may not be equal to the input \code{anchor1} and \code{anchor2}.

Alternatively, GInteractions objects can be constructed by directly supplying the GRanges of the interacting loci to the \code{GInteractions} function.
If \code{regions} is not specified, this will be constructed automatically from the two sets of supplied GRanges.
If \code{regions} is supplied, exact matching will be performed to identify the indices in \code{regions} corresponding to the regions in the supplied anchor GRanges.
Missing values are not tolerated and will cause an error to be raised.
In both cases, any metadata in the input GRanges will be transferred to the \code{mcols} of the output GInteractions object.

All constructors will return a GInteractions object containing all of the specified information.
Sorting of \code{regions} is performed automatically, with re-indexing of all anchor indices to preserve the correct pairings between regions.
If \code{mode="strict"}, a StrictGInterctions object is returned with anchor indices swapped such that \code{anchor1 <= anchor2} for all interactions.
If \code{mode="reverse"}, a ReverseStrictGInterctions object is returned with anchor indices swapped such that \code{anchor1 >= anchor2}.
If both anchors are missing, the constructor will return an empty GInteractions object.
}

\examples{
set.seed(1000)
N <- 30
all.starts <- round(runif(N, 1, 100))
all.ends <- all.starts + round(runif(N, 5, 20))
all.regions <- GRanges(rep(c("chrA", "chrB"), c(N-10, 10)),
    IRanges(all.starts, all.ends))

Np <- 20
all.anchor1 <- sample(N, Np)
all.anchor2 <- sample(N, Np)
gi <- GInteractions(all.anchor1, all.anchor2, all.regions)

# Equivalent construction:
gi <- GInteractions(all.regions[all.anchor1], all.regions[all.anchor2])
gi <- GInteractions(all.regions[all.anchor1], all.regions[all.anchor2], 
        all.regions)

# Putting in metadata, elementMetadata
temp.gi <- gi
metadata(temp.gi)$name <- "My GI object"
mcols(temp.gi)$score <- runif(Np)

# Strict construction
sgi <- GInteractions(all.regions[all.anchor1], all.regions[all.anchor2], 
        all.regions, mode="strict")
rsgi <- GInteractions(all.regions[all.anchor1], all.regions[all.anchor2], 
        all.regions, mode="reverse")
}

\author{
Aaron Lun,
with contributions from Malcolm Perry and Liz Ing-Simmons.
}

\seealso{
\code{\link{interaction-access}},
\code{\link{interaction-subset}},
\code{\link{interaction-compare}},
\code{\link{Vector-class}}
}

