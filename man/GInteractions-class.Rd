\name{GInteractions class}
\docType{class}

\alias{GInteractions-class}
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
\S4method{GInteractions}{numeric,numeric,GRanges}(anchor1, anchor2, regions, metadata=list())

\S4method{GInteractions}{GRanges,GRanges,GenomicRangesORmissing}(anchor1, anchor2, regions, metadata=list())

\S4method{GInteractions}{missing,missing,GenomicRangesORmissing}(anchor1, anchor2, regions, metadata=list())
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
}

\details{
The GInteractions class inherits from the Vector class and has access to all of its data members and methods 
    (e.g, \code{metadata} and \code{elementMetadata} - see \code{\link{Vector-class}} for more details). 
It also contains several additional slots:
\describe{
\item{\code{anchor1}:}{An integer vector specifying the index of the first interacting region.}
\item{\code{anchor2}:}{An integer vector specifying the index of the second interacting region.}
\item{\code{regions}:}{A sorted GRanges object containing the coordinates of all interacting regions.}
}
Each interaction is defined by the corresponding entries in the \code{anchor1} and \code{anchor2} slots, which point to the relevant coordinates in \code{regions} for each locus.

GInteractions objects can be constructed by specifying integer vectors to define the pairwise interactions in the \code{GInteractions} call.
For entry \code{x}, the corresponding interaction is defined between \code{regions[anchor1[x]]} and \code{regions[anchor2[x]]}.
Obviously, coordinates of all of the interacting loci must be specified in the \code{regions} argument.
Any metadata in \code{regions} will be preserved.

Alternatively, GInteractions objects can be constructed by directly supplying the GRanges of the interacting loci to the \code{GInteractions} function.
If \code{regions} is not specified, this will be constructed automatically from the two sets of supplied GRanges.
Both GRanges should have the same set of metadata fields, which will be retained in the constructed \code{regoins}.
If \code{regions} is supplied, exact matching will be performed to identify the indices in \code{regions} corresponding to the regions in the supplied anchor GRanges.
Any metadata in the anchor GRanges will be ignored in favour of metadata in the supplied \code{regions}.
Missing values are not tolerated and will cause an error to be raised.

All constructors will return an GInteractions object containing all of the specified information.
Sorting of \code{regions} is performed automatically, with re-indexing of all anchor indices to preserve the correct pairings between regions.
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
}

\author{
Aaron Lun,
based on suggestions from Malcolm Perry and Liz Ing-Simmons.
}

\seealso{
\code{\link{interaction-access}},
\code{\link{interaction-subset}},
\code{\link{interaction-compare}}
}

