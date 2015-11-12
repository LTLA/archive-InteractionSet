\name{GInteractions class}
\docType{class}

\alias{GInteractions-class}
\alias{GInteractions}
\alias{GInteractions,numeric,numeric-method}
\alias{GInteractions,GRanges,GRanges-method}
\alias{GInteractions,missing,missing-method}

\title{GInteractions class and constructors}

\description{
The GInteractions class stores pairwise genomic interactions, and is intended for use in data analysis from Hi-C or ChIA-PET experiments.
Each row of the GInteractions corresponds to a pairwise interaction between two loci, with indexing to improve computational efficiency.
}

\usage{
\S4method{GInteractions}{numeric,numeric}(anchor1, anchor2, regions, metadata=list())

\S4method{GInteractions}{GRanges,GRanges}(anchor1, anchor2, regions, metadata=list())

\S4method{GInteractions}{missing,missing}(anchor1, anchor2, regions, metadata=list())
}

\arguments{
\item{anchor1, anchor2}{
    For \code{GInteractions,ANY,numeric,numeric-method}, an integer vector containing indices to \code{regions}.
    
    For \code{GInteractions,ANY,GRanges,GRanges-method}, a GRanges object specifying the interacting loci.
    
    Lengths of both \code{anchor1} and \code{anchor2} must be equal.
}

\item{regions}{
    A GRanges object containing the coordinates of the interacting regions.
    This argument is optional for \code{GInteractions,ANY,GRanges,GRanges-method}.
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
Values of \code{anchor1} and \code{anchor2} are set such that the former is not smaller than the latter, to avoid redundancy from different permutations.

GInteractions objects can be constructed by specifying integer vectors in \code{GInteractions,numeric,numeric-method} to define pairwise interactions.
For entry \code{x}, the corresponding interaction is defined between \code{regions[anchor1[x]]} and \code{regions[anchor2[x]]}.
Obviously, coordinates of all of the interacting loci must be specified in the \code{regions} argument.

Alternatively, GInteractions objects can be constructed by directly supplying the GRanges of the interacting loci to \code{GInteractions,GRanges,GRanges-method}.
If \code{regions} is not specified, this will be constructed automatically from the two sets of supplied GRanges.
If \code{regions} is supplied, exact matching will be performed to identify the corresponding entries in \code{anchor1} and \code{anchor2}.
Missing values are not tolerated and will cause an error.

Both methods will return an GInteractions object containing all of the specified information.
Values may be swapped between the two anchor vectors in the returned object, to ensure that \code{anchor1} is not less than \code{anchor2}.
Sorting of \code{regions} is also performed automatically, with re-indexing of all anchor indices to preserve the correct pairings between regions.
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
elementMetadata(temp.gi)$score <- runif(Np)
}

\author{
Aaron Lun,
based on suggestions from Malcolm Perry and Liz Ing-Simmons.
}

\seealso{
\code{\link{interaction-access}},
\code{\link{interaction-subset}},
\code{\link{interaction-sort}}
}

