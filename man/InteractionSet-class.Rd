\name{InteractionSet class}
\docType{class}

\alias{InteractionSet-class}
\alias{InteractionSet}
\alias{InteractionSet,ANY,numeric,numeric-method}
\alias{InteractionSet,ANY,GRanges,GRanges-method}

\title{InteractionSet class and constructors}

\description{
The InteractionSet class stores information about pairwise genomic interactions, and is intended for use in data analysis from Hi-C or ChIA-PET experiments.
Each row of the InteractionSet corresponds to a pairwise interaction between two loci, while each column corresponds to a library or sample.
Each InteractionSet contains one or more assays, designed to hold experimental data about interaction frequencies.
}

\usage{
\S4method{InteractionSet}{ANY,numeric,numeric}(assays, anchor1, anchor2, regions,
    colData=DataFrame(), metadata=list())

\S4method{InteractionSet}{ANY,GRanges,GRanges}(assays, anchor1, anchor2, regions,
    colData=DataFrame(), metadata=list())
}

\arguments{
\item{assays}{
    A numeric matrix or a list or \code{SimpleList} of matrices, containing data for each interaction.
}

\item{anchor1,anchor2}{
    For \code{InteractionSet,ANY,numeric,numeric-method}, an integer vector containing indices to \code{regions}.
    
    For \code{InteractionSet,ANY,GRanges,GRanges-method}, a GRanges object specifying the interacting loci.
}

\item{regions}{
    A GRanges object containing the coordinates of the interacting regions.
    This argument is optional for \code{InteractionSet,ANY,GRanges,GRanges-method}.
} 

\item{colData}{
    An optional \link{DataFrame} describing the samples. 
    Row names, if present, become the column names of the InteractionSet.
}

\item{metadata}{
    An optional \code{list} of arbitrary content describing the overall experiment.
}
}

\details{
The InteractionSet class inherits from the SummarizedExperiment0 class and has access to all of its data members and methods 
    (see \code{\link{SummarizedExperiment0-class}} for more details). 
It also contains several additional slots:
\describe{
\item{\code{anchor1}:}{an integer vector specifying the index of the first interacting region}
\item{\code{anchor2}:}{an integer vector specifying the index of the second interacting region}
\item{\code{regions}:}{a sorted GRanges object containing the coordinates of all interacting regions}
}
Each interaction is defined by the corresponding entries in \code{anchor1} and \code{anchor2}, which point to the relevant coordinates in \code{regions} for each locus.
Values of \code{anchor1} and \code{anchor2} are set such that the former is not smaller than the latter, to avoid redundancy from different permutations.

InteractionSet objects can be constructed by specifying integer vectors in \code{InteractionSet,ANY,numeric,numeric-method} to define pairwise interactions.
For row \code{x}, the corresponding interaction is defined between \code{regions[anchor1[x]]} and \code{regions[anchor2[x]]}.
Obviously, coordinates of all of the interacting loci must be specified in the \code{regions} argument.

Alternatively, InteractionSet objects can be constructed by directly supplying the GRanges of the interacting loci to \code{InteractionSet,ANY,GRanges,GRanges-method}.
If \code{regions} is not specified, this will be constructed automatically from the two sets of supplied GRanges.
If \code{regions} is supplied, exact matching will be performed to identify the corresponding entries in \code{anchor1} and \code{anchor2}.
Missing values are not tolerated and will cause an error.

Both methods will return an InteractionSet object containing all of the specified information.
Values may be swapped between the two anchor vectors in the returned object, to ensure that \code{anchor1} is not less than \code{anchor2}.
Sorting of \code{regions} is also performed automatically, with re-indexing of all anchor indices to preserve the correct pairings between regions.
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
Nlibs <- 4
counts <- matrix(rpois(Np*Nlibs, lambda=10), ncol=Nlibs)
x <- InteractionSet(counts, all.anchor1, all.anchor2, all.regions)

# Equivalent construction:
x <- InteractionSet(counts, all.regions[all.anchor1], 
        all.regions[all.anchor2])
x <- InteractionSet(counts, all.regions[all.anchor1], 
        all.regions[all.anchor2], all.regions)
}

\author{
Aaron Lun
}

\seealso{
\code{\link{InteractionSet-access}},
\code{\link{InteractionSet-subset}},
\code{\link{InteractionSet-sort}}
}

