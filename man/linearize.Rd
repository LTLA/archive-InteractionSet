\name{Linearize interactions}

\alias{linearize}
\alias{linearize,InteractionSet,GRanges-method}
\alias{linearize,InteractionSet,numeric-method}

\title{Linearize 2D interactions into 1D ranges}

\description{
Convert interactions in two-dimensional space to one-dimensional ranges on the linear genome.
}

\usage{
\S4method{linearize}{InteractionSet,numeric}(x, ref, ...)

\S4method{linearize}{InteractionSet,GRanges}(x, ref, ...)
}

\arguments{
\item{x}{An InteractionSet object}

\item{ref}{ 
    A numeric scalar or a GRanges object of length 1, specifying the reference region to use for linearization.
    If numeric, this should be an index pointing to an entry in \code{regions(x)}.
}

\item{...}{
    Other arguments, passed to \code{\link[GenomicRanges]{overlapsAny}} in \code{linearize,InteractionSet,GRanges-method}.
}
}

\details{
This method identifies all interactions with at least one interacting locus overlapping the specified region.
It returns a RangedSummarizedExperiment object with one entry per interaction, 
    where the corresponding \code{rowRanges} is defined as the interacting locus for each interaction that does not overlap the reference region.
If the interacting loci are identical such that both overlap the reference region, one is chosen and placed into \code{rowRanges}.

This effectively converts two-dimesional interaction data into one-dimensional coverage across the linear genome.
It is useful when a particular genomic region is of interest - this can be used as \code{ref}, to examine the behaviour of all other regions relative to it.
For example, Hi-C data in \code{x} can be converted into pseudo-4C contact frequencies after linearization.

If \code{ref} is a GRanges object, it should only overlap a single region in \code{regions(x)}.
A warning will be generated and the first region will be used if multiple overlaps are present.
If \code{ref} is a numeric scalar, it is assumed to refer to an entry in \code{regions(x)} to specify the reference region.
}

\examples{
example(InteractionSet, echo=FALSE)

# With integers
out <- linearize(x, 1)
linearize(x, 10)
linearize(x, 20)

# With ranges
linearize(x, regions(x)[1], type="equal")
linearize(x, regions(x)[10], type="equal")
linearize(x, regions(x)[20], type="equal")
}


