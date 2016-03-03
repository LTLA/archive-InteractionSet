\name{Linearize interactions}

\alias{linearize}
\alias{linearize,GInteractions,GRanges-method}
\alias{linearize,GInteractions,numeric-method}
\alias{linearize,InteractionSet,GRanges-method}
\alias{linearize,InteractionSet,numeric-method}

\title{Linearize 2D interactions into 1D ranges}

\description{
Convert interactions in two-dimensional space to one-dimensional ranges on the linear genome.
}

\usage{
\S4method{linearize}{GInteractions,numeric}(x, ref)

\S4method{linearize}{GInteractions,GRanges}(x, ref, ...)

# Equivalent calls for InteractionSet objects.
}

\arguments{
\item{x}{A GInteractions or InteractionSet object}

\item{ref}{ 
    A numeric scalar or a GRanges object of length 1, specifying the reference region to use for linearization.
    If numeric, this should be an index pointing to an entry in \code{regions(x)}.
}

\item{...}{
    Other arguments, passed to \code{\link{overlapsAny}} in the GRanges methods.
}
}

\details{
This method identifies all interactions with at least one anchor region overlapping the specified region in \code{ref}.
When \code{x} is a GInteractions object, the method returns a GRanges object with one entry per identified interaction, 
    where the coordinates are defined as the \emph{other} anchor region, i.e., the one that does \emph{not} overlap the reference region.
If the anchor regions are identical such that both overlap the reference region, one is chosen and placed into the GRanges.

When \code{x} is an InteractionSet object, a RangedSummarizedExperiment object is returned.
Each entry corresponds to an identified interaction with the non-overlapping anchor region stored in the \code{rowRanges}.
Experimental data associated with each identified interaction is stored in the various assays.

This method effectively converts two-dimesional interaction data into one-dimensional coverage across the linear genome.
It is useful when a particular genomic region is of interest - this can be used as \code{ref}, to examine the behaviour of all other regions relative to it.
For example, Hi-C data in \code{x} can be converted into pseudo-4C contact frequencies after linearization.

If \code{ref} is a GRanges object, it should only overlap a single region in \code{regions(x)}.
A warning will be generated and the first overlapping region in \code{regions(x)} will be used if multiple overlaps are present.
If \code{ref} is a numeric scalar, it is assumed to refer to an entry in \code{regions(x)} to specify the reference region.
}

\value{
A GRanges when \code{x} is a GInteractions object, and a RangedSummarizedExperiment when \code{x} is an InteractionSet object.
}

\examples{
example(InteractionSet, echo=FALSE)

# With integers
out <- linearize(iset, 1)
linearize(iset, 10)
linearize(iset, 20)

# With ranges
linearize(iset, regions(iset)[1], type="equal")
linearize(iset, regions(iset)[10], type="equal")
linearize(iset, regions(iset)[20], type="equal")
}


