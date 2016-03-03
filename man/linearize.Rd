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
\S4method{linearize}{GInteractions,numeric}(x, ref, internal=TRUE)

\S4method{linearize}{GInteractions,GRanges}(x, ref, ..., internal=TRUE)

# Equivalent calls for InteractionSet objects.
}

\arguments{
\item{x}{A GInteractions or InteractionSet object}

\item{ref}{ 
    A numeric vector or a GRanges object, specifying the reference region(s) to use for linearization.
    If numeric, the entries should be indices pointing to a genomic interval in \code{regions(x)}.
}

\item{internal}{
    A logical scalar specifying whether interactions within \code{ref} should be reported. 
}

\item{...}{
    Other arguments, passed to \code{\link{overlapsAny}} in the GRanges methods.
}
}

\details{
This method identifies all interactions with at least one anchor region overlapping the specified region(s) in \code{ref}.
When \code{x} is a GInteractions object, the method returns a GRanges object with one entry per identified interaction, 
    where the coordinates are defined as the \emph{other} anchor region, i.e., the one that does \emph{not} overlap the reference region.

If both of the anchor regions for an interaction overlap the reference regions, the genomic interval spanned by both anchor regions is returned.
This is because it is not clear which region should be defined as the "other" anchor in such circumstances.
Note that this will fail if the reference regions occur across multiple chromosomes.
If \code{internal=FALSE}, interactions with both overlapping anchor regions are removed from the output.

When \code{x} is an InteractionSet object, a RangedSummarizedExperiment object is returned.
Each entry corresponds to an identified interaction with the non-overlapping anchor region stored in the \code{rowRanges}.
Experimental data associated with each identified interaction is stored in the various assays.

This method effectively converts two-dimesional interaction data into one-dimensional coverage across the linear genome.
It is useful when a particular genomic region is of interest - this can be used as \code{ref}, to examine the behaviour of all other regions relative to it.
For example, Hi-C data in \code{x} can be converted into pseudo-4C contact frequencies after linearization.

Disjoint ranges across multiple chromosomes are supported when \code{ref} is a GRanges object.
However, it usually only makes sense to use continguous ranges as a single bait region.
Similarly, if \code{ref} is numeric, it should refer to consecutive entries in \code{regions(x)} to specify the bait region.
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


