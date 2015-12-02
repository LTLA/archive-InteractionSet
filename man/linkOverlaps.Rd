\name{linkOverlaps}
\alias{linkOverlaps}

\alias{linkOverlaps,GInteractions,GRanges,GRanges-method}
\alias{linkOverlaps,InteractionSet,GRanges,GRanges-method}
\alias{linkOverlaps,GInteractions,GRanges,missing-method}
\alias{linkOverlaps,InteractionSet,GRanges,missing-method}

\title{Link overlapping regions}

\description{
Identify interactions that link two sets of regions by having anchor regions overlapping one entry in each set.
}

\usage{
\S4method{linkOverlaps}{GInteractions,GRanges,GRanges}(query, subject1, subject2, ...)
\S4method{linkOverlaps}{GInteractions,GRanges,missing}(query, subject1, subject2, ...)

# The same methods apply for InteractionSet objects as 'query'.
}

\arguments{
\item{query}{A GInteractions or InteractionSet object.}

\item{subject1, subject2}{A GRanges object.}

\item{...}{Additional arguments to be passed to \code{\link{findOverlaps}}.}
}

\details{
This function identifies all interactions in \code{query} where one anchor overlaps an entry in \code{subject1} and the other anchor overlaps an entry in \code{subject2}.
It can be used to identify distinct sets of regions that are linked by interactions in \code{query}.
For example, one might specify genes as \code{subject1} and enhancers as \code{subject2}, to identify all gene-enhancer contacts present in \code{query}.

The function returns a dataframe specifying the index of the interaction in \code{query}; the index of the overlapped region in \code{subject1};
    and the index of the overlapped region in \code{subject2}.
If multiple regions in \code{subject1} and/or \code{subject2} are overlapping the anchor regions of a given interaction,
   all combinations of two regions (one in each \code{subject*} set) are reported for that interaction.

If \code{subject2} is not specified, links within \code{subject1} are identified instead, i.e., \code{subject2} is set to \code{subject1}.
In such cases, the returned dataframe is such that the first subject index is always greater than the second subject index, to avoid redundant permutations.
Note that trival self-interactions are not reported.
}

\seealso{
\code{\link{findOverlaps,GInteractions,GRanges-method}}
}

\author{
Aaron Lun
}

\examples{
example(GInteractions, echo=FALSE)

all.genes <- GRanges("chrA", IRanges(0:9*10, 1:10*10))
all.enhancers <- GRanges("chrB", IRanges(0:9*10, 1:10*10))

out <- linkOverlaps(gi, all.genes, all.enhancers)
head(out)

out <- linkOverlaps(gi, all.genes)
head(out)

# Same methods apply for InteractionSet objects.

example(InteractionSet, echo=FALSE)
out <- linkOverlaps(iset, all.genes, all.enhancers)
out <- linkOverlaps(iset, all.genes)
}
