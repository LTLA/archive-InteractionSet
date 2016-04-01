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
\S4method{linkOverlaps}{GInteractions,GRanges,GRanges}(query, subject1, subject2, ..., use.region="both")
\S4method{linkOverlaps}{GInteractions,GRanges,missing}(query, subject1, subject2, ..., use.region="both")

# The same methods apply for InteractionSet objects as 'query'.
}

\arguments{
\item{query}{A GInteractions or InteractionSet object.}

\item{subject1, subject2}{A GRanges object.}

\item{...}{Additional arguments to be passed to \code{\link{findOverlaps}}.}

\item{use.region}{
A string specifying which \code{query} regions should be used to overlap which \code{subject}.
Ignored if \code{subject2} is missing.
}
}

\details{
This function identifies all interactions in \code{query} where one anchor overlaps an entry in \code{subject1} and the other anchor overlaps an entry in \code{subject2}.
It is designed to be used to identify regions that are linked by interactions in \code{query}.
For example, one might specify genes as \code{subject1} and enhancers as \code{subject2}, to identify all gene-enhancer contacts present in \code{query}.
This is useful when the exact pairings between \code{subject1} and \code{subject2} are undefined.

The function returns a dataframe specifying the index of the interaction in \code{query}; the index of the overlapped region in \code{subject1};
    and the index of the overlapped region in \code{subject2}.
If multiple regions in \code{subject1} and/or \code{subject2} are overlapping the anchor regions of a particular interaction,
    all combinations of two overlapping regions (one from each \code{subject*} set) are reported for that interaction.

By default, \code{use.region="both"} such that overlaps will be considered between any first/second interacting region in \code{query} and either \code{subject1} or \code{subject2}.
If \code{use.region="same"}, overlaps will only be considered between the first interacting region in \code{query} and entries in \code{subject1}, 
    and between the second interacting region and \code{subject2}.
The opposite applies with \code{use.region="reverse"}, where the first and second interacting regions are overlapped with \code{subject2} and \code{subject1} respectively.

If \code{subject2} is not specified, links within \code{subject1} are identified instead, i.e., \code{subject2} is set to \code{subject1}.
In such cases, the returned dataframe is such that the first subject index is always greater than the second subject index, to avoid redundant permutations.
}

\value{
A dataframe of integer indices indicating which elements of \code{query} link which elements of \code{subject1} and \code{subject2}.
}

\seealso{
\code{\link{findOverlaps,GInteractions,Vector-method}}
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
