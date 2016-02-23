\name{pairs}
\alias{pairs}
\alias{pairs,GInteractions-method}
\alias{pairs,InteractionSet-method}

\title{Extract paired ranges as a GRangesList}

\description{
Represent each interaction in a GInteractions or Interaction object as a GRanges of length 2, in a GRangesList of length equal to the number of interactions.
}

\usage{
\S4method{pairs}{GInteractions}(x, id=FALSE)
\S4method{pairs}{InteractionSet}(x, id=FALSE)
}

\arguments{
\item{x}{A GInteractions or InteractionSet object.}

\item{id}{A logical scalar specifying whether indices should be returned instead of regions.}
}

\value{
A GRangesList if \code{id=FALSE}, and a SelfHits object otherwise.
}

\details{
Recall that the GInteractions object stores anchor regions for each interaction in two parallel GRanges, 
    where corresponding entries between the two GRanges constitute the pair of regions for one interaction.
An alternative representation involves storing the two anchors for each interaction in a single GRanges of length 2.
Multiple interactions are then stored as a GRangesList.

The \code{pairs} method takes a GInteractions object and returns a GRangesList object with the above representation.
This may be more useful in some applications where the two interacting regions must be in the same GRanges.
If \code{id=TRUE}, the anchor indices are extracted and returned as a SelfHits object.
This may be useful for graph construction.
}

\examples{
example(GInteractions, echo=FALSE)
pairs(gi)
pairs(gi, id=TRUE)

example(InteractionSet, echo=FALSE)
pairs(iset)
pairs(iset, id=TRUE)
}

\author{
Aaron Lun
}

\seealso{
\code{\link{GInteractions}},
\code{\link{GRangesList}},
\code{\link{SelfHits}}
}

