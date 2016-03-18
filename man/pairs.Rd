\name{pairs}
\alias{pairs}
\alias{pairs,GInteractions-method}
\alias{pairs,InteractionSet-method}
\alias{makeGInteractionsFromGRangesPairs}

\title{Extract paired ranges}

\description{
Represent interactions in a GInteractions or Interaction object as a Pairs, SelfHits or GRangesList object.
}

\usage{
\S4method{pairs}{GInteractions}(x, id=FALSE, as.grlist=FALSE)

# Equivalent call for InteractionSet.

makeGInteractionsFromGRangesPairs(x)
}

\arguments{
\item{x}{
For \code{pairs}, a GInteractions or InteractionSet object 
For \code{makeGInteractionsFromGRangesPairs}, a Pairs object containing two parallel GRanges.
}

\item{id}{A logical scalar specifying whether indices should be returned instead of regions.}

\item{as.grlist}{A logical scalar specifying whether a GRangesList should be returned.}
}

\value{
For \code{pairs}, if \code{id=TRUE}, a SelfHits object is returned. 
Otherwise, if \code{as.grlist=TRUE}, a GRangesList object is returned.
Otherwise, a Pairs object is returned.

For \code{makeGInteractionsFromGRangesPairs}, a GInteractions object is returned.
}

\details{
Recall that the GInteractions object stores anchor regions for each interaction in two parallel GRanges, 
    where corresponding entries between the two GRanges constitute the pair of regions for one interaction.
These parallel ranges can be extracted and stored as a Pairs object for further manipulation.
This is similar to the GRangesList reported by \code{\link{anchors}} with \code{type="both"} and \code{id=FALSE}.
The reverse conversion is performed using \code{makeGInteractionsFromGRangesPairs}.

An alternative representation involves storing the two anchors for each interaction in a single GRanges of length 2.
Multiple interactions are then stored as a GRangesList, along with any metadata and sequence information.
This is returned if \code{as.grlist=FALSE}, may be more useful in some applications where the two interacting regions must be in the same GRanges.
Finally, if \code{id=TRUE}, the anchor indices are extracted and returned as a SelfHits object.
This may be useful for graph construction.
}

\examples{
example(GInteractions, echo=FALSE)
y <- pairs(gi)
y
makeGInteractionsFromGRangesPairs(y)

pairs(gi, id=TRUE)
pairs(gi, as.grlist=TRUE)

example(InteractionSet, echo=FALSE)
pairs(iset)
pairs(iset, id=TRUE)
pairs(iset, as.grlist=TRUE)
}

\author{
Aaron Lun
}

\seealso{
\code{\link{GInteractions}},
\code{\link{Pairs}},
\code{\link{SelfHits}},
\code{\link{GRangesList}}
}

