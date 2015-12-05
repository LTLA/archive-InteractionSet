\name{flip}
\alias{flip}
\alias{flip,GInteractions-method}
\alias{flip,GRangesList-method}

\title{Flip the arrangement of anchor regions}

\description{
Flip a GInteractions object so that each interaction is a GRanges of length 2, in a GRangesList.
The opposing conversion is also provided.
}

\usage{
\S4method{flip}{GInteractions}(x)
\S4method{flip}{GRangesList}(x, ...)
}

\arguments{
\item{x}{A GInteractions or GRangesList object.}

\item{...}{Other arguments to be passed to \code{\link{GInteractions}}.}
}

\details{
The GInteractions object stores anchor regions for each interaction in two parallel GRanges, 
    where corresponding entries between the two GRanges constitute the pair of regions for one interaction.
An alternative representation involves storing the two anchors for each interaction in a single GRanges of length 2.
Multiple interactions are then stored as a GRangesList.

The \code{flip} method takes a GInteractions object and returns a GRangesList object with the above representation.
This may be more useful in some applications where the two interacting regions must be in the same GRanges.
The converse can also be performed, whereby a GRangesList is converted into an appropriate GInteractions object.
In both cases, metadata and names will be transferred from one object to another.
}

\examples{
example(GInteractions, echo=FALSE)

out <- flip(gi)
out

flip(out) # Automated 'regions' construction.

flip(out, regions=regions(gi))
}

\author{
Aaron Lun
}

\seealso{
\code{\link{GInteractions}},
\code{\link{GRangesList}}
}

