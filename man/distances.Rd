\name{Interaction distances}

\alias{pairdist}
\alias{pairdist,InteractionSet-method}
\alias{intrachr}
\alias{intrachr,InteractionSet-method}

\title{Get the linear distance for each interaction}

\description{
Compute the distance between interacting regions on the linear genome, for each pairwise interaction contained in an InteractionSet object.
}

\usage{
\S4method{pairdist}{InteractionSet}(x, type="mid")

\S4method{intrachr}{InteractionSet}(x)
}

\arguments{
\item{x}{An InteractionSet object.}

\item{type}{
    A character string specifying the type of distance to compute.
    Can take values of \code{"mid"}, \code{"gap"}, \code{"span"}, \code{"diag"} or \code{"intra"}.
}
}

\details{
For each interaction in \code{x}, the \code{pairdist} method computes the distance between the two interacting regions.
An integer vector is returned, with values computed according to the specified value of \code{type}:
\describe{
\item{\code{"mid"}:}{
    The distance between the midpoints of the two regions (rounded down to the nearest integer) is returned.
}

\item{\code{"gap"}:}{
    The length of the gap between the closest points of the two regions is computed - 
    negative lengths are returned for overlapping regions, indicating the length of the overlap.
}

\item{\code{"span"}:}{
    The distance between the furthermost points of the two regions is computed.
}

\item{\code{"diag"}:}{
    The difference between the anchor indices is returned.
    This corresponds to a diagonal on the interaction space when bins are used in the \code{regions} slot of \code{x}.
}
}

Interchromosomal interactions are marked with \code{NA}.
Alternatively, if \code{type="intra"}, a logical vector is returned indicating whether the interaction occurs between two regions on the same chromosome.
\code{intrachr(x)} is an alias for \code{pairdist(x, type="intra")}.
}

\examples{
example(InteractionSet, echo=FALSE)

pairdist(x)
pairdist(x, type="gap")
pairdist(x, type="span")
pairdist(x, type="diag")

intrachr(x)
}

\author{
Aaron Lun
}

\seealso{
\code{\link{InteractionSet-class}}
}

