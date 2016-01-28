\name{ContactMatrix distances}
\alias{ContactMatrix-dist}

\alias{pairdist,ContactMatrix-method}
\alias{intrachr,ContactMatrix-method}

\title{Compute linear distances from ContactMatrix objects}

\description{
Methods to compute linear distances between pairs of interacting regions in a ContactMatrix object.
}

\usage{
\S4method{pairdist}{ContactMatrix}(x, type="mid")

\S4method{intrachr}{ContactMatrix}(x)
}

\arguments{
\item{x}{A ContactMatrix object.}

\item{type}{
    A character string specifying the type of distance to compute.
    See \code{?\link{pairdist,InteractionSet-method}} for possible values.
}
}

\value{
An integer or logical matrix of the same dimensions as \code{x}, containing the specified distances.
}

\details{
\code{pairdist,ContactMatrix-method} will return a matrix of integer (or, if \code{type="intra"}, logical) values.
Each entry of this matrix specifies the distance between the interacting loci that are represented by the corresponding row and column.
If \code{type="intra"}, each entry indicates whether the corresponding interaction is intra-chromosomal.
Running \code{intrachr(x)} is equivalent to \code{pairdist(x, type="intra")} for any ContactMatrix object \code{x}.
See \code{\link{pairdist,InteractionSet-method}} for more details on the type of distances that can be computed.
}

\examples{
example(ContactMatrix, echo=FALSE)

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
\code{\link{ContactMatrix-class}},
\code{\link{pairdist,InteractionSet-method}}
}

