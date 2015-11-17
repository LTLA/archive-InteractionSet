\name{Compare interactions}

\alias{match,GInteractions,GInteractions-method}
\alias{match,InteractionSet,GInteractions-method}
\alias{match,GInteractions,InteractionSet-method}
\alias{match,InteractionSet,InteractionSet-method}

\alias{compare,GInteractions,GInteractions-method}

\title{Compare interaction classes}

\description{
Compare or match GInteractions objects with each other, or with InteractionSet objects.
}

\usage{
\S4method{match}{GInteractions,GInteractions}(x, table, nomatch=NA_integer_, incomparables=NULL, ...)

# The same call for 'match' is used for all combinations of GInteractions and
# InteractionSet objects in 'x' or 'table'. For brevity, these will not be
# listed here.

\S4method{compare}{GInteractions,GInteractions}(x, y)
}

\arguments{
\item{x}{
    A GInteractions or InteractionSet object for \code{match}, or a GInteractions object for \code{compare}.
}

\item{table}{A GInteractions or InteractionSet object.}

\item{nomatch, incomparables}{See \code{?\link[base]{match}} for a description of these arguments.}

\item{y}{A GInteractions object, of the same length as \code{x}.}

\item{...}{Additional arguments passed to \code{\link{match,Hits,Hits-method}}.}
}

\details{
The \code{match} function will return an integer vector of length equal to that of \code{x}.
Each entry of the vector corresponds to an interaction in \code{x} and contains the index of the first interaction \code{table} with the same anchor regions.
Interactions in \code{x} without any matches are assigned \code{NA} values by default.
Note that the \code{regions} slot must be identical between \code{x} and \code{table}.

The \code{compare} function will return an integer vector of length equal to \code{x} and \code{y}.
This performs parallel comparisons between corresponding entries in two GInteractions objects, based on the values of the anchor indices 
(\code{anchor1} first, and then \code{anchor2} if \code{anchor1} is tied).
Negative, zero and positive values indicate that the corresponding interaction in \code{x} is 'lesser', equal or 'greater' than the corresponding interaction in \code{y}.
Again, the \code{regions} must be identical between the two objects.
}

\examples{
example(GInteractions, echo=FALSE)

another.gi <- gi[sample(length(gi))]
match(gi, another.gi)
match(gi, another.gi[1:5])

compare(gi, another.gi)
}

\author{
Aaron Lun
}

\seealso{
\code{\link{InteractionSet-class}},
\code{\link{GInteractions-class}},
\code{\link[base]{match}},
\code{\link[S4Vectors]{compare}}
}
