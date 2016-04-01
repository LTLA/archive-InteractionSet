\name{Interaction overlaps}
\alias{Interaction-overlaps}

\alias{findOverlaps}
\alias{findOverlaps,GInteractions,Vector-method}
\alias{findOverlaps,Vector,GInteractions-method}
\alias{findOverlaps,GInteractions,GInteractions-method}
\alias{findOverlaps,GInteractions,missing-method}

\alias{countOverlaps}
\alias{countOverlaps,GInteractions,Vector-method}
\alias{countOverlaps,Vector,GInteractions-method}
\alias{countOverlaps,GInteractions,GInteractions-method}
\alias{countOverlaps,GInteractions,missing-method}

\alias{overlapsAny}
\alias{overlapsAny,GInteractions,Vector-method}
\alias{overlapsAny,Vector,GInteractions-method}
\alias{overlapsAny,GInteractions,GInteractions-method}
\alias{overlapsAny,GInteractions,missing-method}

\alias{subsetByOverlaps}
\alias{subsetByOverlaps,GInteractions,Vector-method}
\alias{subsetByOverlaps,Vector,GInteractions-method}
\alias{subsetByOverlaps,GInteractions,GInteractions-method}

\alias{findOverlaps,InteractionSet,Vector-method}
\alias{findOverlaps,Vector,InteractionSet-method}
\alias{findOverlaps,InteractionSet,InteractionSet-method}
\alias{findOverlaps,InteractionSet,missing-method}

\alias{countOverlaps,InteractionSet,Vector-method}
\alias{countOverlaps,Vector,InteractionSet-method}
\alias{countOverlaps,InteractionSet,InteractionSet-method}
\alias{countOverlaps,InteractionSet,missing-method}

\alias{overlapsAny,InteractionSet,Vector-method}
\alias{overlapsAny,Vector,InteractionSet-method}
\alias{overlapsAny,InteractionSet,InteractionSet-method}
\alias{overlapsAny,InteractionSet,missing-method}

\alias{subsetByOverlaps,InteractionSet,Vector-method}
\alias{subsetByOverlaps,Vector,InteractionSet-method}
\alias{subsetByOverlaps,InteractionSet,InteractionSet-method}

\alias{findOverlaps,InteractionSet,GInteractions-method}
\alias{findOverlaps,GInteractions,InteractionSet-method}

\alias{countOverlaps,InteractionSet,GInteractions-method}
\alias{countOverlaps,GInteractions,InteractionSet-method}

\alias{overlapsAny,InteractionSet,GInteractions-method}
\alias{overlapsAny,GInteractions,InteractionSet-method}

\alias{subsetByOverlaps,InteractionSet,GInteractions-method}
\alias{subsetByOverlaps,GInteractions,InteractionSet-method}

\title{Find overlaps between interactions in one or two dimensions}

\description{
Find overlaps between interactions and linear intervals, between interactions and pairs of intervals, 
    and between interactions and other interactions in a GInteractions or InteractionSet object.
}

\usage{
\S4method{findOverlaps}{GInteractions,Vector}(query, subject, maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    select=c("all", "first", "last", "arbitrary"),
    ignore.strand=FALSE, use.region="both")

\S4method{overlapsAny}{GInteractions,Vector}(query, subject, maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    ignore.strand=FALSE, use.region="both")

\S4method{countOverlaps}{GInteractions,Vector}(query, subject, maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    ignore.strand=FALSE, use.region="both")

\S4method{subsetByOverlaps}{GInteractions,Vector}(query, subject, maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    ignore.strand=FALSE, use.region="both")

# For brevity, only 'GInteractions,Vector-methods' are listed. Methods are
# available for all pairwise combinations of GInteractions, InteractionSet, and
# Vector objects, so long as at least one InteractionSet or GInteractions
# object is present. In all cases, function calls are identical.  'subject' can
# also be missing for all functions except for 'subsetByOverlaps', as long as
# 'query' is a GInteractions or an InteractionSet.
}

\arguments{
\item{query, subject}{
    A Vector, GInteractions or InteractionSet object, depending on the specified method.
    At least one of these must be a GInteractions or InteractionSet object.
    Also, \code{subject} can be missing if \code{query} is a GInteractions or InteractionSet object.
}

\item{maxgap, minoverlap, type}{
    See \code{?\link[GenomicRanges]{findOverlaps}} in the \pkg{GenomicRanges} package.
}

\item{select, ignore.strand}{
    See \code{?\link[GenomicRanges]{findOverlaps}} in the \pkg{GenomicRanges} package.
}

\item{use.region}{
A string specifying the regions to be used to identify overlaps.
}
}

\value{
For \code{findOverlaps}, a Hits object is returned if \code{select="all"}, and an integer vector of subject indices otherwise.

For \code{countOverlaps} and \code{overlapsAny}, an integer or logical vector is returned, respectively.

For \code{subsetByOverlaps}, a subsetted object of the same class as \code{query} is returned.
}

\section{Overview of overlaps for GInteractions}{
For all methods taking a Vector and an GInteractions, the Vector is assumed to represent some region on the linear genome (e.g., GRanges) or set of such regions (GRangesList).
An overlap will be defined between the interval and an GInteractions interaction if either anchor region of the latter overlaps the former.
This is considered to be a one-dimensional overlap, i.e., on the linear genome.

For methods between two GInteractions objects, a two-dimensional overlap will be computed between the anchor regions of the two objects.
An overlap is defined if each anchor region of the first object overlaps at least one anchor region of the second object,
    and each anchor region of the second object overlaps at least one anchor region of the first object, i.e., there are overlapping areas in the two-dimensional interaction space. 
If \code{subject} is missing, overlaps will be computed between interactions in \code{query}.
}

\section{Description of overlap methods}{
When \code{select="all"}, \code{findOverlaps} returns a Hits object containing overlapping pairs of queries and subjects (or more specifically, their indices in the supplied objects 
    - see \code{?\link[GenomicRanges]{findOverlaps}} for more details).
For other values of \code{select}, an integer vector is returned with one entry for each element of \code{query},
    which specifies the index of the chosen (first, last or arbitrary) overlapping feature in \code{subject} for that query. 
Queries with no overlaps at all are assigned \code{NA} values.

For the other methods, \code{countOverlaps} returns an integer vector indicating the number of elements in \code{subject} that were overlapped by each element in \code{query}.
\code{overlapsAny} returns a logical vector indicating which elements in \code{query} were overlapped by at least one element in \code{subject}. 
\code{subsetByOverlaps} returns a subsetted \code{query} containing only those elements overlapped by at least one element in \code{subject}.
}

\section{Choice of regions to define overlaps}{
For one-dimensional overlaps, \code{use.region="both"} by default such that overlaps with either interacting region are considered.
If \code{use.region="first"}, overlaps are only considered between the interval and the first interacting region.
Similarly, if \code{use.region="second"}, only the second interaction region is used.

For two-dimensional overlaps, \code{use.region="both"} by default such that the order of first/second interacting/target regions is ignored.
This means that an overlap will be considered between, e.g., the first interacting region and the second target region.
If \code{use.region="same"}, overlaps are only considered between the set of first regions for the GInteractions object and the first regions for the target pairs, and similarly for the second regions.
If \code{use.region="reverse"}, overlaps are only considered between the set of first regions for the GInteractions object and the second regions for the target pairs, and vice versa.
These options tend only to be useful if the order of first/second regions is informative.
}

\section{Details for InteractionSet}{
The behaviour of each method for InteractionSet objects is largely the same as that described for GInteractions objects.
For a given InteractionSet object \code{x}, the corresponding method is called on the GInteractions object in the \code{interactions} slot of \code{x}.
The return value is identical to that from calling the method on \code{interactions(x)}, except for \code{subsetByOverlaps} for InteractionSet queries
    (which returns a subsetted InteractionSet object, containing only those rows/interactions overlapping the \code{subject}).
}

\examples{
example(GInteractions, echo=FALSE)

# Making a larger object, for more overlaps.
Np <- 100
N <- length(regions(gi))
all.anchor1 <- sample(N, Np, replace=TRUE)
all.anchor2 <- sample(N, Np, replace=TRUE)
gi <- GInteractions(all.anchor1, all.anchor2, regions(gi))

# GRanges overlaps:
of.interest <- resize(sample(regions(gi), 2), width=1, fix="center")
findOverlaps(of.interest, gi)
findOverlaps(gi, of.interest)
findOverlaps(gi, of.interest, select="first")
overlapsAny(gi, of.interest)
overlapsAny(of.interest, gi)
countOverlaps(gi, of.interest)
countOverlaps(of.interest, gi)
subsetByOverlaps(gi, of.interest)
subsetByOverlaps(of.interest, gi)

# GRangesList overlaps:
pairing <- GRangesList(first=regions(gi)[1:3], second=regions(gi)[4:6], 
    third=regions(gi)[7:10], fourth=regions(gi)[15:17])
findOverlaps(pairing, gi)
findOverlaps(gi, pairing)
findOverlaps(gi, pairing, select="last")
overlapsAny(gi, pairing)
overlapsAny(pairing, gi)
countOverlaps(gi, pairing)
countOverlaps(pairing, gi)
subsetByOverlaps(gi, pairing)
subsetByOverlaps(pairing, gi)

# GInteractions overlaps (split into two):
first.half <- gi[1:(Np/2)]
second.half <- gi[Np/2+1:(Np/2)]
findOverlaps(first.half, second.half)
findOverlaps(first.half, second.half, select="arbitrary")
overlapsAny(first.half, second.half)
countOverlaps(first.half, second.half)
subsetByOverlaps(first.half, second.half)

findOverlaps(gi)
countOverlaps(gi)
overlapsAny(gi) # trivial result

#################
# Same can be done for an InteractionSet object:

Nlibs <- 4
counts <- matrix(rpois(Nlibs*Np, lambda=10), ncol=Nlibs)
colnames(counts) <- seq_len(Nlibs)
iset <- InteractionSet(counts, gi)

findOverlaps(of.interest, iset)
findOverlaps(iset, pairing)
findOverlaps(iset[1:(Np/2),], iset[Np/2+1:(Np/2),])

# Obviously returns InteractionSet objects instead
subsetByOverlaps(of.interest, iset)
subsetByOverlaps(iset, pairing)
subsetByOverlaps(iset[1:(Np/2),], iset[Np/2+1:(Np/2),])

# Self-overlaps
findOverlaps(iset)
countOverlaps(iset)
overlapsAny(iset) # trivial result
}

\author{
Aaron Lun
}

\seealso{
\code{\link{InteractionSet-class}},
\code{\link[GenomicRanges]{findOverlaps}},
\code{\link{linkOverlaps}}
}


