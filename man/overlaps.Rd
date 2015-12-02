\name{Interaction overlaps}
\alias{Interaction-overlaps}

\alias{findOverlaps,GInteractions,GRanges-method}
\alias{findOverlaps,GRanges,GInteractions-method}
\alias{findOverlaps,GInteractions,GRangesList-method}
\alias{findOverlaps,GRangesList,GInteractions-method}
\alias{findOverlaps,GInteractions,GInteractions-method}

\alias{countOverlaps,GInteractions,GRanges-method}
\alias{countOverlaps,GRanges,GInteractions-method}
\alias{countOverlaps,GInteractions,GRangesList-method}
\alias{countOverlaps,GRangesList,GInteractions-method}
\alias{countOverlaps,GInteractions,GInteractions-method}

\alias{overlapsAny,GInteractions,GRanges-method}
\alias{overlapsAny,GRanges,GInteractions-method}
\alias{overlapsAny,GInteractions,GRangesList-method}
\alias{overlapsAny,GRangesList,GInteractions-method}
\alias{overlapsAny,GInteractions,GInteractions-method}

\alias{subsetByOverlaps,GInteractions,GRanges-method}
\alias{subsetByOverlaps,GRanges,GInteractions-method}
\alias{subsetByOverlaps,GInteractions,GRangesList-method}
\alias{subsetByOverlaps,GRangesList,GInteractions-method}
\alias{subsetByOverlaps,GInteractions,GInteractions-method}

\alias{findOverlaps,InteractionSet,GRanges-method}
\alias{findOverlaps,GRanges,InteractionSet-method}
\alias{findOverlaps,InteractionSet,GRangesList-method}
\alias{findOverlaps,GRangesList,InteractionSet-method}
\alias{findOverlaps,InteractionSet,InteractionSet-method}

\alias{countOverlaps,InteractionSet,GRanges-method}
\alias{countOverlaps,GRanges,InteractionSet-method}
\alias{countOverlaps,InteractionSet,GRangesList-method}
\alias{countOverlaps,GRangesList,InteractionSet-method}
\alias{countOverlaps,InteractionSet,InteractionSet-method}

\alias{overlapsAny,InteractionSet,GRanges-method}
\alias{overlapsAny,GRanges,InteractionSet-method}
\alias{overlapsAny,InteractionSet,GRangesList-method}
\alias{overlapsAny,GRangesList,InteractionSet-method}
\alias{overlapsAny,InteractionSet,InteractionSet-method}

\alias{subsetByOverlaps,InteractionSet,GRanges-method}
\alias{subsetByOverlaps,GRanges,InteractionSet-method}
\alias{subsetByOverlaps,InteractionSet,GRangesList-method}
\alias{subsetByOverlaps,GRangesList,InteractionSet-method}
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
\S4method{findOverlaps}{GInteractions,GRanges}(query, subject, maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    select=c("all", "first", "last", "arbitrary"),
    algorithm=c("nclist", "intervaltree"), ignore.strand=TRUE)

\S4method{overlapsAny}{GInteractions,GRanges}(query, subject, maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    algorithm=c("nclist", "intervaltree"), ignore.strand=TRUE)

\S4method{countOverlaps}{GInteractions,GRanges}(query, subject, maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    algorithm=c("nclist", "intervaltree"), ignore.strand=TRUE)

\S4method{subsetByOverlaps}{GInteractions,GRanges}(query, subject, maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    algorithm=c("nclist", "intervaltree"), ignore.strand=TRUE)

# For brevity, only GInteractions,GRanges-methods are listed. Methods are
# available for all pairwise combinations of GInteractions, InteractionSet,
# GRanges and GRangesList objects, so long as at least one InteractionSet or
# GInteractions object is present. In all cases, function calls are identical.
}

\arguments{
\item{query, subject}{
    A GRanges object, a GRangesList object of length 2, a GInteractions object or an InteractionSet object, depending on the specified method.
    At least one of these must be a GInteractions or InteractionSet object.
}

\item{maxgap, minoverlap, type}{
    See \code{?\link[GenomicRanges]{findOverlaps}} in the \pkg{GenomicRanges} package.
}

\item{select, algorithm, ignore.strand}{
    See \code{?\link[GenomicRanges]{findOverlaps}} in the \pkg{GenomicRanges} package.
}
}

\section{Details for GInteractions}{
For all methods taking a GRanges and an GInteractions, 
    an overlap will be defined between a GRanges interval and an GInteractions interaction if either interacting region of the latter overlaps that interval.

If a GRangesList is supplied, it must be of length 2.
Corresponding entries in the constituent GRanges objects are treated as pairs of regions (referred to here as target regions, to distinguish them from the interacting regions).
An overlap is defined between each pair of target regions and an GInteractions interaction if each interacting region overlaps at least one target region,
    and each target region overlaps at least one interacting region, i.e., there are overlapping areas in the two-dimensional interaction space. 

Overlaps between two GInteractions objects are computed in a similar manner.
Interactions from one object are used to define pairs of target regions, which are overlapped with the interactions from the other object in two-dimensional space.

When \code{select="all"}, \code{findOverlaps} returns a Hits object containing overlapping pairs of queries and subjects (or more specifically, their indices in the supplied objects 
    - see \code{?\link[GenomicRanges]{findOverlaps}} for more details).
For other values of \code{select}, an integer vector is returned with one entry for each element of \code{query},
    which specifies the index of the chosen (first, last or arbitrary) overlapping feature in \code{subject} for that query. 
Queries with no overlaps at all are assigned \code{NA} values.

For the other methods, \code{countOverlaps} returns an integer vector indicating the number of elements in \code{subject} that were overlapped by each element in \code{query}.
\code{overlapsAny} returns a logical vector indicating which elements in \code{query} were overlapped by at least one element in \code{subject}. 
\code{subsetByOverlaps} returns a subsetted \code{query} containing only those elements overlapped by at least one element in \code{subject}.
For GRangesList queries, the constituent GRanges objects are subsetted rather than the GRangesList itself.

Note that, unlike the \code{findOverlaps} method for \code{GRanges} objects, strandedness is ignored by default.
This is due to the typical use of GInteractions objects to represent interactions between unstranded genomic regions rather than stranded transcripts.
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
pairing <- GRangesList(first=regions(gi), second=regions(gi))
findOverlaps(pairing, gi)
findOverlaps(gi, pairing)
findOverlaps(gi, pairing, select="last")
overlapsAny(gi, pairing)
overlapsAny(pairing, gi)
countOverlaps(gi, pairing)
countOverlaps(pairing, gi)
subsetByOverlaps(gi, pairing)
subsetByOverlaps(pairing, gi)

# InteractionSet overlaps (split into two):
first.half <- gi[1:(Np/2)]
second.half <- gi[Np/2+1:(Np/2)]
findOverlaps(first.half, second.half)
findOverlaps(first.half, second.half, select="arbitrary")
overlapsAny(first.half, second.half)
countOverlaps(first.half, second.half)
subsetByOverlaps(first.half, second.half)

#################
# Same can be done for an InteractionSet object:

Nlibs <- 4
counts <- matrix(rpois(Nlibs*Np, lambda=10), ncol=Nlibs)
iset <- InteractionSet(counts, gi)

findOverlaps(of.interest, iset)
findOverlaps(iset, pairing)
findOverlaps(iset[1:(Np/2),], iset[Np/2+1:(Np/2),])

# Obviously returns InteractionSet objects instead
subsetByOverlaps(of.interest, iset)
subsetByOverlaps(iset, pairing)
subsetByOverlaps(iset[1:(Np/2),], iset[Np/2+1:(Np/2),])
}

\author{
Aaron Lun
}

\seealso{
\code{\link{InteractionSet-class}},
\code{\link[GenomicRanges]{findOverlaps}},
\code{\link{linkOverlaps}}
}


