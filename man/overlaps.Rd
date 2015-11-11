\name{Interaction overlaps}

\alias{findOverlaps}
\alias{findOverlaps,InteractionSet,GRanges-method}
\alias{findOverlaps,GRanges,InteractionSet-method}
\alias{findOverlaps,InteractionSet,GRangesList-method}
\alias{findOverlaps,GRangesList,InteractionSet-method}
\alias{findOverlaps,InteractionSet,InteractionSet-method}

\alias{countOverlaps}
\alias{countOverlaps,InteractionSet,GRanges-method}
\alias{countOverlaps,GRanges,InteractionSet-method}
\alias{countOverlaps,InteractionSet,GRangesList-method}
\alias{countOverlaps,GRangesList,InteractionSet-method}
\alias{countOverlaps,InteractionSet,InteractionSet-method}

\alias{overlapsAny}
\alias{overlapsAny,InteractionSet,GRanges-method}
\alias{overlapsAny,GRanges,InteractionSet-method}
\alias{overlapsAny,InteractionSet,GRangesList-method}
\alias{overlapsAny,GRangesList,InteractionSet-method}
\alias{overlapsAny,InteractionSet,InteractionSet-method}

\alias{subsetByOverlaps}
\alias{subsetByOverlaps,InteractionSet,GRanges-method}
\alias{subsetByOverlaps,GRanges,InteractionSet-method}
\alias{subsetByOverlaps,InteractionSet,GRangesList-method}
\alias{subsetByOverlaps,GRangesList,InteractionSet-method}
\alias{subsetByOverlaps,InteractionSet,InteractionSet-method}

\title{Find overlaps between interactions in one or two dimensions}

\description{
Find overlaps between interactions and linear intervals, between interactions and pairs of intervals, and between interactions and other interactions.
}

\usage{
\S4method{findOverlaps}{InteractionSet,GRanges}(query, subject, maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    select=c("all", "first", "last", "arbitrary"),
    algorithm=c("nclist", "intervaltree"), 
    ignore.strand=TRUE)
\S4method{findOverlaps}{GRanges,InteractionSet}(query, subject, maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    select=c("all", "first", "last", "arbitrary"),
    algorithm=c("nclist", "intervaltree"),
    ignore.strand=TRUE)
\S4method{findOverlaps}{InteractionSet,GRangesList}(query, subject, maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    select=c("all", "first", "last", "arbitrary"),
    algorithm=c("nclist", "intervaltree"),
    ignore.strand=TRUE)
\S4method{findOverlaps}{GRangesList,InteractionSet}(query, subject, maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    select=c("all", "first", "last", "arbitrary"),
    algorithm=c("nclist", "intervaltree"),
    ignore.strand=TRUE)
\S4method{findOverlaps}{InteractionSet,InteractionSet}(query, subject, maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    select=c("all", "first", "last", "arbitrary"),
    algorithm=c("nclist", "intervaltree"),
    ignore.strand=TRUE)

\S4method{countOverlaps}{InteractionSet,GRanges}(query, subject, maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    algorithm=c("nclist", "intervaltree"),
    ignore.strand=TRUE)
\S4method{countOverlaps}{GRanges,InteractionSet}(query, subject, maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    algorithm=c("nclist", "intervaltree"),
    ignore.strand=TRUE)
\S4method{countOverlaps}{InteractionSet,GRangesList}(query, subject, maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    algorithm=c("nclist", "intervaltree"),
    ignore.strand=TRUE)
\S4method{countOverlaps}{GRangesList,InteractionSet}(query, subject, maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    algorithm=c("nclist", "intervaltree"),
    ignore.strand=TRUE)
\S4method{countOverlaps}{InteractionSet,InteractionSet}(query, subject, maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    algorithm=c("nclist", "intervaltree"),
    ignore.strand=TRUE)


\S4method{overlapsAny}{InteractionSet,GRanges}(query, subject, maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    algorithm=c("nclist", "intervaltree"),
    ignore.strand=TRUE)
\S4method{overlapsAny}{GRanges,InteractionSet}(query, subject, maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    algorithm=c("nclist", "intervaltree"),
    ignore.strand=TRUE)
\S4method{overlapsAny}{InteractionSet,GRangesList}(query, subject, maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    algorithm=c("nclist", "intervaltree"),
    ignore.strand=TRUE)
\S4method{overlapsAny}{GRangesList,InteractionSet}(query, subject, maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    algorithm=c("nclist", "intervaltree"),
    ignore.strand=TRUE)
\S4method{overlapsAny}{InteractionSet,InteractionSet}(query, subject, maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    algorithm=c("nclist", "intervaltree"),
    ignore.strand=TRUE)

\S4method{subsetByOverlaps}{InteractionSet,GRanges}(query, subject, maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    algorithm=c("nclist", "intervaltree"),
    ignore.strand=TRUE)
\S4method{subsetByOverlaps}{GRanges,InteractionSet}(query, subject, maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    algorithm=c("nclist", "intervaltree"),
    ignore.strand=TRUE)
\S4method{subsetByOverlaps}{InteractionSet,GRangesList}(query, subject, maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    algorithm=c("nclist", "intervaltree"),
    ignore.strand=TRUE)
\S4method{subsetByOverlaps}{GRangesList,InteractionSet}(query, subject, maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    algorithm=c("nclist", "intervaltree"),
    ignore.strand=TRUE)
\S4method{subsetByOverlaps}{InteractionSet,InteractionSet}(query, subject, maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    algorithm=c("nclist", "intervaltree"),
    ignore.strand=TRUE)
}

\arguments{
\item{query, subject}{
    A GRanges object, a GRangesList object of length 2, or an InteractionSet object, depending on the specified method.
    One of these must be an InteractionSet object.
}

\item{maxgap, minoverlap, type}{
    See \code{?\link[GenomicRanges]{findOverlaps}} in the \pkg{GenomicRanges} package.
}

\item{select, algorithm, ignore.strand}{
    See \code{?\link[GenomicRanges]{findOverlaps}} in the \pkg{GenomicRanges} package.
}
}

\details{
For all methods taking a GRanges and an InteractionSet, 
    an overlap will be defined between a GRanges interval and an InteractionSet interaction if either interacting region of the latter overlaps that interval.

If a GRangesList is supplied, it must be of length 2.
Corresponding entries in the constituent GRanges objects are treated as pairs of regions (referred to here as target regions, to distinguish them from the interacting regions).
An overlap is defined between each pair of target regions and an InteractionSet interaction if each interacting region overlaps at least one target region,
    and each target region overlaps at least one interacting region, i.e., there are overlapping areas in the two-dimensional interaction space. 

Overlaps between two InteractionSet objects are computed in a similar manner.
Interactions from one object are used to define pairs of target regions, which are overlapped with the interactions from the other object as described for the GRangesList.

When \code{select="all"}, \code{findOverlaps} returns a Hits object containing overlapping pairs of queries and subjects (or more specifically, their indices in the supplied objects - 
    see \code{?\link[GenomicRanges]{findOverlaps}} for more details).
For other values of \code{select}, an integer vector is returned that specifies the index of the chosen (first, last or arbitrary) overlapping feature in \code{subject} for each element of \code{query}.
Queries with no overlaps at all are assigned \code{NA} values.

For the other methods, \code{countOverlaps} returns an integer vector indicating the number of elements in \code{subject} that were overlapped by each element in \code{query}.
\code{overlapsAny} returns a logical vector indicating which elements in \code{query} were overlapped by at least one element in \code{subject}. 
\code{subsetByOverlaps} returns a subsetted \code{query} containing only those elements overlapped by at least one element in \code{subject}.
For GRangesList queries, the constituent GRanges objects are subsetted rather than the GRangesList itself.

Note that, unlike \code{findOverlaps} for \code{GRanges} methods, strandedness is ignored by default.
This is due to the typical use of InteractionSet objects for interactions between unstranded genomic regions rather than stranded transcripts.
}

\examples{
example(InteractionSet, echo=FALSE)

# Making a larger object, for more overlaps.
Np <- 100
all.anchor1 <- sample(N, Np, replace=TRUE)
all.anchor2 <- sample(N, Np, replace=TRUE)
Nlibs <- 4
counts <- matrix(rpois(Np*Nlibs, lambda=10), ncol=Nlibs)
x <- InteractionSet(counts, all.anchor1, all.anchor2, all.regions)

# GRanges overlaps:
of.interest <- resize(sample(regions(x), 2), width=1, fix="center")
findOverlaps(of.interest, x)
findOverlaps(x, of.interest)
findOverlaps(x, of.interest, select="first")
overlapsAny(x, of.interest)
overlapsAny(of.interest, x)
countOverlaps(x, of.interest)
countOverlaps(of.interest, x)
subsetByOverlaps(x, of.interest)
subsetByOverlaps(of.interest, x)

# GRangesList overlaps:
of.interest <- GRangesList(first=regions(x), second=regions(x))
findOverlaps(of.interest, x)
findOverlaps(x, of.interest)
findOverlaps(x, of.interest, select="last")
overlapsAny(x, of.interest)
overlapsAny(of.interest, x)
countOverlaps(x, of.interest)
countOverlaps(of.interest, x)
subsetByOverlaps(x, of.interest)
subsetByOverlaps(of.interest, x)

# InteractionSet overlaps (split into two):
first.half <- x[1:Np/2]
second.half <- x[Np/2+1:Np/2]
findOverlaps(first.half, second.half)
findOverlaps(first.half, second.half, select="arbitrary")
overlapsAny(first.half, second.half)
countOverlaps(first.half, second.half)
subsetByOverlaps(first.half, second.half)
}

\author{
Aaron Lun
}

\seealso{
\code{\link{InteractionSet-class}},
\code{\link[GenomicRanges]{findOverlaps}}
}


