\name{Interaction overlaps}

\alias{overlapsAny,ContactMatrix,GRanges-method}

\title{Find overlaps between GRanges and a ContactMatrix}

\description{
Find overlaps between a set of linear intervals in a GRanges object, and the set of regions representing the rows or columns of a ContactMatrix.
}

\usage{
\S4method{overlapsAny}{ContactMatrix,GRanges}(query, subject, maxgap=0L, minoverlap=1L,
    type=c("any", "start", "end", "within", "equal"),
    algorithm=c("nclist", "intervaltree"),
    ignore.strand=TRUE)
}

\arguments{
\item{query}{A ContactMatrix object.}

\item{subject}{A GRanges object.}

\item{maxgap, minoverlap, type}{
    See \code{?\link[GenomicRanges]{findOverlaps}} in the \pkg{GenomicRanges} package.
}

\item{algorithm, ignore.strand}{
    See \code{?\link[GenomicRanges]{findOverlaps}} in the \pkg{GenomicRanges} package.
}
}

\details{
Overlaps are identified between the row regions of \code{query} and the regions in \code{subject}.
This is repeated for the column regions of \code{query}.
A list of two logical vectors is returned, specifying the rows and columns in \code{query} that are overlapped by any region in \code{subject}.

These vectors can be directly used to subset \code{query}.
Alternatively, they can be used in \code{outer} to generate a logical matrix for masking.
See the Examples section for some examples (obviously).
}

\examples{
example(ContactMatrix, echo=FALSE)

of.interest <- resize(sample(regions(x), 2), width=1, fix="center")
olap <- overlapsAny(x, of.interest)
olap
x[olap$row,]
x[,olap$column]
x[olap$row,olap$column]

keep <- outer(olap$row, olap$column, "|") # OR mask
temp <- as.matrix(x)
temp[!keep] <- NA

keep <- outer(olap$row, olap$column, "&") # AND mask
temp <- as.matrix(x)
temp[!keep] <- NA
}

\author{
Aaron Lun
}

\seealso{
\code{\link{ContactMatrix-class}},
\code{\link[GenomicRanges]{findOverlaps}}
}


