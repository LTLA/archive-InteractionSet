###############################################################
# This defines the findOverlaps method for InteractionSet objects.

.get_used <- function(query) {
    all1 <- anchors(query, type="first", id=TRUE)
    all2 <- anchors(query, type="second", id=TRUE)
    used <- logical(length(regions(query)))
    used[all1] <- TRUE
    used[all2] <- TRUE
    which(used)
}

.fast_overlap <- function(query, subject, ...) {
    rquery <- regions(query)
    subset <- .get_used(query)
    if (length(subset)!=length(rquery)) { rquery <- rquery[subset] }
    olap <- findOverlaps(rquery, subject, select="all", ...)
    return(list(query.dex=subset[queryHits(olap)], subject.dex=subjectHits(olap)))
}

.get_query_bounds <- function(olap, N) {
    qrle <- rle(olap$query.dex)
    first.query <- rep(1L, N)
    last.query <- integer(N)
    cum.end <- cumsum(qrle$lengths)
    first.query[qrle$values] <- cum.end - qrle$lengths + 1L
    last.query[qrle$values] <- cum.end
    return(list(first=first.query, last=last.query))
}

.reselect <- function(hits, select) {
    if (select=="all") { return(hits) }

    # Hits are sorted by query, so we can do this:    
    all.qs <- rle(queryHits(hits))
    actual.q <- all.qs$values
    stopifnot(!is.unsorted(actual.q))
    indices <- rep(as.integer(NA), queryLength(hits))

    last.of.q <- cumsum(all.qs$lengths)
    first.of.q <- last.of.q - all.qs$lengths + 1L
    if (select=="first") { 
        indices[actual.q] <- subjectHits(hits)[first.of.q]
    } else if (select=="last") { 
        indices[actual.q] <- subjectHits(hits)[last.of.q]
    } else {
        # 'arbitrary'. Who knows what this gives?
        indices[actual.q] <- subjectHits(hits)[round((first.of.q + last.of.q)/2L)]
    }
    return(indices)
}

setMethod("findOverlaps", c(query="InteractionSet", subject="GRanges"), 
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=TRUE) {
        olap <- .fast_overlap(query, subject, maxgap=maxgap, minoverlap=minoverlap,
                    type=type, algorithm=algorithm, ignore.strand=ignore.strand)
        a1 <- anchors(query, type="first", id=TRUE)
        a2 <- anchors(query, type="second", id=TRUE)

        # Getting all combinations of overlaps (zero-indexing for C code).
        qbounds <- .get_query_bounds(olap, length(regions(query)))
        out <- .Call("expand_olaps", a1 - 1L, a2 - 1L, qbounds$first - 1L, qbounds$last, 
                     olap$subject.dex - 1L, length(subject), PACKAGE="InteractionSet")
        if (is.character(out)) { stop(out) }

        # Cleaning up (1-indexing, and resorting).
        final <- Hits(out[[1]]+1L, out[[2]]+1L, nrow(query), length(subject))
        select <- match.arg(select)
        return(.reselect(final, select=select))
    }
)

setMethod("findOverlaps", c(query="GRanges", subject="InteractionSet"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=TRUE) {
        final <- t(findOverlaps(subject, query, maxgap=maxgap, minoverlap=minoverlap, select="all",
                       type=type, algorithm=algorithm, ignore.strand=ignore.strand))
        select <- match.arg(select)
        return(.reselect(final, select=select))
    }
)

###############################################################

.paired_overlap_finder <- function(query, subject, cxxfun, ...) 
# This is split into a separate function, because we'll re-use 
# the code to run 'overlapsAny'.
{
    stopifnot(length(subject)==2L)
    npairs <- length(subject[[1]])
    stopifnot(npairs==length(subject[[2]]))
    
    olap1 <- .fast_overlap(query, subject[[1]], ...)
    olap2 <- .fast_overlap(query, subject[[2]], ...)
    a1 <- anchors(query, type="first", id=TRUE)
    a2 <- anchors(query, type="second", id=TRUE)
    
    # Getting all 2D overlaps.
    qbounds1 <- .get_query_bounds(olap1, length(regions(query)))
    qbounds2 <- .get_query_bounds(olap2, length(regions(query)))
    out <- .Call(cxxfun, a1 - 1L, a2 - 1L, 
                 qbounds1$first - 1L, qbounds1$last, olap1$subject.dex - 1L, 
                 qbounds2$first - 1L, qbounds2$last, olap2$subject.dex - 1L,
                 npairs, PACKAGE="InteractionSet")
    if (is.character(out)) { stop(out) }
    return(out)
}

setMethod("findOverlaps", c(query="InteractionSet", subject="GRangesList"), 
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=TRUE) {
        out <- .paired_overlap_finder(query, subject, "expand_paired_olaps", 
                    maxgap=maxgap, minoverlap=minoverlap, type=type, 
                    algorithm=algorithm, ignore.strand=ignore.strand)
        final <- Hits(out[[1]]+1L, out[[2]]+1L, nrow(query), length(subject[[1]])) # Cleaning up (1-indexing and resorting).
        select <- match.arg(select)
        return(.reselect(final, select=select))
    }
)

setMethod("findOverlaps", c(query="GRangesList", subject="InteractionSet"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=TRUE) {
        final <- t(findOverlaps(subject, query, maxgap=maxgap, minoverlap=minoverlap, select="all",
                    type=type, algorithm=algorithm, ignore.strand=ignore.strand))
        select <- match.arg(select)
        return(.reselect(final, select=select))
    }
)

###############################################################

.reindex_by_anchor <- function(ref.list, anchor.dex) 
# 'ref.list' stores the map from indices of region(subject) -> multiple indices of region(query) 
# We want to get the map from indices of anchor(subject) -> multiple indices of region(query)
{
    ref.list <- ref.list[anchor.dex]
    all.queries <- unlist(ref.list)
    all.anchors <- rep(seq_along(anchor.dex), lengths(ref.list))
    o <- order(all.queries, all.anchors)
    return(list(query.dex=all.queries[o], subject.dex=all.anchors[o]))
}

.paired_overlap_finder2 <- function(query, subject, cxxfun, ...) 
# Again, splitting the code for re-use in 'overlapsAny'.
{
    used2 <- .get_used(subject)
    olap <- .fast_overlap(query, regions(subject)[used2], ...)
    olap$subject.dex <- used2[olap$subject.dex]
    sub.list <- split(olap$query.dex, olap$subject.dex)
    ref.list <- rep(list(integer(0)), length(regions(subject)))
    ref.list[as.integer(names(sub.list))] <- sub.list

    # Reconstructing, as if we had done .fast_overlap on the anchor ranges directly.
    as1 <- anchors(subject, type="first", id=TRUE)
    olap1 <- .reindex_by_anchor(ref.list, as1)
    as2 <- anchors(subject, type="second", id=TRUE)
    olap2 <- .reindex_by_anchor(ref.list, as2)
   
    aq1 <- anchors(query, type="first", id=TRUE)
    aq2 <- anchors(query, type="second", id=TRUE)

    # Getting all 2D overlaps.
    npairs <- nrow(subject)
    qbounds1 <- .get_query_bounds(olap1, length(regions(query)))
    qbounds2 <- .get_query_bounds(olap2, length(regions(query)))
    out <- .Call(cxxfun, aq1 - 1L, aq2 - 1L, 
                 qbounds1$first - 1L, qbounds1$last, olap1$subject.dex - 1L, 
                 qbounds2$first - 1L, qbounds2$last, olap2$subject.dex - 1L,
                 npairs, PACKAGE="InteractionSet")
    if (is.character(out)) { stop(out) }
    return(out)
}

setMethod("findOverlaps", c(query="InteractionSet", subject="InteractionSet"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=TRUE) {
        out <- .paired_overlap_finder2(query, subject, "expand_paired_olaps", 
                    maxgap=maxgap, minoverlap=minoverlap, type=type, 
                    algorithm=algorithm, ignore.strand=ignore.strand)
        final <- Hits(out[[1]]+1L, out[[2]]+1L, nrow(query), nrow(subject))
        select <- match.arg(select)
        return(.reselect(final, select=select))
    }
)

###############################################################
# This defines the overlapsAny method.

setMethod("overlapsAny", c(query="InteractionSet", subject="GRanges"), 
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=TRUE) {
        rquery <- regions(query)
        keep <- logical(length(rquery))
        subset <- .get_used(query)
        if (length(subset)!=length(rquery)) { rquery <- rquery[subset] }
        keep[subset] <- overlapsAny(rquery, subject, maxgap=maxgap, 
                            minoverlap=minoverlap, type=type, 
                            algorithm=algorithm, ignore.strand=ignore.strand)
        return(keep[anchors(query, type="first", id=TRUE)] | 
               keep[anchors(query, type="second", id=TRUE)])
    }
)

setMethod("overlapsAny", c(query="GRanges", subject="InteractionSet"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=TRUE) {
        subset <- .get_used(subject)
        return(overlapsAny(query, regions(subject)[subset], maxgap=maxgap, 
                minoverlap=minoverlap, type=type, algorithm=algorithm, 
                ignore.strand=ignore.strand)) 
    }
)

setMethod("overlapsAny", c(query="InteractionSet", subject="GRangesList"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=TRUE) {
        out <- .paired_overlap_finder(query, subject, "queryhit_paired_olaps", 
                    maxgap=maxgap, minoverlap=minoverlap, type=type, 
                    algorithm=algorithm, ignore.strand=ignore.strand)
        return(out>0L)
    }
)

setMethod("overlapsAny", c(query="GRangesList", subject="InteractionSet"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=TRUE) {
        out <- .paired_overlap_finder(subject, query, "subjecthit_paired_olaps", 
                    maxgap=maxgap, minoverlap=minoverlap, type=type, 
                    algorithm=algorithm, ignore.strand=ignore.strand)
        return(out>0L)
    }
)

setMethod("overlapsAny", c(query="InteractionSet", subject="InteractionSet"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=TRUE) {
        out <- .paired_overlap_finder2(query, subject, "queryhit_paired_olaps", 
                    maxgap=maxgap, minoverlap=minoverlap, type=type, 
                    algorithm=algorithm, ignore.strand=ignore.strand)
        return(out>0L)        
    }
)

###############################################################
# This defines the subsetByOverlaps method.

setMethod("subsetByOverlaps", c(query="InteractionSet", subject="Vector"), 
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=TRUE) {
        query[overlapsAny(query, subject, maxgap=maxgap, minoverlap=minoverlap, 
                type=type, algorithm=algorithm, ignore.strand=ignore.strand),] 
    }
)

setMethod("subsetByOverlaps", c(query="GRangesList", subject="InteractionSet"), 
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=TRUE) {
        keep <- overlapsAny(query, subject, maxgap=maxgap, minoverlap=minoverlap, 
                    type=type, algorithm=algorithm, ignore.strand=ignore.strand)
        query[[1]] <- query[[1]][keep]
        query[[2]] <- query[[2]][keep]
        return(query)
    }
)

setMethod("subsetByOverlaps", c(query="GRanges", subject="InteractionSet"), 
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=TRUE) {
        query[overlapsAny(query, subject, maxgap=maxgap, minoverlap=minoverlap, 
                type=type, algorithm=algorithm, ignore.strand=ignore.strand),] 
    }
)

###############################################################
# End




