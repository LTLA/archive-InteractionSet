###############################################################
# This defines the findOverlaps method for InteractionSet objects.

.get_used <- function(iset) {
    all1 <- anchors(iset, type="first", id=TRUE)
    all2 <- anchors(iset, type="second", id=TRUE)
    used <- logical(length(regions(iset)))
    used[all1] <- TRUE
    used[all2] <- TRUE
    which(used)
}

.fast_overlap <- function(iset, ranges, ..., IS.query=TRUE) {
    regs <- regions(iset)
    subset <- .get_used(iset)
    if (length(subset)!=length(regs)) { regs <- regs[subset] }

    # Need this behaviour, as type="within" will vary depending on manner of query vs. subject.        
    if (IS.query) { 
        olap <- findOverlaps(regs, ranges, select="all", ...)
        iset.dex <- subset[queryHits(olap)]
        ranges.dex <- subjectHits(olap)
    } else {
        olap <- findOverlaps(ranges, regs, select="all", ...)
        iset.dex <- subset[subjectHits(olap)]
        ranges.dex <- queryHits(olap)
        o <- order(iset.dex, ranges.dex)
        iset.dex <- iset.dex[o]
        ranges.dex <- ranges.dex[o]
    }

    return(list(iset.dex=iset.dex, ranges.dex=ranges.dex))
}

.get_iset_bounds <- function(olap, N) {
    current.rle <- rle(olap$iset.dex)
    first.in.rle <- rep(1L, N)
    last.in.rle <- integer(N)
    cum.end <- cumsum(current.rle$lengths)
    first.in.rle[current.rle$values] <- cum.end - current.rle$lengths + 1L
    last.in.rle[current.rle$values] <- cum.end
    return(list(first=first.in.rle, last=last.in.rle))
}

.linear_olap_finder <- function(iset, ranges, cxxfun, ..., IS.query=TRUE) {
    olap <- .fast_overlap(iset, ranges, ..., IS.query=IS.query)
    a1 <- anchors(iset, type="first", id=TRUE)
    a2 <- anchors(iset, type="second", id=TRUE)

    # Getting all combinations of overlaps (zero-indexing for C code).
    bounds <- .get_iset_bounds(olap, length(regions(iset)))
    out <- .Call(cxxfun, a1 - 1L, a2 - 1L, bounds$first - 1L, bounds$last, 
                 olap$ranges.dex - 1L, length(ranges))
    if (is.character(out)) { stop(out) }
    return(out)
}

setMethod("findOverlaps", c(query="InteractionSet", subject="GRanges"), 
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=TRUE) {

        out <- .linear_olap_finder(query, subject, cxx_expand_olaps, 
                    maxgap=maxgap, minoverlap=minoverlap,
                    type=type, algorithm=algorithm, 
                    ignore.strand=ignore.strand, IS.query=TRUE)
        final <- Hits(out[[1]]+1L, out[[2]]+1L, nrow(query), length(subject))
        return(selectHits(final, select=match.arg(select))) 
    }
)

setMethod("findOverlaps", c(query="GRanges", subject="InteractionSet"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=TRUE) {

        out <- .linear_olap_finder(subject, query, cxx_expand_olaps, 
                    maxgap=maxgap, minoverlap=minoverlap,
                    type=type, algorithm=algorithm, 
                    ignore.strand=ignore.strand, IS.query=FALSE)
        final <- Hits(out[[2]]+1L, out[[1]]+1L, length(query), nrow(subject))
        final <- sort(final) 
        return(selectHits(final, select=match.arg(select))) 
    }
)

###############################################################

.paired_overlap_finder <- function(iset, pairings, cxxfun, ..., IS.query=TRUE) 
# This is split into a separate function, because we'll re-use 
# the code to run 'overlapsAny'.
{
    if (length(pairings)!=2L) { stop("input GRangesList must be of length 2") }
    npairs <- length(pairings[[1]])
    if (npairs!=length(pairings[[2]])) { stop("component GRanges in the GRangesList must be of the same length") }
    
    olap1 <- .fast_overlap(iset, pairings[[1]], ..., IS.query=IS.query)
    olap2 <- .fast_overlap(iset, pairings[[2]], ..., IS.query=IS.query)
    a1 <- anchors(iset, type="first", id=TRUE)
    a2 <- anchors(iset, type="second", id=TRUE)
    
    # Getting all 2D overlaps.
    bounds1 <- .get_iset_bounds(olap1, length(regions(iset)))
    bounds2 <- .get_iset_bounds(olap2, length(regions(iset)))
    out <- .Call(cxxfun, a1 - 1L, a2 - 1L, 
                 bounds1$first - 1L, bounds1$last, olap1$ranges.dex - 1L, 
                 bounds2$first - 1L, bounds2$last, olap2$ranges.dex - 1L,
                 npairs)
    if (is.character(out)) { stop(out) }
    return(out)
}

setMethod("findOverlaps", c(query="InteractionSet", subject="GRangesList"), 
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=TRUE) {

        out <- .paired_overlap_finder(query, subject, cxx_expand_paired_olaps,
                    maxgap=maxgap, minoverlap=minoverlap, type=type, 
                    algorithm=algorithm, ignore.strand=ignore.strand, IS.query=TRUE)
       
        final <- Hits(out[[1]]+1L, out[[2]]+1L, nrow(query), length(subject[[1]])) # Cleaning up (1-indexing).
        return(selectHits(final, select=match.arg(select)))
    }
)

setMethod("findOverlaps", c(query="GRangesList", subject="InteractionSet"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=TRUE) {

        out <- .paired_overlap_finder(subject, query, cxx_expand_paired_olaps,
                    maxgap=maxgap, minoverlap=minoverlap, type=type, 
                    algorithm=algorithm, ignore.strand=ignore.strand, IS.query=FALSE)

        final <- Hits(out[[2]]+1L, out[[1]]+1L, length(query[[1]]), nrow(subject)) 
        final <- sort(final) 
        return(selectHits(final, select=match.arg(select)))
    }
)

###############################################################

.reindex_by_anchor <- function(ref.list, all.anchors) 
# 'ref.list' stores the map from indices of region(subject) -> multiple indices of region(query) 
# We want to get the map from indices of anchor(subject) -> multiple indices of region(query)
{
    ref.list <- ref.list[all.anchors]
    left.indices <- as.integer(unlist(ref.list)) # coerces to 'integer(0)' if all.anchors is empty.
    anchor.indices <- rep(seq_along(all.anchors), lengths(ref.list))
    o <- order(left.indices, anchor.indices)
    return(list(iset.dex=left.indices[o], ranges.dex=anchor.indices[o]))
}

.paired_overlap_finder2 <- function(iset.left, iset.right, cxxfun, ...) 
# Again, splitting the code for re-use in 'overlapsAny'.
{
    used2 <- .get_used(iset.right)
    olap <- .fast_overlap(iset.left, regions(iset.right)[used2], ..., IS.query=TRUE)
    olap$ranges.dex <- used2[olap$ranges.dex]
    sub.list <- split(olap$iset.dex, olap$ranges.dex)
    ref.list <- rep(list(integer(0)), length(regions(iset.right)))
    ref.list[as.integer(names(sub.list))] <- sub.list

    # Reconstructing, as if we had done .fast_overlap on the anchor ranges directly.
    as1 <- anchors(iset.right, type="first", id=TRUE)
    olap1 <- .reindex_by_anchor(ref.list, as1)
    as2 <- anchors(iset.right, type="second", id=TRUE)
    olap2 <- .reindex_by_anchor(ref.list, as2)
   
    aq1 <- anchors(iset.left, type="first", id=TRUE)
    aq2 <- anchors(iset.left, type="second", id=TRUE)

    # Getting all 2D overlaps.
    npairs <- nrow(iset.right)
    bounds1 <- .get_iset_bounds(olap1, length(regions(iset.left)))
    bounds2 <- .get_iset_bounds(olap2, length(regions(iset.left)))
    out <- .Call(cxxfun, aq1 - 1L, aq2 - 1L, 
                 bounds1$first - 1L, bounds1$last, olap1$ranges.dex - 1L, 
                 bounds2$first - 1L, bounds2$last, olap2$ranges.dex - 1L,
                 npairs)
    if (is.character(out)) { stop(out) }
    return(out)
}

setMethod("findOverlaps", c(query="InteractionSet", subject="InteractionSet"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=TRUE) {
        out <- .paired_overlap_finder2(query, subject, cxx_expand_paired_olaps, 
                    maxgap=maxgap, minoverlap=minoverlap, type=type, 
                    algorithm=algorithm, ignore.strand=ignore.strand)
        final <- Hits(out[[1]]+1L, out[[2]]+1L, nrow(query), nrow(subject))
        return(selectHits(final, select=match.arg(select)))
    }
)

###############################################################
# This defines the countOverlaps method.

setMethod("countOverlaps", c(query="InteractionSet", subject="GRanges"), 
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=TRUE) {

        .linear_olap_finder(query, subject, cxx_queryhit_olaps,
            maxgap=maxgap, minoverlap=minoverlap,
            type=type, algorithm=algorithm, 
            ignore.strand=ignore.strand, IS.query=TRUE)
    }
)

setMethod("countOverlaps", c(query="GRanges", subject="InteractionSet"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=TRUE) {

        .linear_olap_finder(subject, query, cxx_subjecthit_olaps,
            maxgap=maxgap, minoverlap=minoverlap,
            type=type, algorithm=algorithm, 
            ignore.strand=ignore.strand, IS.query=FALSE)
    }
)

setMethod("countOverlaps", c(query="InteractionSet", subject="GRangesList"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=TRUE) {
        out <- .paired_overlap_finder(query, subject, cxx_queryhit_paired_olaps, 
                    maxgap=maxgap, minoverlap=minoverlap, type=type, 
                    algorithm=algorithm, ignore.strand=ignore.strand, IS.query=TRUE)
        return(out)
    }
)

setMethod("countOverlaps", c(query="GRangesList", subject="InteractionSet"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=TRUE) {
        out <- .paired_overlap_finder(subject, query, cxx_subjecthit_paired_olaps, 
                    maxgap=maxgap, minoverlap=minoverlap, type=type, 
                    algorithm=algorithm, ignore.strand=ignore.strand, IS.query=FALSE)
        return(out)
    }
)

setMethod("countOverlaps", c(query="InteractionSet", subject="InteractionSet"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=TRUE) {
        out <- .paired_overlap_finder2(query, subject, cxx_queryhit_paired_olaps, 
                    maxgap=maxgap, minoverlap=minoverlap, type=type, 
                    algorithm=algorithm, ignore.strand=ignore.strand)
        return(out)
    }
)

###############################################################
# This defines the overlapsAny method.

for (siglist in list(
        c(query="InteractionSet", subject="GRanges"), 
        c(query="GRanges", subject="InteractionSet"),
        c(query="InteractionSet", subject="GRangesList"),
        c(query="GRangesList", subject="InteractionSet"),
        c(query="InteractionSet", subject="InteractionSet")
    )) { 
    setMethod("overlapsAny", siglist, function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=TRUE) {
        return(countOverlaps(query, subject,
                    maxgap=maxgap, minoverlap=minoverlap, type=type, 
                    algorithm=algorithm, ignore.strand=ignore.strand) > 0L)
    })
}

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
