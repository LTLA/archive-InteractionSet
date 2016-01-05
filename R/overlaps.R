###############################################################
# This defines the findOverlaps method for GInteractions objects.

.get_used <- function(gi) 
# Gets the indices for all gi@regions that are actually used as anchors.
{
    all1 <- anchors(gi, type="first", id=TRUE)
    all2 <- anchors(gi, type="second", id=TRUE)
    used <- logical(length(regions(gi)))
    used[all1] <- TRUE
    used[all2] <- TRUE
    which(used)
}

.fast_overlap <- function(gi, ranges, ..., gi.is.query=TRUE) 
# Overlaps gi@regions with ranges, but only for those regions used as anchors.
{
    regs <- regions(gi)
    subset <- .get_used(gi)
    if (length(subset)!=length(regs)) { regs <- regs[subset] }

    # Need this behaviour, as type="within" will vary depending on manner of query vs. subject.        
    if (gi.is.query) { 
        olap <- findOverlaps(regs, ranges, select="all", ...)
        gi.dex <- subset[queryHits(olap)]
        ranges.dex <- subjectHits(olap)
    } else {
        olap <- findOverlaps(ranges, regs, select="all", ...)
        gi.dex <- subset[subjectHits(olap)]
        ranges.dex <- queryHits(olap)
        o <- order(gi.dex, ranges.dex)
        gi.dex <- gi.dex[o]
        ranges.dex <- ranges.dex[o]
    }

    return(list(gi.dex=gi.dex, ranges.dex=ranges.dex))
}

.get_olap_bounds <- function(olap, N) 
# Gets the first and last index of the overlap vector for each GI index.
{
    current.rle <- rle(olap$gi.dex)
    first.in.rle <- rep(1L, N)
    last.in.rle <- integer(N)
    cum.end <- cumsum(current.rle$lengths)
    first.in.rle[current.rle$values] <- cum.end - current.rle$lengths + 1L
    last.in.rle[current.rle$values] <- cum.end
    return(list(first=first.in.rle, last=last.in.rle))
}

.linear_olap_finder <- function(gi, ranges, cxxfun, ..., gi.is.query=TRUE) 
# Identifies linear overlaps, with differing C++ function depending on whether
# all overlaps are desired, or counting should be performed, etc.
{
    olap <- .fast_overlap(gi, ranges, ..., gi.is.query=gi.is.query)
    a1 <- anchors(gi, type="first", id=TRUE)
    a2 <- anchors(gi, type="second", id=TRUE)

    # Getting all combinations of overlaps (zero-indexing for C code).
    bounds <- .get_olap_bounds(olap, length(regions(gi)))
    out <- .Call(cxxfun, a1 - 1L, a2 - 1L, bounds$first - 1L, bounds$last, 
                 olap$ranges.dex - 1L, length(ranges))
    if (is.character(out)) { stop(out) }
    return(out)
}

setMethod("findOverlaps", c(query="GInteractions", subject="GRanges"), 
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=FALSE) {

        out <- .linear_olap_finder(query, subject, cxx_expand_olaps, 
                    maxgap=maxgap, minoverlap=minoverlap, type=type, 
                    ignore.strand=ignore.strand, gi.is.query=TRUE)
        final <- Hits(out[[1]]+1L, out[[2]]+1L, nrow(query), length(subject))
        return(selectHits(final, select=match.arg(select))) 
    }
)

setMethod("findOverlaps", c(query="GRanges", subject="GInteractions"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=FALSE) {

        out <- .linear_olap_finder(subject, query, cxx_expand_olaps, 
                    maxgap=maxgap, minoverlap=minoverlap, type=type, 
                    ignore.strand=ignore.strand, gi.is.query=FALSE)
        final <- Hits(out[[2]]+1L, out[[1]]+1L, length(query), nrow(subject))
        final <- sort(final) 
        return(selectHits(final, select=match.arg(select))) 
    }
)

###############################################################

.paired_overlap_finder <- function(gi, pairings, cxxfun, ..., gi.is.query=TRUE) 
# Identifies overlaps between interactions in GIs and paired regions,
# with different C++ functions for different behaviours as before.
{
    if (length(pairings)!=2L) { stop("input GRangesList must be of length 2") }
    npairs <- length(pairings[[1]])
    if (npairs!=length(pairings[[2]])) { stop("component GRanges in the GRangesList must be of the same length") }
    
    olap1 <- .fast_overlap(gi, pairings[[1]], ..., gi.is.query=gi.is.query)
    olap2 <- .fast_overlap(gi, pairings[[2]], ..., gi.is.query=gi.is.query)
    a1 <- anchors(gi, type="first", id=TRUE)
    a2 <- anchors(gi, type="second", id=TRUE)
    
    # Getting all 2D overlaps.
    bounds1 <- .get_olap_bounds(olap1, length(regions(gi)))
    bounds2 <- .get_olap_bounds(olap2, length(regions(gi)))
    out <- .Call(cxxfun, a1 - 1L, a2 - 1L, 
                 bounds1$first - 1L, bounds1$last, olap1$ranges.dex - 1L, 
                 bounds2$first - 1L, bounds2$last, olap2$ranges.dex - 1L,
                 npairs)
    if (is.character(out)) { stop(out) }
    return(out)
}

setMethod("findOverlaps", c(query="GInteractions", subject="GRangesList"), 
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=FALSE) {

        out <- .paired_overlap_finder(query, subject, cxx_expand_paired_olaps,
                    maxgap=maxgap, minoverlap=minoverlap, type=type, 
                    ignore.strand=ignore.strand, gi.is.query=TRUE)
       
        final <- Hits(out[[1]]+1L, out[[2]]+1L, nrow(query), length(subject[[1]])) # Cleaning up (1-indexing).
        return(selectHits(final, select=match.arg(select)))
    }
)

setMethod("findOverlaps", c(query="GRangesList", subject="GInteractions"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=FALSE) {

        out <- .paired_overlap_finder(subject, query, cxx_expand_paired_olaps,
                    maxgap=maxgap, minoverlap=minoverlap, type=type, 
                    ignore.strand=ignore.strand, gi.is.query=FALSE)

        final <- Hits(out[[2]]+1L, out[[1]]+1L, length(query[[1]]), nrow(subject)) 
        final <- sort(final) 
        return(selectHits(final, select=match.arg(select)))
    }
)

###############################################################

.reindex_by_anchor <- function(ref.list, all.anchors) 
# 'ref.list' stores the map from indices of region(subject) -> multiple indices of region(query) 
# This function generates the map from indices of anchor(subject) -> multiple indices of region(query)
{
    ref.list <- ref.list[all.anchors]
    left.indices <- as.integer(unlist(ref.list)) # coerces to 'integer(0)' if all.anchors is empty.
    anchor.indices <- rep(seq_along(all.anchors), lengths(ref.list))
    o <- order(left.indices, anchor.indices)
    return(list(gi.dex=left.indices[o], ranges.dex=anchor.indices[o]))
}

.paired_overlap_finder2 <- function(gi.left, gi.right, cxxfun, ...) 
# Identifies overlaps between two GI objects (left and right),
# with different C++ functions for different behaviours as before.
{
    used2 <- .get_used(gi.right)
    olap <- .fast_overlap(gi.left, regions(gi.right)[used2], ..., gi.is.query=TRUE)
    olap$ranges.dex <- used2[olap$ranges.dex]
    sub.list <- split(olap$gi.dex, olap$ranges.dex)
    ref.list <- rep(list(integer(0)), length(regions(gi.right)))
    ref.list[as.integer(names(sub.list))] <- sub.list

    # Reconstructing, as if we had done .fast_overlap on the anchor ranges directly.
    as1 <- anchors(gi.right, type="first", id=TRUE)
    olap1 <- .reindex_by_anchor(ref.list, as1)
    as2 <- anchors(gi.right, type="second", id=TRUE)
    olap2 <- .reindex_by_anchor(ref.list, as2)
   
    aq1 <- anchors(gi.left, type="first", id=TRUE)
    aq2 <- anchors(gi.left, type="second", id=TRUE)

    # Getting all 2D overlaps.
    npairs <- nrow(gi.right)
    bounds1 <- .get_olap_bounds(olap1, length(regions(gi.left)))
    bounds2 <- .get_olap_bounds(olap2, length(regions(gi.left)))
    out <- .Call(cxxfun, aq1 - 1L, aq2 - 1L, 
                 bounds1$first - 1L, bounds1$last, olap1$ranges.dex - 1L, 
                 bounds2$first - 1L, bounds2$last, olap2$ranges.dex - 1L,
                 npairs)
    if (is.character(out)) { stop(out) }
    return(out)
}

setMethod("findOverlaps", c(query="GInteractions", subject="GInteractions"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=FALSE) {
        out <- .paired_overlap_finder2(query, subject, cxx_expand_paired_olaps, 
                    maxgap=maxgap, minoverlap=minoverlap, type=type, 
                    ignore.strand=ignore.strand)
        final <- Hits(out[[1]]+1L, out[[2]]+1L, nrow(query), nrow(subject))
        return(selectHits(final, select=match.arg(select)))
    }
)

setMethod("findOverlaps", c(query="GInteractions", subject="missing"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=FALSE) {
        findOverlaps(query, query, maxgap=maxgap, minoverlap=minoverlap, 
                     type=type, select=select, ignore.strand=ignore.strand)
    }
)

###############################################################
# This defines the countOverlaps method.

setMethod("countOverlaps", c(query="GInteractions", subject="GRanges"), 
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=FALSE) {

        .linear_olap_finder(query, subject, cxx_queryhit_olaps,
            maxgap=maxgap, minoverlap=minoverlap, type=type,  
            ignore.strand=ignore.strand, gi.is.query=TRUE)
    }
)

setMethod("countOverlaps", c(query="GRanges", subject="GInteractions"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=FALSE) {

        .linear_olap_finder(subject, query, cxx_subjecthit_olaps,
            maxgap=maxgap, minoverlap=minoverlap, type=type, 
            ignore.strand=ignore.strand, gi.is.query=FALSE)
    }
)

setMethod("countOverlaps", c(query="GInteractions", subject="GRangesList"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=FALSE) {
        out <- .paired_overlap_finder(query, subject, cxx_queryhit_paired_olaps, 
                    maxgap=maxgap, minoverlap=minoverlap, type=type, 
                    ignore.strand=ignore.strand, gi.is.query=TRUE)
        return(out)
    }
)

setMethod("countOverlaps", c(query="GRangesList", subject="GInteractions"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=FALSE) {
        out <- .paired_overlap_finder(subject, query, cxx_subjecthit_paired_olaps, 
                    maxgap=maxgap, minoverlap=minoverlap, type=type, 
                    ignore.strand=ignore.strand, gi.is.query=FALSE)
        return(out)
    }
)

setMethod("countOverlaps", c(query="GInteractions", subject="GInteractions"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=FALSE) {
        out <- .paired_overlap_finder2(query, subject, cxx_queryhit_paired_olaps, 
                    maxgap=maxgap, minoverlap=minoverlap, type=type, 
                    ignore.strand=ignore.strand)
        return(out)
    }
)

setMethod("countOverlaps", c(query="GInteractions", subject="missing"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=FALSE) {
        countOverlaps(query, query, maxgap=maxgap, minoverlap=minoverlap,
                        type=type, ignore.strand=ignore.strand)
    }
)

###############################################################
# This defines the overlapsAny method.

setMethod("overlapsAny", c(query="GInteractions", subject="GRanges"), 
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=FALSE) {
        # Slightly faster than 'countOverlaps > 0', as there's no need to enumerate.
        rquery <- regions(query)
        keep <- logical(length(rquery))
        subset <- .get_used(query)
        if (length(subset)!=length(rquery)) { rquery <- rquery[subset] }
        keep[subset] <- overlapsAny(rquery, subject, maxgap=maxgap, minoverlap=minoverlap, type=type, 
                                    ignore.strand=ignore.strand)
        return(keep[anchors(query, type="first", id=TRUE)] | keep[anchors(query, type="second", id=TRUE)])
    }
)

setMethod("overlapsAny", c(query="GRanges", subject="GInteractions"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=FALSE) {
        subset <- .get_used(subject)
        overlapsAny(query, regions(subject)[subset], maxgap=maxgap, minoverlap=minoverlap, type=type, 
                    ignore.strand=ignore.strand) 
    }
)

setMethod("overlapsAny", c(query="GInteractions", subject="missing"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=FALSE) {
        return(!logical(length(query)))
    }
)

for (siglist in list(
#        c(query="GInteractions", subject="GRanges"), 
#        c(query="GRanges", subject="GInteractions"),
        c(query="GInteractions", subject="GRangesList"),
        c(query="GRangesList", subject="GInteractions"),
        c(query="GInteractions", subject="GInteractions")
    )) { 
    setMethod("overlapsAny", siglist, function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=FALSE) {
        return(countOverlaps(query, subject,
                    maxgap=maxgap, minoverlap=minoverlap, type=type, 
                    ignore.strand=ignore.strand) > 0L)
    })
}

###############################################################
# This defines the subsetByOverlaps method.

for (siglist in list(
        c(query="GInteractions", subject="GRanges"), 
        c(query="GRanges", subject="GInteractions"),
        c(query="GInteractions", subject="GRangesList"),
        c(query="GInteractions", subject="GInteractions")
    )) { 
    setMethod("subsetByOverlaps", siglist, function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=FALSE) {
        query[overlapsAny(query, subject, maxgap=maxgap, minoverlap=minoverlap, 
                type=type, ignore.strand=ignore.strand),] 
    })
}

setMethod("subsetByOverlaps", c(query="GRangesList", subject="GInteractions"), 
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=FALSE) {
        keep <- overlapsAny(query, subject, maxgap=maxgap, minoverlap=minoverlap, 
                    type=type, ignore.strand=ignore.strand)
        query[[1]] <- query[[1]][keep]
        query[[2]] <- query[[2]][keep]
        return(query)
    }
)

###############################################################
# Defining corresponding functions for InteractionSet objects.
    
olap.fun.gen <- function(first, second, fun, other.args) 
# A function to generate overlap functions for InteractionSet objects.
# Basically replaces query/subject with interactions(...) as appropriate.
# Avoids the need to manually write out each function again.
{
    if (is.na(second)) {
        if (!first) { stop("first should be IS if second is NA") }
        internals <- "query=interactions(query)"
    } else if (first && second) {
        internals <- "query=interactions(query), subject=interactions(subject)"
    } else if (first) { 
        internals <- "query=interactions(query), subject=subject"
    } else {
        internals <- "query=query, subject=interactions(subject)"
    }
    
    all.others <- names(other.args)
    if ("..."  %in% all.others) {  # Just in case.
        ellipsis <- TRUE
        all.others <- setdiff(all.others, "...") 
    } else { 
        ellipsis <- FALSE 
    }
    
    combined <- paste(paste(all.others, "=", all.others), collapse=", ")
    if (ellipsis) { combined <- paste(combined, ", ...") }
    full.call <- sprintf("%s(%s, %s)", fun, internals, combined)
                                           
    out <- function() { }
    formals(out) <- c(alist(query=, subject=), other.args)
    body(out) <- parse(text=full.call)
    return(out)
}

for (siglist in list(
        c(query="InteractionSet", subject="GRanges"), 
        c(query="GRanges", subject="InteractionSet"),
        c(query="InteractionSet", subject="GRangesList"),
        c(query="GRangesList", subject="InteractionSet"),
        c(query="InteractionSet", subject="InteractionSet"),
        c(query="InteractionSet", subject="GInteractions"),
        c(query="GInteractions", subject="InteractionSet"),
        c(query="InteractionSet", subject="missing")
    )) {
    first.IS <- siglist[["query"]]=="InteractionSet"
    second.IS <- siglist[["subject"]]=="InteractionSet"
    if (!second.IS) { 
        if (siglist[["subject"]]=="missing") { second.IS <- NA }
    }

    setMethod("overlapsAny", siglist, olap.fun.gen(first.IS, second.IS, "overlapsAny", 
             alist(maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=FALSE)))
    setMethod("countOverlaps", siglist, olap.fun.gen(first.IS, second.IS, "countOverlaps", 
             alist(maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=FALSE)))
    setMethod("findOverlaps", siglist, olap.fun.gen(first.IS, second.IS, "findOverlaps", 
             alist(maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=FALSE)))

    if (is.na(second.IS)) { next }
    if (first.IS) {
        # Special treatment here, otherwise it'll return a GInteractions.
        setMethod("subsetByOverlaps", siglist, function(query, subject, maxgap=0L, minoverlap=1L, 
            type=c("any", "start", "end", "within", "equal"),
            ignore.strand=FALSE) {
            query[overlapsAny(query, subject, maxgap=maxgap, minoverlap=minoverlap,
                              type=type, ignore.strand=ignore.strand),]
        })
    } else {
        setMethod("subsetByOverlaps", siglist, olap.fun.gen(first.IS, second.IS, "subsetByOverlaps", 
            alist(maxgap=0L, minoverlap=1L, 
            type=c("any", "start", "end", "within", "equal"),
            ignore.strand=FALSE)))
    }
}

###############################################################
# Defining overlapsAny for ContactMatrix objects.

setMethod("overlapsAny", c("ContactMatrix", "GRanges"), 
    function(query, subject, maxgap=0L, minoverlap=1L,
        type=c("any", "start", "end", "within", "equal"),
        ignore.strand=FALSE) {
        a1 <- anchors(query, id=TRUE, type="row")
        a2 <- anchors(query, id=TRUE, type="column")
        
        is.used <- union(a1, a2)
        is.overlapped <- logical(length(regions(query)))
        is.overlapped[is.used] <- overlapsAny(regions(query)[is.used], subject, maxgap=maxgap,
                                        minoverlap=minoverlap, type=type, ignore.strand=ignore.strand)
        return(list(row=is.overlapped[a1], column=is.overlapped[a2]))
})

# Use outer(output$row, output$column, "|" or "&") to get the logical area in the interaction space.
# Not sure it makes a great deal of sense to define 'findOverlaps' here.

for (siglist in list(
        c(query="ContactMatrix", subject="GRangesList"), 
        c(query="ContactMatrix", subject="GInteractions"), 
        c(query="ContactMatrix", subject="InteractionSet")
    )) {
    setMethod("overlapsAny", siglist,
        function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=FALSE) {
    
        # It's possible to do this more efficiently by avoiding instantiation of the full object.
        # But it would require a total re-implementation at the C++ level, which is a pain.
        row.a <- rep(anchors(query, type="row", id=TRUE), ncol(query))
        col.a <- rep(anchors(query, type="column", id=TRUE), each=nrow(query))
        new.query <- GInteractions(row.a, col.a, regions(query)) 
        out <- overlapsAny(new.query, subject, maxgap=maxgap, minoverlap=minoverlap, type=type,
                           ignore.strand=ignore.strand)
        dim(out) <- dim(query)
        return(out)
    })
}

# Haven't defined the converse methods, as it's not clear whether you want to consider the entire
# interaction space in the ContactMatrix, or just the non-NA entries. 

###############################################################
# End
