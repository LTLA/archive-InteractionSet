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

.decode_region_mode <- function(use.region, possibilities=c("both", "first", "second"))
# Gets the index of possibilities, for easy entry into C++.
# 1 -> both, 2 -> first, 3 -> second.
{
    use.region <- match.arg(use.region, possibilities)
    match(use.region, possibilities)
}

.linear_olap_finder <- function(gi, ranges, cxxfun, ..., gi.is.query=TRUE, use.region="both")
# Identifies linear overlaps, with differing C++ function depending on whether
# all overlaps are desired, or counting should be performed, etc.
{
    olap <- .fast_overlap(gi, ranges, ..., gi.is.query=gi.is.query)
    a1 <- anchors(gi, type="first", id=TRUE)
    a2 <- anchors(gi, type="second", id=TRUE)

    # Getting all combinations of overlaps (zero-indexing for C code).
    bounds <- .get_olap_bounds(olap, length(regions(gi)))
    out <- .Call(cxxfun, a1 - 1L, a2 - 1L, bounds$first - 1L, bounds$last, 
                 olap$ranges.dex - 1L, length(ranges), .decode_region_mode(use.region))
    if (is.character(out)) { stop(out) }
    return(out)
}

setMethod("findOverlaps", c(query="GInteractions", subject="Vector"), 
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=FALSE, use.region="both") {
        type <- match.arg(type)
        select <- match.arg(select)
        out <- .linear_olap_finder(query, subject, cxx_expand_olaps, 
                    maxgap=maxgap, minoverlap=minoverlap, type=type, 
                    ignore.strand=ignore.strand, gi.is.query=TRUE,
                    use.region=use.region)
        final <- Hits(out[[1]]+1L, out[[2]]+1L, length(query), length(subject), sort.by.query=TRUE)
        return(selectHits(final, select=match.arg(select))) 
    }
)

setMethod("findOverlaps", c(query="Vector", subject="GInteractions"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=FALSE, use.region="both") {
        type <- match.arg(type)
        select <- match.arg(select)
        out <- .linear_olap_finder(subject, query, cxx_expand_olaps, 
                    maxgap=maxgap, minoverlap=minoverlap, type=type, 
                    ignore.strand=ignore.strand, gi.is.query=FALSE,
                    use.region=use.region)
        final <- Hits(out[[2]]+1L, out[[1]]+1L, length(query), length(subject), sort.by.query=TRUE)
        final <- sort(final) 
        return(selectHits(final, select=match.arg(select))) 
    }
)

###############################################################

.paired_overlap_finder2 <- function(gi.left, gi.right, cxxfun, ..., use.region="both") 
# Identifies overlaps between two GI objects (left and right),
# with different C++ functions for different behaviours as before.
{
    used2 <- .get_used(gi.right)
    olap <- .fast_overlap(gi.left, regions(gi.right)[used2], ..., gi.is.query=TRUE)
    left.a1 <- anchors(gi.left, type="first", id=TRUE)
    left.a2 <- anchors(gi.left, type="second", id=TRUE)
    left.bounds <- .get_olap_bounds(olap, length(regions(gi.left)))
    olap$ranges.dex <- used2[olap$ranges.dex]

    right.a1 <- anchors(gi.right, type="first", id=TRUE)
    o1 <- order(right.a1)
    right.bounds1 <- .get_olap_bounds(list(gi.dex=right.a1[o1]), length(regions(gi.right)))
    right.a2 <- anchors(gi.right, type="second", id=TRUE)
    o2 <- order(right.a2)
    right.bounds2 <- .get_olap_bounds(list(gi.dex=right.a2[o2]), length(regions(gi.right)))
   
    # Getting all 2D overlaps.
    npairs <- length(gi.right)
    out <- .Call(cxxfun, left.a1 - 1L, left.a2 - 1L, 
                 left.bounds$first - 1L, left.bounds$last, olap$ranges.dex - 1L, 
                 right.bounds1$first - 1L, right.bounds1$last, o1 - 1L,
                 right.bounds2$first - 1L, right.bounds2$last, o2 - 1L,
                 npairs, .decode_region_mode(use.region, c("both", "same", "reverse")))
    if (is.character(out)) { stop(out) }
    return(out)
}

setMethod("findOverlaps", c(query="GInteractions", subject="GInteractions"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=FALSE, use.region="both") {
        type <- match.arg(type)
        select <- match.arg(select)
        out <- .paired_overlap_finder2(query, subject, cxx_expand_paired_olaps, 
                    maxgap=maxgap, minoverlap=minoverlap, type=type, 
                    ignore.strand=ignore.strand, use.region=use.region)
        final <- Hits(out[[1]]+1L, out[[2]]+1L, length(query), length(subject), sort.by.query=TRUE)
        return(selectHits(final, select=match.arg(select)))
    }
)

setMethod("findOverlaps", c(query="GInteractions", subject="missing"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=FALSE, use.region="both") {
        type <- match.arg(type)
        select <- match.arg(select)
        findOverlaps(query, query, maxgap=maxgap, minoverlap=minoverlap, 
                     type=type, select=select, ignore.strand=ignore.strand,
                     use.region=use.region)
    }
)

###############################################################
# This defines the countOverlaps method.

setMethod("countOverlaps", c(query="GInteractions", subject="Vector"), 
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=FALSE, use.region="both") {
        type <- match.arg(type)
        .linear_olap_finder(query, subject, cxx_queryhit_olaps,
            maxgap=maxgap, minoverlap=minoverlap, type=type,  
            ignore.strand=ignore.strand, gi.is.query=TRUE,
            use.region=use.region)
    }
)

setMethod("countOverlaps", c(query="Vector", subject="GInteractions"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=FALSE, use.region="both") {
        type <- match.arg(type)
        .linear_olap_finder(subject, query, cxx_subjecthit_olaps,
            maxgap=maxgap, minoverlap=minoverlap, type=type, 
            ignore.strand=ignore.strand, gi.is.query=FALSE,
            use.region=use.region)
    }
)

setMethod("countOverlaps", c(query="GInteractions", subject="GInteractions"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=FALSE, use.region="both") {
        type <- match.arg(type)
        out <- .paired_overlap_finder2(query, subject, cxx_queryhit_paired_olaps, 
                    maxgap=maxgap, minoverlap=minoverlap, type=type, 
                    ignore.strand=ignore.strand, use.region=use.region)
        return(out)
    }
)

setMethod("countOverlaps", c(query="GInteractions", subject="missing"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=FALSE, use.region='both') {
        type <- match.arg(type)
        countOverlaps(query, query, maxgap=maxgap, minoverlap=minoverlap,
                        type=type, ignore.strand=ignore.strand,
                        use.region=use.region)
    }
)

###############################################################
# This defines the overlapsAny method.

setMethod("overlapsAny", c(query="GInteractions", subject="Vector"), 
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=FALSE, use.region="both") {
        type <- match.arg(type)

        # Slightly faster than 'countOverlaps > 0', as there's no need to enumerate.
        rquery <- regions(query)
        keep <- logical(length(rquery))
        subset <- .get_used(query)
        if (length(subset)!=length(rquery)) { rquery <- rquery[subset] }
        keep[subset] <- overlapsAny(rquery, subject, maxgap=maxgap, minoverlap=minoverlap, 
            type=type, ignore.strand=ignore.strand)

        akeep1 <- keep[anchors(query, type="first", id=TRUE)]
        akeep2 <- keep[anchors(query, type="second", id=TRUE)]
        out <- .decode_region_mode(use.region)
        if (out==1L) { 
            return(akeep1 | akeep2)
        } else if (out==2L) {
            return(akeep1)
        } else {
            return(akeep2)
        }
    }
)

setMethod("overlapsAny", c(query="Vector", subject="GInteractions"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=FALSE, use.region='both') {

        type <- match.arg(type)
        out <- .decode_region_mode(use.region)
        if (out==1L) { 
            subset <- .get_used(subject)
        } else if (out==2L) {
            subset <- unique(anchors(subject, type="first", id=TRUE))
        } else { 
            subset <- unique(anchors(subject, type="second", id=TRUE))
        }
        overlapsAny(query, regions(subject)[subset], maxgap=maxgap, minoverlap=minoverlap, type=type, 
                    ignore.strand=ignore.strand) 
    }
)

setMethod("overlapsAny", c(query="GInteractions", subject="missing"),
    function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=FALSE, use.region="both") {
        return(!logical(length(query)))
    }
)

setMethod("overlapsAny", c(query="GInteractions", subject="GInteractions"),
    function(query, subject, maxgap=0L, minoverlap=1L,
            type=c("any", "start", "end", "within", "equal"),
            ignore.strand=FALSE, use.region="both") {
        type <- match.arg(type)
        return(countOverlaps(query, subject, maxgap=maxgap, minoverlap=minoverlap, 
                             type=type, ignore.strand=ignore.strand, 
                             use.region=use.region) > 0L)
})

###############################################################
# This defines the subsetByOverlaps method.

for (siglist in list(
        c(query="GInteractions", subject="Vector"), 
        c(query="Vector", subject="GInteractions"),
        c(query="GInteractions", subject="GInteractions")
    )) { 
    setMethod("subsetByOverlaps", siglist, function(query, subject, maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=FALSE, use.region='both') {
        type <- match.arg(type)
        query[overlapsAny(query, subject, maxgap=maxgap, minoverlap=minoverlap, 
                type=type, ignore.strand=ignore.strand, use.region=use.region),] 
    })
}

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

    full.call <- ""
    if ("select" %in% all.others) {
        full.call <- paste0(full.call, "select <- match.arg(select)\n")        
    } 
    if ("type" %in% all.others) {
        full.call <- paste0(full.call, "type <- match.arg(type)\n")        
    }
    full.call <- sprintf("{ %s%s(%s, %s) }", full.call, fun, internals, combined)
                                           
    out <- function() { }
    formals(out) <- c(alist(query=, subject=), other.args)
    body(out) <- parse(text=full.call)
    return(out)
}

for (siglist in list(
        c(query="InteractionSet", subject="Vector"), 
        c(query="Vector", subject="InteractionSet"),
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
             ignore.strand=FALSE, use.region='both')))
    setMethod("countOverlaps", siglist, olap.fun.gen(first.IS, second.IS, "countOverlaps", 
             alist(maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=FALSE, use.region='both')))
    setMethod("findOverlaps", siglist, olap.fun.gen(first.IS, second.IS, "findOverlaps", 
             alist(maxgap=0L, minoverlap=1L, 
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=FALSE, use.region='both')))

    if (is.na(second.IS)) { next }
    if (first.IS) {
        # Special treatment here, otherwise it'll return a GInteractions.
        setMethod("subsetByOverlaps", siglist, function(query, subject, maxgap=0L, minoverlap=1L, 
            type=c("any", "start", "end", "within", "equal"),
            ignore.strand=FALSE, use.region='both') {
            query[overlapsAny(query, subject, maxgap=maxgap, minoverlap=minoverlap,
                              type=type, ignore.strand=ignore.strand, use.region=use.region),]
        })
    } else {
        setMethod("subsetByOverlaps", siglist, olap.fun.gen(first.IS, second.IS, "subsetByOverlaps", 
            alist(maxgap=0L, minoverlap=1L, 
            type=c("any", "start", "end", "within", "equal"),
            ignore.strand=FALSE, use.region='both')))
    }
}

###############################################################
# Defining overlapsAny for ContactMatrix objects.

setMethod("overlapsAny", c("ContactMatrix", "Vector"), 
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
        c(query="ContactMatrix", subject="GInteractions"), 
        c(query="ContactMatrix", subject="InteractionSet")
    )) {
    setMethod("overlapsAny", siglist,
        function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"),
             ignore.strand=FALSE, use.region='both') {
        # It's possible to do this more efficiently by avoiding instantiation of the full object.
        # But it would require a total re-implementation at the C++ level, which is a pain.
        row.a <- rep(anchors(query, type="row", id=TRUE), ncol(query))
        col.a <- rep(anchors(query, type="column", id=TRUE), each=nrow(query))
        new.query <- GInteractions(row.a, col.a, regions(query)) 
        out <- overlapsAny(new.query, subject, maxgap=maxgap, minoverlap=minoverlap, type=type,
                           ignore.strand=ignore.strand, use.region=use.region)
        dim(out) <- dim(query)
        return(out)
    })
}

# Haven't defined the converse methods, as it's not clear whether you want to consider the entire
# interaction space in the ContactMatrix, or just the non-NA entries. 

###############################################################
# End
