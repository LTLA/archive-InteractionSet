# This defines the linkOverlaps method. The idea is to use interactions
# to link any two sets of regions together; to figure out which ones are
# linked, and to identify the interactions linking them.

.linkOverlap <- function(query, subject1, subject2, ..., use.region="both") {
    a1 <- anchors(query, id=TRUE, type="first")
    a2 <- anchors(query, id=TRUE, type="second")
    nregs <- length(regions(query))

    olap1 <- .fast_overlap(query, subject1, ..., gi.is.query=TRUE)
    bounds1 <- .get_olap_bounds(olap1, nregs)
    nregs1 <- length(subject1)
    if (missing(subject2)) { 
        olap2 <- olap1
        bounds2 <- bounds1
        is.same <- TRUE
        nregs2 <- nregs1
    } else {
        olap2 <- .fast_overlap(query, subject2, ..., gi.is.query=TRUE)
        bounds2 <- .get_olap_bounds(olap2, nregs)
        is.same <- FALSE
        nregs2 <- length(subject2)
    }

    out <- .Call(cxx_expand_pair_links, a1 - 1L, a2 - 1L, 
                 bounds1$first - 1L, bounds1$last, olap1$ranges.dex - 1L, nregs1,
                 bounds2$first - 1L, bounds2$last, olap2$ranges.dex - 1L, nregs2,
                 is.same, .decode_region_mode(use.region, c("both", "same", "reverse")))
    if (is.character(out)) { stop(out) }
    return(data.frame(query=out[[1]]+1L, subject1=out[[2]]+1L, subject2=out[[3]]+1L))    
}

setGeneric("linkOverlaps", function(query, subject1, subject2, ...) standardGeneric('linkOverlaps'))
setMethod("linkOverlaps", c("GInteractions", "GRanges", "GRanges"), .linkOverlap)
setMethod("linkOverlaps", c("GInteractions", "GRanges", "missing"), .linkOverlap)
setMethod("linkOverlaps", c("InteractionSet", "GRanges", "GRanges"), .linkOverlap)
setMethod("linkOverlaps", c("InteractionSet", "GRanges", "missing"), .linkOverlap)

