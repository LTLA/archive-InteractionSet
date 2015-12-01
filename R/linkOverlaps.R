# This defines the linkOverlaps method. The idea is to use interactions
# to link any two sets of regions together; to figure out which ones are
# linked, and to identify the interactions linking them.

.linkOverlap <- function(x, r1, r2, ...) {
    npairs <- length(x)
    a1 <- anchors(x, id=TRUE, type="first")
    a2 <- anchors(x, id=TRUE, type="second")
    nregs <- length(regions(x))

    olap1 <- .fast_overlap(x, r1, ..., IS.query=TRUE)
    bounds1 <- .get_iset_bounds(olap1, nregs)
    if (missing(r2)) { 
        olap2 <- olap1
        bounds2 <- bounds1
        is.same <- TRUE
    } else {
        olap2 <- .fast_overlap(x, r2, ..., IS.query=TRUE)
        bounds2 <- .get_iset_bounds(olap2, nregs)
        is.same <- FALSE
    }

    out <- .Call(cxx_expand_pair_links, a1 - 1L, a2 - 1L, 
                 bounds1$first - 1L, bounds1$last, olap1$ranges.dex - 1L, 
                 bounds2$first - 1L, bounds2$last, olap2$ranges.dex - 1L,
                 npairs, is.same)
    if (is.character(out)) { stop(out) }
    return(data.frame(query=out[[1]]+1L, subject1=out[[2]]+1L, subject2=out[[3]]+1L))    
}

setGeneric("linkOverlaps", function(query, subject1, subject2, ...) { standardGeneric('linkOverlaps') })
setMethod("linkOverlaps", c("GInteractions", "GRanges", "GRanges"), function(query, subject1, subject2, ...) {
    .linkOverlap(query, subject1, subject2, ...)
})

setMethod("linkOverlaps", c("GInteractions", "GRanges", "missing"), function(query, subject1, subject2, ...) {
    .linkOverlap(query, subject1, subject2, ...)
})

setMethod("linkOverlaps", c("InteractionSet", "GRanges", "GRanges"), function(query, subject1, subject2, ...) {
    .linkOverlap(query, subject1, subject2, ...)
})

setMethod("linkOverlaps", c("InteractionSet", "GRanges", "missing"), function(query, subject1, subject2, ...) {
    .linkOverlap(query, subject1, subject2, ...)
})


