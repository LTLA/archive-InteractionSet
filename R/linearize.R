.linearize <- function(x, ref, internal=TRUE) {
    ref <- as.integer(ref)
    keep.a1 <- anchors(x, type="first", id=TRUE) %in% ref
    keep.a2 <- anchors(x, type="second", id=TRUE) %in% ref

    if (!internal) {
        keep <- keep.a1 != keep.a2
    } else {
        keep <- keep.a1 | keep.a2
    }
    
    x <- x[keep,]
    new.ranges <- anchors(x, type="first") 
    keep.a1 <- keep.a1[keep]
    new.ranges[keep.a1] <- anchors(x[keep.a1,], type="second") # Replace anchor1 matches to "ref" with (mostly) non-ref anchor2.

    if (internal) {
        keep.a2 <- keep.a2[keep]
        is.same <- keep.a1 & keep.a2
        expanded <- range(pairs(x[is.same], as.grlist=TRUE)) # replace internal combos with ranges.
        if (any(lengths(expanded)>1L)) { stop("multi-chromosome sets of 'ref' are not supported") }
        new.ranges[is.same] <- unlist(expanded)
    }

    metadata(new.ranges) <- metadata(x)    
    mcols(new.ranges) <- cbind(mcols(new.ranges), mcols(x))   
    return(list(ranges=new.ranges, keep=keep))
}

setMethod("linearize", c("GInteractions", "numeric"), function(x, ref, internal=TRUE) {
    out <- .linearize(x, ref, internal=internal)
    return(out$ranges)
})

.chooseRegion <- function(x, ref, ...) {
    which(overlapsAny(regions(x), ref, ...))
}

setMethod("linearize", c("GInteractions", "GRanges"), function(x, ref, ..., internal=TRUE) {
    i <- .chooseRegion(x, ref, ...)          
    linearize(x, i, internal=internal)
})

.ISetToRSE <- function(x, out) {
    new("RangedSummarizedExperiment", x[out$keep,], rowRanges=out$ranges)
}

setMethod("linearize", c("InteractionSet", "numeric"), function(x, ref, internal=TRUE) {
    y <- interactions(x)
    out <- .linearize(y, ref, internal=internal)
    .ISetToRSE(x, out)
})

setMethod("linearize", c("InteractionSet", "GRanges"), function(x, ref, ..., internal=TRUE) {
    y <- interactions(x)
    i <- .chooseRegion(x, ref, ...)          
    out <- .linearize(y, i, internal=internal)
    .ISetToRSE(x, out)
})

