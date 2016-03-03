setGeneric("linearize", function(x, ref, ...) standardGeneric("linearize"))

.linearize <- function(x, ref) {
    ref <- as.integer(ref)
    keep.a1 <- anchors(x, type="first", id=TRUE) == ref
    keep.a2 <- anchors(x, type="second", id=TRUE) == ref

    keep <- keep.a1 | keep.a2
    x <- x[keep,]
    new.ranges <- anchors(x, type="first") 
    keep.a1 <- keep.a1[keep]
    new.ranges[keep.a1] <- anchors(x[keep.a1,], type="second") # Replace anchor1 matches to "ref" with (mostly) non-ref anchor2.

    metadata(new.ranges) <- metadata(x)    
    mcols(new.ranges) <- cbind(mcols(new.ranges), mcols(x))   
    return(list(ranges=new.ranges, keep=keep))
}

setMethod("linearize", c("GInteractions", "numeric"), function(x, ref) {
    out <- .linearize(x, ref)
    return(out$ranges)
})

.chooseRegion <- function(x, ref, ...) {
    i <- which(overlapsAny(regions(x), ref, ...))
    if (length(i)==0L) { 
        i <- 0L # Dummy value, returns empty RSE.
    } else if (length(i)>1L) { 
        warning("multiple matching reference regions, using the first region only")
        i <- i[1]
    }
    return(i)
}

setMethod("linearize", c("GInteractions", "GRanges"), function(x, ref, ...) {
    i <- .chooseRegion(x, ref, ...)          
    linearize(x, i)   
})

.ISetToRSE <- function(x, out) {
    new("RangedSummarizedExperiment", x[out$keep,], rowRanges=out$ranges)
}

setMethod("linearize", c("InteractionSet", "numeric"), function(x, ref) {
    y <- interactions(x)
    out <- .linearize(y, ref)
    .ISetToRSE(x, out)
})

setMethod("linearize", c("InteractionSet", "GRanges"), function(x, ref, ...) {
    y <- interactions(x)
    i <- .chooseRegion(x, ref, ...)          
    out <- .linearize(y, i)
    .ISetToRSE(x, out)
})

