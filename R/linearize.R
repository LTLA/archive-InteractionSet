setGeneric("linearize", function(x, ref, ...) standardGeneric("linearize"))
setMethod("linearize", c("InteractionSet", "numeric"), function(x, ref, ...) {
    ref <- as.integer(ref)
    keep.a1 <- anchors(x, type="first", id=TRUE) == ref
    keep.a2 <- anchors(x, type="second", id=TRUE) == ref

    keep <- keep.a1 | keep.a2
    x <- x[keep,]
    new.ranges <- anchors(x, type="first") 
    keep.a1 <- keep.a1[keep]
    new.ranges[keep.a1] <- anchors(x[keep.a1,], type="second") # Replace anchor1 matches to "ref" with (mostly) non-ref anchor2.

    out <- SummarizedExperiment(assays=assays(x), rowRanges=new.ranges,
        colData=colData(x), metadata=metadata(x))
    elementMetadata(out) <- elementMetadata(x)
    return(out)    
})

setMethod("linearize", c("InteractionSet", "GRanges"), function(x, ref, ...) {
    i <- which(overlapsAny(regions(x), ref, ...))
    if (length(i)==0L) { 
        i <- 0L # Dummy value, returns empty RSE.
    } else if (length(i)>1L) { 
        warning("multiple matching reference regions, using the first region only")
        i <- i[1]
    }
    linearize(x, i)   
})
