setGeneric("swapAnchors", function(x, ...) { standardGeneric("swapAnchors") })
setMethod("swapAnchors", "GInteractions", function(x, mode=c("order", "reverse", "all")) {
    mode <- match.arg(mode)
    if (mode=="order") {
        out <- .enforce_order(x@anchor1, x@anchor2)
    } else if (mode=="reverse") {
        out <- rev(.enforce_order(x@anchor1, x@anchor2))
    } else {
        out <- list(x@anchor2, x@anchor1)
    }
    x@anchor1 <- out[[1]]        
    x@anchor2 <- out[[2]]
    validObject(x)
    return(x)
})

setMethod("swapAnchors", "InteractionSet", function(x, mode=c("order", "reverse", "all")) {
    x@interactions <- swapAnchors(x@interactions, mode)
    return(x)
})

