setMethod("swapAnchors", "GInteractions", function(x, mode=c("order", "reverse", "all")) {
    mode <- match.arg(mode)
    if (mode=="order") {
        out <- .enforce_order(anchor1(x), anchor2(x))
    } else if (mode=="reverse") {
        out <- rev(.enforce_order(anchor1(x), anchor2(x)))
    } else {
        out <- list(anchor2(x), anchor1(x))
    }
    x@anchor1 <- out[[1]]        
    x@anchor2 <- out[[2]]
    validObject(x)
    return(x)
})

setMethod("swapAnchors", "InteractionSet", function(x, mode=c("order", "reverse", "all")) {
    interactions(x) <- swapAnchors(interactions(x), mode)
    return(x)
})

