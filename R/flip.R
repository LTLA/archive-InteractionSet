# Flipping to a GRangesList (and back):

setGeneric("flip", function(x, ...) { standardGeneric("flip") })

setMethod("flip", "GInteractions", function(x) {
    in.order <- as.vector(do.call(rbind, anchors(x, id=TRUE)))
    status <- rep(c('first', 'second'), length(x))
    pairings <- rep(seq_along(x), each=2)
    all.regions <- regions(x)[in.order]
    names(all.regions) <- status

    out <- split(all.regions, pairings)
    mcols(out) <- mcols(x)
    names(out) <- names(x)
    metadata(out) <- metadata(x)
    return(out)
})

setMethod("flip", "GRangesList", function(x, ...) {
    if (!all(lengths(x)==2L)) { 
        stop("'x' can only contain GRanges of length 2") 
    }
    everything <- unname(unlist(x))
    all.first <- seq(1, length(x)*2L, by=2)
    all.second <- seq(2, length(x)*2L, by=2)

    out <- GInteractions(everything[all.first], everything[all.second], ...)
    mcols(out) <- mcols(x)
    metadata(out) <- metadata(x)
    names(out) < names(x)
    return(out)
})

