setGeneric("pairs", function(x, ...) standardGeneric("pairs"))

setAs("GInteractions", "SelfHits", function(from) {
      SelfHits(from=from@anchor1, to=from@anchor2, nnode=length(regions(from)), sort.by.query=FALSE)
})

.flipGI <- function(x) {
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
}

setAs("GInteractions", "GRangesList", function(from) .flipGI(from))


setMethod("pairs", "GInteractions", function(x, id=FALSE) {
    if (id) {
        return(as(x, "SelfHits"))
    } else {
        return(as(x, "GRangesList"))
    }
})

# Probably not to be used, as GRangesList may not always have two entries.
.unflipGI <- function(x, ...) {
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
}

