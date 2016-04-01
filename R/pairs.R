setGeneric("pairs", function(x, ...) standardGeneric("pairs"))

setAs("GInteractions", "SelfHits", function(from) {
      SelfHits(from=from@anchor1, to=from@anchor2, nnode=length(regions(from)), sort.by.query=FALSE)
})

.flipGI <- function(x) {
    in.order <- as.vector(do.call(rbind, anchors(x, id=TRUE)))
    all.regions <- regions(x)[in.order]
    names(all.regions) <- rep(c('first', 'second'), length(x))

    out_breakpoints <- seq.int(2L, by=2L, length.out=length(x))
    out_partitioning <- PartitioningByEnd(out_breakpoints, names=names(x))
    out <- relist(all.regions, out_partitioning)
    
    mcols(out) <- mcols(x)
    metadata(out) <- metadata(x)
    return(out)
}

setAs("GInteractions", "GRangesList", function(from) .flipGI(from))

setAs("GInteractions", "Pairs", function(from) {
      Pairs(anchors(from, type="first"), anchors(from, type="second"),
            names=names(from), mcols(from))
})

setMethod("pairs", "GInteractions", function(x, id=FALSE, as.grlist=FALSE) {
    if (id) {
        return(as(x, "SelfHits"))
    } else if (as.grlist) {
        return(as(x, "GRangesList"))
    } else {
        return(as(x, "Pairs"))
    }
})

setMethod("pairs", "InteractionSet", function(x, id=FALSE, as.grlist=FALSE) {
    pairs(interactions(x), id=id, as.grlist=as.grlist)
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

makeGInteractionsFromGRangesPairs <- function(x) {
    if (!is(x, "Pairs")) { 
        stop("'x' must be a Pairs object")
    }
    if (!is(first(x), "GRanges") || !is(second(x), "GRanges")) {
        stop("both paired elements must be GRanges")
    }
    out <- GInteractions(anchor1=first(x), anchor2=second(x))
    mcols(out) <- mcols(x)
    names(out) <- names(x)
    return(out)
}

