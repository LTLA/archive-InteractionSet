.boundingBox <- function(x, f) 
# As the name suggests, it computes a bounding box for groups of interactions in 'x'.
# Groups are specified by a factor 'f'.
#
# written by Aaron Lun
# created 2 December 2015
{
    if (missing(f)) { 
        f <- rep(1L, length(x))
    }
    if (length(f)!=length(x)) { 
        stop("length of 'f' must be equal to number of interactions")
    } 
    o <- order(f)
    f <- f[o]

    a1 <- anchors(x, type="first", id=TRUE)
    a2 <- anchors(x, type="second", id=TRUE)
    a1 <- a1[o] - 1L # Play nice with zero indexing in C++.
    a2 <- a2[o] - 1L
    is.first <- !duplicated(f)
    ref.fac <- as.character(f[is.first])
    f <- cumsum(is.first) - 1L

    chrs <- seqnames(regions(x))
    ref.chr <- levels(chrs)
    chrs <- as.integer(chrs)
    starts <- start(regions(x))
    ends <- end(regions(x))

    # Using indices for easy comparison inside C++.
    bound1 <- .Call(cxx_get_box_bounds, f, ref.fac, a1, chrs, starts, ends)
    if (is.character(bound1)) { stop(bound1) }
    bound2 <- .Call(cxx_get_box_bounds, f, ref.fac, a2, chrs, starts, ends)
    if (is.character(bound2)) { stop(bound2) }

    gr1 <- GRanges(ref.chr[bound1[[2]]], IRanges(bound1[[3]], bound1[[4]]), seqinfo=seqinfo(x)) 
    gr2 <- GRanges(ref.chr[bound2[[2]]], IRanges(bound2[[3]], bound2[[4]]), seqinfo=seqinfo(x))
    out <- GInteractions(gr1, gr2)
    names(out) <- ref.fac[bound1[[1]]+1L]
    return(out)
}

setMethod("boundingBox", "GInteractions", .boundingBox)
setMethod("boundingBox", "InteractionSet", .boundingBox)

