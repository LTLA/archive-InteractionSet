.get_dist_output <- function(regs, ai1, ai2, type) {
    type <- match.arg(type, c("mid", "gap", "span", "diag", "intra"))
    chr <- as.character(seqnames(regs))

    # To get sensible distances
    swapped <- .enforce_order(ai1, ai2) 
    larger <- swapped[[2]]
    smaller <- swapped[[1]]

    # Protection when all inter's.
    is.same <- chr[larger]==chr[smaller]
    if (type=="intra") { return(is.same) }
    output <- rep(as.integer(NA), length(larger))
    if (!any(is.same)) { return(output) }

    st <- start(regs)
    en <- end(regs)
    larger <- larger[is.same]
    smaller <- smaller[is.same]
    all.as <- st[larger]
    all.ae <- en[larger]
    all.ts <- st[smaller]
    all.te <- en[smaller]

    if (type=="gap") {
        output[is.same] <- pmax(all.as, all.ts) - pmin(all.ae, all.te) - 1L
    } else if (type=="span") {
        output[is.same] <- pmax(all.ae, all.te) - pmin(all.as, all.ts) + 1L
    } else if (type=="mid") {
        output[is.same] <- as.integer(abs(all.as + all.ae - all.ts - all.te)/2L) # Need 'abs', in case later range has earlier midpoint (e.g., if nested).
    } else if (type=="diag") {
        output[is.same] <- larger - smaller
    }
    return(output)
}

setMethod("pairdist", "GInteractions", function(x, type="mid") 
# Outputs an integer vector specifying the distance between the interacting bins,
# depending on the type of distance specified.
{
    ai1 <- anchors(x, type="first", id=TRUE)
    ai2 <- anchors(x, type="second", id=TRUE)
    .get_dist_output(regions(x), ai1, ai2, type)
})

setMethod("pairdist", "InteractionSet", function(x, type="mid") { 
    pairdist(interactions(x), type=type) 
})

setMethod("pairdist", "ContactMatrix", function(x, type="mid") 
# Outputs an integer vector specifying the distance between the interacting bins,
# depending on the type of distance specified.
{
    ai1 <- rep(anchors(x, type="row", id=TRUE), ncol(x))
    ai2 <- rep(anchors(x, type="column", id=TRUE), each=nrow(x))
    out <- .get_dist_output(regions(x), ai1, ai2, type)
    dim(out) <- dim(as.matrix(x))
    return(out)
})

setMethod("intrachr", "GInteractions", function(x) { pairdist(x, type="intra") })
setMethod("intrachr", "InteractionSet", function(x) { pairdist(x, type="intra") })
setMethod("intrachr", "ContactMatrix", function(x) { pairdist(x, type="intra") })

