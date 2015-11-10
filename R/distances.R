setGeneric("pairdist", function(x, ...) { standardGeneric("pairdist") })

.get_dist_output <- function(regs, ai1, ai2, type) {
    type <- match.arg(type, c("mid", "gap", "span", "diag", "intra"))
    chr <- as.character(seqnames(regs))

    # Protection when all inter's.
    is.same <- chr[ai1]==chr[ai2]
    if (type=="intra") { return(is.same) }
    output <- rep(as.integer(NA), length(ai1))
    if (!any(is.same)) { return(output) }

    st <- start(regs)
    en <- end(regs)
    ai1 <- ai1[is.same]
    ai2 <- ai2[is.same]
    all.as <- st[ai1]
    all.ae <- en[ai1]
    all.ts <- st[ai2]
    all.te <- en[ai2]

    if (type=="gap") {
        output[is.same] <- pmax(all.as, all.ts) - pmin(all.ae, all.te) - 1L
    } else if (type=="span") {
        output[is.same] <- pmax(all.ae, all.te) - pmin(all.as, all.ts) + 1L
    } else if (type=="mid") {
        output[is.same] <- as.integer((all.as + all.ae - all.ts - all.te)/2L)
    } else if (type=="diag") {
        output[is.same] <- ai1 - ai2
    }
    return(output)
}

setMethod("pairdist", "InteractionSet", function(x, type="mid") 
# Outputs an integer vector specifying the distance between the interacting bins,
# depending on the type of distance specified.
{
    ai1 <- anchors(x, type="first", id=TRUE)
    ai2 <- anchors(x, type="second", id=TRUE)
    .get_dist_output(regions(x), ai1, ai2, type)
})

setMethod("pairdist", "ContactMatrix", function(x, type="mid") 
# Outputs an integer vector specifying the distance between the interacting bins,
# depending on the type of distance specified.
{
    ai1 <- rep(anchors(x, type="row", id=TRUE), ncol(x))
    ai2 <- rep(anchors(x, type="column", id=TRUE), each=nrow(x))
    swapped <- .enforce_order(ai1, ai2) # To get sensible distances
    out <- .get_dist_output(regions(x), swapped$anchor1, swapped$anchor2, type)
    dim(out) <- dim(x@matrix)
    return(out)
})

setGeneric("intrachr", function(x) { standardGeneric("intrachr") })
setMethod("intrachr", "InteractionSet", function(x) { pairdist(x, type="intra") })
setMethod("intrachr", "ContactMatrix", function(x) { pairdist(x, type="intra") })

