setGeneric("pairdist", function(x, ...) { standardGeneric("pairdist") })
setMethod("pairdist", "InteractionSet", function(x, type="mid") 
# Outputs an integer vector specifying the distance between the interacting bins,
# depending on the type of distance specified.
{
    type <- match.arg(type, c("mid", "gap", "span", "diag", "intra")) 
    chr <- as.character(seqnames(regions(x)))
    ai1 <- anchors(x, type="first", id=TRUE)
    ai2 <- anchors(x, type="second", id=TRUE)

    # Protection when all inter's.
    is.same <- chr[ai1]==chr[ai2]
    if (type=="intra") { return(is.same) }
    output <- rep(as.integer(NA), nrow(x))
    if (!any(is.same)) { return(output) }

    st <- start(regions(x))
    en <- end(regions(x))
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
})

setGeneric("intrachr", function(x) { standardGeneric("intrachr") })
setMethod("intrachr", "InteractionSet", function(x) { pairdist(x, type="intra") })

