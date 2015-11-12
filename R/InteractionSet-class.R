###############################################################
# Defines the InteractionSet class, based on the SummarizedExperiment base class.
# This allows us to avoid re-defining various standard functions.

setClass("InteractionSet", 
    contains="SummarizedExperiment0",
    representation(
        interactions="GInteractions"
    ),
    prototype(
        interactions=GInteractions()
    )
)

setValidity("InteractionSet", function(object) {
    if (nrow(object@assays)!=length(object@interactions)) {
        return("'assay' nrow differs from length of anchor vectors")
    } 
    msg <- validObject(object@interactions)
    if (is.character(msg)) { 
        return(msg)
    }
    return(TRUE)
})

setMethod("show", signature("InteractionSet"), function(object) {
    callNextMethod()
    cat(sprintf("regions: %i\n", length(object@regions)))
})

###############################################################
# Constructors

.new_InteractionSet <- function(assays, interactions, colData, metadata) {
    elementMetadata <- new("DataFrame", nrows=length(interactions))
    if (!is(assays, "Assays")) { 
        assays <- Assays(assays)
    }
    if (ncol(colData)==0L) { colData <- new("DataFrame", nrows=ncol(assays)) } # If empty, we fill it.

    new("InteractionSet", 
        interactions=interactions,
        colData=colData,
        assays=assays,
        elementMetadata=elementMetadata,
        metadata=as.list(metadata))
}

setGeneric("InteractionSet", function(assays, interactions, ...) { standardGeneric("InteractionSet") })
setMethod("InteractionSet", c("ANY", "GInteractions"),
    function(assays, interactions, colData=DataFrame(), metadata=list()) {
        .new_InteractionSet(assays, interactions, colData=colData, metadata=metadata)
   }
)

setMethod("InteractionSet", c("missing", "missing"),
    function(assays, interactions, colData=DataFrame(), metadata=list()) {
        .new_InteractionSet(list(), GInteractions(), colData=colData, metadata=metadata)
   }
)

###############################################################
# Subsetting

setMethod("[", c("InteractionSet", "ANY", "ANY"), function(x, i, j, ..., drop=TRUE) {
    a1 <- x@anchor1[i]
    a2 <- x@anchor2[i]
    ans <- callNextMethod()
    ans@anchor1 <- a1
    ans@anchor2 <- a2
    return(ans)
})

setMethod("subset", "InteractionSet", function(x, i, j) {
    x[i, j]
})

###############################################################
# Combining

setMethod("cbind", "InteractionSet", function(..., deparse.level=1) {
    args <- unname(list(...))
    ans <- args[[1]]
    unequal.element.meta <- FALSE

    for (x in args[-1]) {
        if (!identical(regions(x), regions(ans))) { stop("regions must be identical in 'cbind'") }
        if (!identical(anchors(x, id=TRUE), anchors(ans, id=TRUE))) { 
            stop("anchors must be identical in 'cbind'") } 
        if (!identical(elementMetadata(x), elementMetadata(ans))) { unequal.element.meta <- TRUE }
    }
    if (unequal.element.meta) { warning("'elementMetadata' of first object used in 'cbind'") }
    
    # Need to dig into SummarizedExperiment, as callNextMethod fails (undefined dispatch with ...?)
    SummarizedExperiment:::.cbind.SummarizedExperiment(args)
})

setMethod("rbind", "InteractionSet", function(..., deparse.level=1) {
    args <- unname(list(...))
    ans <- args[[1]]
    all1 <- anchors(args[[1]], type="first", id=TRUE)
    all2 <- anchors(args[[1]], type="second", id=TRUE)

    for (x in args[-1]) { 
        if (!identical(regions(x), regions(ans))) { stop("regions must be identical in 'rbind'") }
        all1 <- c(all1, anchors(x, type="first", id=TRUE))
        all2 <- c(all2, anchors(x, type="second", id=TRUE))
    }

    args[[1]]@anchor1 <- all1
    args[[1]]@anchor2 <- all2
    SummarizedExperiment:::.rbind.SummarizedExperiment(args)
})

# "c" is slightly different from "rbind", in that it allows different regions to be combined.
setMethod("c", "InteractionSet", function(x, ..., recursive = FALSE) {
    incoming <- list(x, ...)
    all.regions <- lapply(incoming, FUN=regions)
    collated <- do.call(.collate_GRanges, all.regions)

    for (i in seq_along(incoming)) {
        cur.anchors <- anchors(incoming[[i]], id=TRUE)
        incoming[[i]]@regions <- collated$ranges
        anchors(incoming[[i]]) <- list(collated$indices[[i]][cur.anchors$first],
                                       collated$indices[[i]][cur.anchors$second])
    }

    do.call(rbind, incoming)
})

###############################################################
# Other methods

setMethod("order", "InteractionSet", function(..., na.last=TRUE, decreasing=FALSE) {
    all.ids <- unlist(lapply(list(...), anchors, id=TRUE), recursive=FALSE)
    do.call(order, c(all.ids, list(na.last=na.last, decreasing=decreasing)))
})

setMethod("sort", "InteractionSet", function(x, decreasing=FALSE, ...) {
    x[order(x, decreasing=decreasing),]
})

setMethod("duplicated", "InteractionSet", function(x, incomparables=FALSE, fromLast=FALSE, ...) 
# Stable sort required here: first entry in 'x' is always non-duplicate if fromLast=FALSE,
# and last entry is non-duplicate if fromLast=TRUE.
{
    if (!nrow(x)) { return(logical(0)) }
    a1 <- anchors(x, id=TRUE, type="first")
    a2 <- anchors(x, id=TRUE, type="second")
    o <- order(a1, a2) 
    if (fromLast) { o <- rev(o) }
    is.dup <- c(FALSE, diff(a1[o])==0L & diff(a2[o])==0L)
    is.dup[o] <- is.dup
    return(is.dup)
})

setMethod("unique", "InteractionSet", function(x, incomparables=FALSE, fromLast=FALSE, ...) {
    x[!duplicated(x, incomparables=incomparables, fromLast=fromLast, ...),]
})

# Not sure how much sense it makes to provide GRanges methods on the InteractionSet,
# as these'll be operating on 'regions' rather than 'anchors'.
#setMethod("seqinfo", "InteractionSet", function(x) {
#     seqinfo(x@regions)
#})
#
#setReplaceMethod("seqinfo", "InteractionSet", function(x, value) {
#    seqinfo(x@regions) <- value
#    validObject(x)
#    return(x)
#})

setMethod("split", "InteractionSet", function(x, f, drop=FALSE, ...) {
    splitAsList(x, f, drop=drop)
})

###############################################################
# End
