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

setValidity2("InteractionSet", function(object) {
    if (nrow(object@assays)!=length(object@interactions)) {
        return("'assays' nrow differs from length of anchor vectors")
    } 
    msg <- validObject(object@interactions)
    if (is.character(msg)) { 
        return(msg)
    }
    return(TRUE)
})

setMethod("parallelSlotNames", "InteractionSet", function(x) {
    c("interactions", callNextMethod()) 
})

setMethod("show", signature("InteractionSet"), function(object) {
    callNextMethod()
    cat(sprintf("regions: %i\n", length(regions(object@interactions))))
})

###############################################################
# Constructors

.new_InteractionSet <- function(assays, interactions, colData, ...) {
    # Avoid need to specify column names externally (but respecting it if it is).
    no.names <- FALSE
    if (missing(colData)) {
        assays <- as(Assays(assays), "SimpleList")
        if (length(assays)) {
            nms <- colnames(assays[[1]]) 
            if (is.null(nms)) { 
                no.names <- TRUE
                colnames(assays[[1]]) <- seq_len(ncol(assays[[1]]))
            }
        }
    }
    se0 <- SummarizedExperiment(assays, colData=colData, ...)
    if (no.names) { colnames(se0) <- NULL }
    new("InteractionSet", se0, interactions=interactions)
}

setGeneric("InteractionSet", function(assays, interactions, ...) { standardGeneric("InteractionSet") })
setMethod("InteractionSet", c("ANY", "GInteractions"), function(assays, interactions, ...) { 
        .new_InteractionSet(assays, interactions, ...)
   }
)

setMethod("InteractionSet", c("missing", "missing"), function(assays, interactions, ...) {
        .new_InteractionSet(list(), GInteractions(), ...)
   }
)

###############################################################
# Subsetting

# Need to define these because SummarizedExperiment0 doesn't use extract/replaceROWS directly;
# they divert to these functions anyway.
setMethod("[", c("InteractionSet", "ANY", "ANY"), function(x, i, j, ..., drop=TRUE) {
    if (!missing(i)) { x@interactions <- x@interactions[i] }
    callNextMethod()
})

setMethod("[<-", c("InteractionSet", "ANY", "ANY", "InteractionSet"), function(x, i, j, ..., value) {
    if (!missing(i)) { x@interactions[i] <- value@interactions }
    callNextMethod(x=x, i=i, j=j, ..., value=value)
})

setMethod("subset", "InteractionSet", function(x, i, j) {
    x[i, j]
})

###############################################################
# Combining

setMethod("cbind", "InteractionSet", function(..., deparse.level=1) {
    args <- unname(list(...))
    ans <- args[[1]]
    for (x in args[-1]) {
        if (!identical(interactions(x), interactions(ans))) { 
            # Possible to cbind for different metadata here, but I doubt that gets much use.
            stop("interactions must be identical in 'cbind'") 
        }
    }
    
    base <- do.call(cbind, lapply(args, function(x) { as(x, "SummarizedExperiment0") }))
    new("InteractionSet", base, interactions=interactions(ans))
})

setMethod("rbind", "InteractionSet", function(..., deparse.level=1) {
    args <- unname(list(...))
    base <- do.call(rbind, lapply(args, function(x) { as(x, "SummarizedExperiment0") }))
    new("InteractionSet", base, interactions=do.call(rbind, lapply(args, FUN=interactions)))
})

setMethod("c", "InteractionSet", function(x, ..., recursive = FALSE) {
    rbind(x, ...)
})

# "combine" is slightly different from "rbind", in that it allows different regions to be combined.
setMethod("combine", c("InteractionSet", "InteractionSet"), function(x, y, ...) {
    args <- list(x, y, ...)
    base <- do.call(rbind, lapply(args, function(x) { as(x, "SummarizedExperiment0") }))
    new("InteractionSet", base, interactions=do.call(combine, lapply(args, FUN=interactions)))
})

###############################################################
# Other methods

setMethod("order", "InteractionSet", function(..., na.last=TRUE, decreasing=FALSE) {
    all.ids <- unlist(lapply(list(...), anchors, id=TRUE), recursive=FALSE)
    do.call(order, c(all.ids, list(na.last=na.last, decreasing=decreasing)))
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

###############################################################
# End
