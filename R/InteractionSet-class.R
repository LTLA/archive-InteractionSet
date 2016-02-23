###############################################################
# Defines the InteractionSet class, based on the SummarizedExperiment base class.
# This allows us to avoid re-defining various standard functions.

setClass("InteractionSet", 
    contains="SummarizedExperiment",
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
    if (!is.null(object@NAMES)) {
        return("'NAMES' slot must always be NULL")
    }
    if (ncol(object@elementMetadata) != 0L) {
        return("'elementMetadata' slot must always have zero columns")
    }
    return(TRUE)
})

setMethod("parallelSlotNames", "InteractionSet", function(x) {
    c("interactions", callNextMethod()) 
})

setMethod("show", signature("InteractionSet"), function(object) {
    callNextMethod()
    cat(sprintf("type: %s\n", class(object@interactions)))
    cat(sprintf("regions: %i\n", length(regions(object@interactions))))
})

###############################################################
# Constructors

.new_InteractionSet <- function(assays, interactions, ...) {
    se0 <- SummarizedExperiment(assays, ...)
    new("InteractionSet", se0, interactions=interactions)
}

setGeneric("InteractionSet", function(assays, interactions, ...) standardGeneric("InteractionSet"))
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

# Need to define these because SummarizedExperiment doesn't use extract/replaceROWS directly;
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
    ref.inters <- interactions(ans)
    ref.length <- length(ans)

    for (x in args[-1]) {
        if (ref.length!=length(x) || any(interactions(x)!=ref.inters)) { 
            # Possible to cbind for different metadata here, but I doubt that gets much use.
            stop("interactions must be identical in 'cbind'") 
        }
    }
    
    base <- do.call(cbind, lapply(args, function(x) { as(x, "SummarizedExperiment") }))
    new("InteractionSet", base, interactions=interactions(ans))
})

setMethod("rbind", "InteractionSet", function(..., deparse.level=1) {
    args <- unname(list(...))
    base <- do.call(rbind, lapply(args, function(x) { as(x, "SummarizedExperiment") }))
    new("InteractionSet", base, interactions=do.call(rbind, lapply(args, FUN=interactions)))
})

setMethod("c", "InteractionSet", function(x, ..., recursive = FALSE) {
    rbind(x, ...)
})

###############################################################
# Other methods

setMethod("order", "InteractionSet", function(..., na.last=TRUE, decreasing=FALSE) {
    do.call(order, c(lapply(list(...), interactions), list(na.last=na.last, decreasing=decreasing)))
})

setMethod("duplicated", "InteractionSet", function(x, incomparables=FALSE, fromLast=FALSE, ...) {
    duplicated(interactions(x))
})

###############################################################
# End
