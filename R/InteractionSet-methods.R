###############################################################

setValidity2("InteractionSet", function(object) {
    if (nrow(object)!=length(interactions(object))) {
        return("'interactions' length is not equal to the number of rows")
    } 
    if (!is.null(object@NAMES)) { # Using direct slot access, otherwise diverts to slots in GInteractions.
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
    cat(sprintf("type: %s\n", class(interactions(object))))
    cat(sprintf("regions: %i\n", length(regions(interactions(object)))))
})

###############################################################
# Constructors

.new_InteractionSet <- function(assays, interactions, ...) {
    se0 <- SummarizedExperiment(assays, ...)
    new("InteractionSet", se0, interactions=interactions)
}

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
    if (!missing(i)) { unchecked_interactions(x) <- interactions(x)[i] }
    callNextMethod()
})

setMethod("[<-", c("InteractionSet", "ANY", "ANY", "InteractionSet"), function(x, i, j, ..., value) {
    if (!missing(i)) { 
        I <- interactions(x)
        I[i] <- interactions(value) 
        unchecked_interactions(x) <- I
    }
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
