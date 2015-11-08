###############################################################
# Defines the InteractionSet class, based on the SummarizedExperiment base class.
# This allows us to avoid re-defining various standard functions.

setClass("InteractionSet", 
     contains="SummarizedExperiment0",
     representation(
         anchor1="integer",
         anchor2="integer",
         regions="GRanges"
     ),
     prototype(
         anchor1=integer(0),
         anchor2=integer(0),
         regions=GRanges()
     )
)

setValidity("InteractionSet", function(object) {
    if (nrow(object@assays)!=length(object@anchor1)) {
        return("'assay' nrow differs from length of first anchor vector")
    } 
    if (nrow(object@assays)!=length(object@anchor2)) { 
        return("'assay' nrow differs from length of second anchor vector")
    }
    if (ncol(object@assays)!=nrow(object@colData)) {
        return("'assay' ncol differs from 'colData' nrow")
    }

    if (!all(object@anchor1 >= 1L) || !all(object@anchor2 >= 1L)) {
        return('all anchors must be positive integers')
    } 
    nregs <- length(object@regions)
    if ( !all(object@anchor1 <= nregs) || !all(object@anchor2 <= nregs)
        || !all(is.finite(object@anchor1)) || !all(is.finite(object@anchor2)) ) {
        return('all anchors must refer to valid regions')
    } 
    if (!all(object@anchor1 >= object@anchor2)) { 
        return('first anchors cannot be less than the second anchor')
    }
    return(TRUE)
})

setMethod("show", signature("InteractionSet"), function(object) {
    callNextMethod()
    cat(sprintf("regions: %i\n", length(object@regions)))
})

###############################################################
# Constructors

.enforce_order <- function(anchor1, anchor2) {
    swap <- anchor2 > anchor1
    if (any(swap)) { 
        temp <- anchor1[swap]
        anchor1[swap] <- anchor2[swap]
        anchor2[swap] <- temp
    }
    return(list(anchor1=anchor1, anchor2=anchor2))   
}

.resort_regions <- function(anchor1, anchor2, regions) {
    o <- order(regions)
    if (any(diff(o)!=1L)) { 
        new.pos <- seq_along(o)
        new.pos[o] <- new.pos
        anchor1 <- new.pos[anchor1]
        anchor2 <- new.pos[anchor2]
        regions <- regions[o]
    } 
    out <- .enforce_order(anchor1, anchor2)
    return(list(anchor1=out$anchor1, anchor2=out$anchor2, regions=regions)) 
}

.new_InteractionSet <- function(assays, anchor1, anchor2, regions, colData, metadata) {
    elementMetadata <- new("DataFrame", nrows=length(anchor1))
    if (!is(assays, "Assays")) { 
        assays <- Assays(assays)
    }
    if (ncol(colData)==0L) { colData <- new("DataFrame", nrows=ncol(assays)) } # If empty, we fill it.

    # Checking odds and ends.
    anchor1 <- as.integer(anchor1)
    anchor2 <- as.integer(anchor2)
    out <- .resort_regions(anchor1, anchor2, regions)
    anchor1 <- out$anchor1
    anchor2 <- out$anchor2
    regions <- out$regions

    new("InteractionSet", 
        anchor1=anchor1,
        anchor2=anchor2,
        regions=regions,
        colData=colData,
        assays=assays,
        elementMetadata=elementMetadata,
        metadata=as.list(metadata))
}

setGeneric("InteractionSet", function(assays, anchor1, anchor2, ...) { standardGeneric("InteractionSet") })
setMethod("InteractionSet", c("ANY", "numeric", "numeric"),
   function(assays, anchor1, anchor2, regions, colData=DataFrame(), metadata=list()) {
       .new_InteractionSet(assays, anchor1=anchor1, anchor2=anchor2, 
            regions=regions, colData=colData, metadata=metadata)
   }
)

.collate_GRanges <- function(...) {
    incoming <- list(...)
    obj.dex <- rep(seq_along(incoming), lengths(incoming))
    combined <- do.call(c, incoming)
    refdex <- seq_along(combined)
    
    # Sorting and re-indexing.
    o <- order(combined)
    new.pos <- seq_along(combined)
    new.pos[o] <- new.pos
    refdex <- new.pos[refdex]
    combined <- combined[o]

    # Removing duplicates and re-indexing.
    is.first <- !duplicated(combined)
    new.pos <- cumsum(is.first)
    combined <- combined[is.first]
    refdex <- new.pos[refdex]    
    return(list(indices=split(refdex, obj.dex), ranges=combined))
}

setMethod("InteractionSet", c("ANY", "GRanges", "GRanges"), 
   function(assays, anchor1, anchor2, regions, colData=DataFrame(), metadata=list()) {

        # Making unique regions to save space (metadata is ignored)
        if (missing(regions)) {
            collated <- .collate_GRanges(anchor1, anchor2)
            regions <- collated$ranges
            anchor1 <- collated$indices[[1]]
            anchor2 <- collated$indices[[2]]
        } else {
            anchor1 <- match(anchor1, regions)
            anchor2 <- match(anchor2, regions)
            if (any(is.na(anchor1)) || any(is.na(anchor2))) {
                stop("anchor regions missing in specified 'regions'")
            }
        }

       .new_InteractionSet(assays, anchor1=anchor1, anchor2=anchor2, 
            regions=regions, colData=colData, metadata=metadata)
   }
)

###############################################################
# Getters 

setGeneric("anchors", function(x, ...) { standardGeneric("anchors") })
setMethod("anchors", signature("InteractionSet"), 
    function(x, type="both", id=FALSE) {
        type <- match.arg(type, c("both", "first", "second"))
        if (id) { 
            if (type=="both") {
                return(list(first=x@anchor1, second=x@anchor2))
            } else if (type=="first") { 
                return(x@anchor1)
            } else { 
                return(x@anchor2)
            }
        } else {
            if (type=="both") { 
                return(GRangesList("first"=x@regions[x@anchor1],
                                   "second"=x@regions[x@anchor2]))
            } else if (type=="first") { 
                return(x@regions[x@anchor1])
            } else {
                return(x@regions[x@anchor2])
            }
        }
    }
)

setGeneric("regions", function(x, ...) { standardGeneric("regions") })
setMethod("regions", signature("InteractionSet"), function(x) {
    x@regions
})

###############################################################
# Setters

setGeneric("regions<-", function(x, value) { standardGeneric("regions<-") }) 
setReplaceMethod("regions", "InteractionSet", function(x, value) {
    out <- .resort_regions(x@anchor1, x@anchor2, value)
    x@anchor1 <- out$anchor1 
    x@anchor2 <- out$anchor2
    x@regions <- out$regions
    validObject(x) 
    return(x)
})

setGeneric("anchors<-", function(x, ..., value) { standardGeneric("anchors<-") })
setReplaceMethod("anchors", "InteractionSet", function(x, value) {
    out <- .enforce_order(as.integer(value[[1]]), as.integer(value[[2]]))
    x@anchor1 <- out$anchor1
    x@anchor2 <- out$anchor2
    validObject(x)
    return(x)
})

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
