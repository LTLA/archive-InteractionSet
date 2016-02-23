##############################################
# Defines the ContactMatrix class.

setClass("ContactMatrix",
    contains="Annotated", 
    slots=list(
        matrix="Matrix", 
        anchor1="integer",
        anchor2="integer",
        regions="GRanges"
    )		
)

setValidity2("ContactMatrix", function(object) {
    if (is.unsorted(object@regions)) {
        return("'regions' should be sorted")
    }
    msg <- .check_inputs(object@anchor1, object@anchor2, object@regions, same.length=FALSE)
    if (is.character(msg)) { return(msg) }

    if (nrow(object@matrix)!=length(object@anchor1)) { 
        return("'matrix' nrow must be equal to length of 'anchor1'")
    }
    if (ncol(object@matrix)!=length(object@anchor2)) {
        return("'matrix' ncol must be equal to length of 'anchor2'")
    }
    return(TRUE)
}) 

setMethod("show", signature("ContactMatrix"), function(object) {
    cat("class:", class(object), "\n")
    cat("dim:", dim(object@matrix), "\n")
    cat("type:", class(object@matrix), "\n")

    rnames <- rownames(object)
    if (!is.null(rnames)) scat("rownames(%d): %s\n", rnames)
    else scat("rownames: NULL\n")

    cnames <- colnames(object)
    if (!is.null(cnames)) scat("colnames(%d): %s\n", cnames)
    else cat("colnames: NULL\n")

    expt <- names(metadata(object))
    if (is.null(expt))
        expt <- character(length(metadata(object)))
    scat("metadata(%d): %s\n", expt)

    cat(sprintf("regions: %i\n", length(object@regions)))
})

##############################################
# Constructor:

.new_ContactMatrix <- function(matrix, anchor1, anchor2, regions, metadata) {
    anchor1 <- as.integer(anchor1)
    anchor2 <- as.integer(anchor2)
    
    msg <- .check_inputs(anchor1, anchor2, regions, same.length=FALSE)
    if (is.character(msg)) { stop(msg) }
    out <- .resort_regions(anchor1, anchor2, regions)

    if (!is(matrix, "Matrix")) { 
        matrix <- Matrix(matrix)
    }
    new("ContactMatrix", matrix=matrix, anchor1=out$anchor1, anchor2=out$anchor2, 
        regions=out$regions, metadata=metadata)
}

setGeneric("ContactMatrix", function(matrix, anchor1, anchor2, regions, ...) standardGeneric("ContactMatrix"))
setMethod("ContactMatrix", c("ANY", "numeric", "numeric", "GRanges"), 
    function(matrix, anchor1, anchor2, regions, metadata=list()) { 
        .new_ContactMatrix(matrix, anchor1, anchor2, regions, metadata)
    }
)

setMethod("ContactMatrix", c("ANY", "GRanges", "GRanges", "GenomicRangesORmissing"), 
    function(matrix, anchor1, anchor2, regions, metadata=list()) { 

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
        
        .new_ContactMatrix(matrix, anchor1, anchor2, regions, metadata)
    }
)

setMethod("ContactMatrix", c("missing", "missing", "missing", "GenomicRangesORmissing"),
    function(matrix, anchor1, anchor2, regions, metadata=list()) {
        if (missing(regions)) { regions <- GRanges() }
        .new_ContactMatrix(Matrix(0L, 0, 0), integer(0), integer(0), regions, metadata)
    } 
)

##############################################
# Subsetting

setMethod("[", c("ContactMatrix", "ANY", "ANY"), function(x, i, j, ..., drop=TRUE) {
    if (!missing(i)) { 
        x@anchor1 <- x@anchor1[i]
        x@matrix <- x@matrix[i,,drop=FALSE]
    }
    if (!missing(j)) {
        x@anchor2 <- x@anchor2[j]
        x@matrix <- x@matrix[,j,drop=FALSE]
    }
    return(x)
}) 

setMethod("[<-", c("ContactMatrix", "ANY", "ANY", "ContactMatrix"), function(x, i, j, ..., value) {
    if (!identical(regions(value), regions(x))) { 
        stop("replacement and original 'regions' must be identical")
    }
    if (!missing(i) && !missing(j)) {
        if (!identical(x@anchor1[i], value@anchor1)) {
            stop("cannot modify row indices for a subset of columns")
        }
        if (!identical(x@anchor2[j], value@anchor2)) {
            stop("cannot modify column indices for a subset of rows")
        }
        x@matrix[i,j] <- value@matrix
    } else if (!missing(i)) { 
        x@anchor1[i] <- value@anchor1
        x@matrix[i,] <- value@matrix
    } else if (!missing(j)) { 
        x@anchor2[j] <- value@anchor2
        x@matrix[,j] <- value@matrix
    }
    return(x)
})

setMethod("subset", "ContactMatrix", function(x, i, j) {
    x[i,j]
})

##############################################
# Combining

setMethod("cbind", "ContactMatrix", function(..., deparse.level=1) {
    incoming <- list(...)
    ref <- incoming[[1]]
    all.regions <- lapply(incoming, FUN=regions)
    all.anchor1 <- lapply(incoming, FUN=slot, name="anchor1")
    all.anchor2 <- lapply(incoming, FUN=slot, name="anchor2")

    # Checking if regions are the same; collating if not.
    unified <- .coerce_to_union(all.regions, all.anchor1, all.anchor2)
    ref@regions <- unified$region
    ref@anchor1 <- unified$anchor1[[1]]

    for (x in unified$anchor1[-1]) {
        if (!identical(ref@anchor1, x)) {
            stop("row anchor indices must be identical for 'cbind'")
        }    
    }

    ref@matrix <- do.call(cbind, lapply(incoming, as.matrix))
    ref@anchor2 <- unlist(unified$anchor2)
    return(ref)
})

setMethod("rbind", "ContactMatrix", function(..., deparse.level=1) {
    incoming <- list(...)
    ref <- incoming[[1]]
    all.regions <- lapply(incoming, FUN=regions)
    all.anchor1 <- lapply(incoming, FUN=slot, name="anchor1")
    all.anchor2 <- lapply(incoming, FUN=slot, name="anchor2")

    # Checking if regions are the same; collating if not.
    unified <- .coerce_to_union(all.regions, all.anchor1, all.anchor2)
    ref@regions <- unified$region
    ref@anchor2 <- unified$anchor2[[1]]

    for (x in unified$anchor2[-1]) { 
        if (!identical(ref@anchor2, x)) { 
            stop("column anchor indices must be identical for 'rbind'")
        }    
    }
    
    ref@matrix <- do.call(rbind, lapply(incoming, as.matrix))
    ref@anchor1 <- unlist(unified$anchor1)
    return(ref)
})

setMethod("t", "ContactMatrix", function(x) { 
    x@matrix <- t(x@matrix)
    tmp <- x@anchor1
    x@anchor1 <- x@anchor2
    x@anchor2 <- tmp
    return(x)
})

##############################################
# Sorting and ordering

setMethod("order", "ContactMatrix", function(..., na.last=TRUE, decreasing=FALSE) {
    incoming <- list(...)
    all.rows <- lapply(incoming, anchors, type="row", id=TRUE)
    all.columns <- lapply(incoming, anchors, type="column", id=TRUE)
    list(row=do.call(order, c(all.rows, na.last=na.last, decreasing=decreasing)),
         column=do.call(order, c(all.columns, na.last=na.last, decreasing=decreasing)))
})

setMethod("sort", "ContactMatrix", function(x, decreasing=FALSE, ...) {
    out <- order(x, decreasing=decreasing)
    x[out$row, out$column]
})

setMethod("duplicated", "ContactMatrix", function(x, incomparables=FALSE, ...) {
    r1 <- duplicated(x@anchor1, incomparables=incomparables, ...)
    r2 <- duplicated(x@anchor2, incomparables=incomparables, ...)
    return(list(row=r1, column=r2))
})

setMethod("unique", "ContactMatrix", function(x, incomparables=FALSE, ...) {
    is.dup <- duplicated(x, incomparables=incomparables, ...)
    return(x[!is.dup$row,!is.dup$column])
})

##############################################
# End

