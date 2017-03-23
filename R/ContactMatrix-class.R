##############################################
# Setting validity and show methods.

setValidity2("ContactMatrix", function(object) {
    if (is.unsorted(regions(object))) {
        return("'regions' should be sorted")
    }
    msg <- .check_inputs(anchor1(object), anchor2(object), regions(object), same.length=FALSE)
    if (is.character(msg)) { return(msg) }

    if (nrow(as.matrix(object))!=length(anchor1(object))) { 
        return("'matrix' nrow must be equal to length of 'anchor1'")
    }
    if (ncol(as.matrix(object))!=length(anchor2(object))) {
        return("'matrix' ncol must be equal to length of 'anchor2'")
    }
    return(TRUE)
}) 

setMethod("show", signature("ContactMatrix"), function(object) {
    cat("class:", class(object), "\n")
    cat("dim:", dim(as.matrix(object)), "\n")
    cat("type:", class(as.matrix(object)), "\n")

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

    cat(sprintf("regions: %i\n", length(regions(object))))
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
        x@anchor1 <- anchor1(x)[i]
        x@matrix <- as.matrix(x)[i,,drop=FALSE]
    }
    if (!missing(j)) {
        x@anchor2 <- anchor2(x)[j]
        x@matrix <- as.matrix(x)[,j,drop=FALSE]
    }
    return(x)
}) 

setMethod("[<-", c("ContactMatrix", "ANY", "ANY", "ContactMatrix"), function(x, i, j, ..., value) {
    if (!identical(regions(value), regions(x))) { 
        stop("replacement and original 'regions' must be identical")
    }
    if (!missing(i) && !missing(j)) {
        if (!identical(anchor1(x)[i], anchor1(value))) {
            stop("cannot modify row indices for a subset of columns")
        }
        if (!identical(anchor2(x)[j], anchor2(value))) {
            stop("cannot modify column indices for a subset of rows")
        }
        x@matrix[i,j] <- as.matrix(value)
    } else if (!missing(i)) { 
        x@anchor1[i] <- anchor1(value)
        x@matrix[i,] <- as.matrix(value)
    } else if (!missing(j)) { 
        x@anchor2[j] <- anchor2(value)
        x@matrix[,j] <- as.matrix(value)
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
        if (!identical(anchor1(ref), x)) {
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
        if (!identical(anchor2(ref), x)) { 
            stop("column anchor indices must be identical for 'rbind'")
        }    
    }
    
    ref@matrix <- do.call(rbind, lapply(incoming, as.matrix))
    ref@anchor1 <- unlist(unified$anchor1)
    return(ref)
})

setMethod("t", "ContactMatrix", function(x) { 
    x@matrix <- t(as.matrix(x))
    tmp <- anchor1(x)
    x@anchor1 <- anchor2(x)
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
    r1 <- duplicated(anchor1(x), incomparables=incomparables, ...)
    r2 <- duplicated(anchor2(x), incomparables=incomparables, ...)
    return(list(row=r1, column=r2))
})

setMethod("unique", "ContactMatrix", function(x, incomparables=FALSE, ...) {
    is.dup <- duplicated(x, incomparables=incomparables, ...)
    return(x[!is.dup$row,!is.dup$column])
})

##############################################
# End

