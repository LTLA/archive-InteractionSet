##############################################
# Setting validity and show methods.

setValidity2("ContactMatrix", function(object) {
    if (is.unsorted(regions(object))) {
        return("'regions' should be sorted")
    }
    msg <- .check_inputs(anchor1(object), anchor2(object), regions(object), same.length=FALSE)
    if (is.character(msg)) { return(msg) }

    mat <- as.matrix(object)
    if (length(dim(mat))!=2L) {
        return("'matrix' slot should contain a matrix-like object with 2 dimensions")
    }
    if (nrow(mat)!=length(anchor1(object))) { 
        return("'matrix' nrow must be equal to length of 'anchor1'")
    }
    if (ncol(mat)!=length(anchor2(object))) {
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
        .new_ContactMatrix(base::matrix(0L, 0, 0), integer(0), integer(0), regions, metadata)
    } 
)

##############################################
# Subsetting

setMethod("[", c("ContactMatrix", "ANY", "ANY"), function(x, i, j, ..., drop=TRUE) {
    if (!missing(i)) { 
        unchecked_anchor1(x) <- anchor1(x)[i]
        unchecked_matrix(x) <- as.matrix(x)[i,,drop=FALSE]
    }
    if (!missing(j)) {
        unchecked_anchor2(x) <- anchor2(x)[j]
        unchecked_matrix(x) <- as.matrix(x)[,j,drop=FALSE]
    }
    return(x)
}) 

setMethod("[<-", c("ContactMatrix", "ANY", "ANY", "ContactMatrix"), function(x, i, j, ..., value) {
    if (!identical(regions(value), regions(x))) { 
        stop("replacement and original 'regions' must be identical")
    }
    
    a1 <- anchor1(x)
    a2 <- anchor2(x)
    m <- as.matrix(x)

    if (!missing(i) && !missing(j)) {
        if (!identical(anchor1(x)[i], anchor1(value))) {
            stop("cannot modify row indices for a subset of columns")
        }
        if (!identical(anchor2(x)[j], anchor2(value))) {
            stop("cannot modify column indices for a subset of rows")
        }
        m[i,j] <- as.matrix(value)
    } else if (!missing(i)) { 
        a1[i] <- anchor1(value)
        m[i,] <- as.matrix(value)
    } else if (!missing(j)) { 
        a2[j] <- anchor2(value)
        m[,j] <- as.matrix(value)
    }

    unchecked_anchor1(x) <- a1
    unchecked_anchor2(x) <- a2
    unchecked_matrix(x) <- m
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
    unchecked_regions(ref) <- unified$region
    unchecked_anchor1(ref) <- unified$anchor1[[1]]

    for (x in unified$anchor1[-1]) {
        if (!identical(anchor1(ref), x)) {
            stop("row anchor indices must be identical for 'cbind'")
        }    
    }

    unchecked_matrix(ref) <- do.call(cbind, lapply(incoming, as.matrix))
    unchecked_anchor2(ref) <- unlist(unified$anchor2)
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
    unchecked_regions(ref) <- unified$region
    unchecked_anchor2(ref) <- unified$anchor2[[1]]

    for (x in unified$anchor2[-1]) { 
        if (!identical(anchor2(ref), x)) { 
            stop("column anchor indices must be identical for 'rbind'")
        }    
    }
    
    unchecked_matrix(ref) <- do.call(rbind, lapply(incoming, as.matrix))
    unchecked_anchor1(ref) <- unlist(unified$anchor1)
    return(ref)
})

setMethod("t", "ContactMatrix", function(x) { 
    unchecked_matrix(x) <- t(as.matrix(x))
    tmp <- anchor1(x)
    unchecked_anchor1(x) <- anchor2(x)
    unchecked_anchor2(x) <- tmp
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

