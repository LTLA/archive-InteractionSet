##############################################
# Defines the ContactMatrix class.

setClass("ContactMatrix", 
    slots=list(
        matrix="matrix", 
        anchor1="integer",
        anchor2="integer",
        regions="GRanges"
    )		
)

setValidity("ContactMatrix", function(object) {
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
    cat(sprintf("regions: %i\n", length(object@regions)))
})

##############################################
# Constructor:

.new_ContactMatrix <- function(matrix, anchor1, anchor2, regions) {
    anchor1 <- as.integer(anchor1)
    anchor2 <- as.integer(anchor2)
    
    msg <- .check_inputs(anchor1, anchor2, regions, same.length=FALSE)
    if (is.character(msg)) { stop(msg) }
    out <- .resort_regions(anchor1, anchor2, regions, enforce.order=FALSE)

    new("ContactMatrix", matrix=matrix, anchor1=out$anchor1, anchor2=out$anchor2, regions=out$regions)
}

setGeneric("ContactMatrix", function(matrix, anchor1, anchor2, ...) { standardGeneric("ContactMatrix") })
setMethod("ContactMatrix", c("ANY", "numeric", "numeric"), .new_ContactMatrix)

setMethod("ContactMatrix", c("ANY", "GRanges", "GRanges"), 
    function(matrix, anchor1, anchor2, regions) { 

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
        
        .new_ContactMatrix(matrix, anchor1, anchor2, regions)
    }
)

##############################################
# Matrix dimensions

setMethod("dim", "ContactMatrix", function(x) { 
    dim(x@matrix)
})

setMethod("length", "ContactMatrix", function(x) { 
    length(x@matrix)
})

setMethod("dimnames", "ContactMatrix", function(x) {
    dimnames(x@matrix)
})

setReplaceMethod("dimnames", "ContactMatrix", function(x, value) {
    dimnames(x@matrix) <- value
    return(x)
})

setMethod("as.matrix", "ContactMatrix", function(x) {
    return(x@matrix)
}) 

##############################################
# Subsetting

setMethod("[", c("ContactMatrix", "ANY", "ANY"), function(x, i, j, ..., drop=TRUE) {
    if (!missing(i)) { 
        x@anchor1 <- x@anchor1[i]
    }
    if (!missing(j)) {
        x@anchor2 <- x@anchor2[j]
    }
    x@matrix <- x@matrix[i,j,drop=FALSE]
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
    for (x in incoming[-1]) {
        if (!identical(regions(ref), regions(x))) { 
            stop("'regions' must be identical for 'cbind'")
        }
        if (!identical(anchors(ref, type="row", id=TRUE),
                       anchors(x, type="row", id=TRUE))) {
            stop("row anchor indices must be identical for 'cbind'")
        }    
    }
    
    ref@matrix <- do.call(cbind, lapply(incoming, as.matrix))
    ref@anchor2 <- unlist(lapply(incoming, anchors, id=TRUE, type="column"))
    return(ref)
})

setMethod("rbind", "ContactMatrix", function(..., deparse.level=1) {
    incoming <- list(...)
    ref <- incoming[[1]]
    for (x in incoming[-1]) {
        if (!identical(regions(ref), regions(x))) { 
            stop("'regions' must be identical for 'rbind'")
        }
        if (!identical(anchors(ref, type="column", id=TRUE),
                       anchors(x, type="column", id=TRUE))) {
            stop("column anchor indices must be identical for 'rbind'")
        }    
    }
    
    ref@matrix <- do.call(rbind, lapply(incoming, as.matrix))
    ref@anchor1 <- unlist(lapply(incoming, anchors, id=TRUE, type="row"))
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
# overlapsAny

setMethod("overlapsAny", c("ContactMatrix", "GRanges"), 
    function(query, subject, type=c("any", "start", "end", "within", "equal"),
        algorithm=c("nclist", "intervaltree"), ignore.strand=TRUE) {
        a1 <- anchors(query, id=TRUE, type="row")
        a2 <- anchors(query, id=TRUE, type="column")
        
        is.used <- union(a1, a2)
        is.overlapped <- logical(length(regions(query)))
        is.overlapped[is.used] <- overlapsAny(regions(query)[is.used], subject, type=type,
                                       algorithm=algorithm, ignore.strand=ignore.strand)
        return(list(row=is.overlapped[a1], column=is.overlapped[a2]))
})

# Use outer(output$row, output$column, "|" or "&") to get the logical area in the interaction space.
# Not sure it makes a great deal of sense to define 'findOverlaps' here.

##############################################
# End


# Testing:
# cm <- ContactMatrix(matrix(runif(16), 4, 4), sample(4), sample(4), GRanges("chr1", IRanges(1:4, 1:4)))

