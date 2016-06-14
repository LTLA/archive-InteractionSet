# Inflate from InteractionSet to ContactMatrix.

setGeneric("inflate", function(x, ...) standardGeneric("inflate"))

.make_to_indices <- function(regs, i, ...) {
    nregs <- length(regs)
    if (is.numeric(i)) { 
        i <- as.integer(i)
        if (any(!is.finite(i)) || any(i<=0L) || any(i > nregs)) { 
            stop("indices must be positive integers no greater than 'length(regions(x))'") 
        }
        return(i)
    } else if (is.character(i)) { 
        return(which(seqnames(regs) %in% i))
    } else if (is(i, "GRanges")) {
        return(which(overlapsAny(regs, i, ...)))
    } else {
        stop("invalid value for row/column selection")
    }
}

setMethod("inflate", "GInteractions", function(x, rows, columns, fill, swap=TRUE, sparse=FALSE, ...) {
    row.chosen <- .make_to_indices(regions(x), rows, ...)
    col.chosen <- .make_to_indices(regions(x), columns, ...)
    fill <- rep(fill, length.out=length(x))
     
    # Removing duplicated rows and resorting (we'll put them back in later)
    ro <- order(row.chosen)
    co <- order(col.chosen)
    row.chosen <- row.chosen[ro]
    col.chosen <- col.chosen[co]
    rnd <- !duplicated(row.chosen)
    cnd <- !duplicated(col.chosen)
    row.chosen <- row.chosen[rnd]   
    col.chosen <- col.chosen[cnd]

    # Duplicated interactions can't be handled.
    dx <- duplicated(x)
    if (any(dx)) { 
        warning("duplicated interactions in 'x' are removed")
        x <- x[!dx,]
        fill <- fill[!dx]
    }

    # Matching.
    a1 <- anchors(x, type="first", id=TRUE)
    a2 <- anchors(x, type="second", id=TRUE)
    ar1 <- match(a1, row.chosen)
    ac1 <- match(a1, col.chosen)
    ar2 <- match(a2, row.chosen)
    ac2 <- match(a2, col.chosen)

    # Filling.
    nR <- length(row.chosen)
    nC <- length(col.chosen)
    relevantA <- !is.na(ar1) & !is.na(ac2)
    relevantB <- !is.na(ar2) & !is.na(ac1)

    if (!sparse) { 
        out.mat <- Matrix(as(NA, typeof(fill)), nR, nC)
        if (any(relevantA)) { out.mat[(ac2[relevantA] - 1L) * nR + ar1[relevantA]] <- fill[relevantA] } # Preserve class as dsyMatrix if empty.
    } else {
        out.mat <- sparseMatrix(i=ar1[relevantA], j=ac2[relevantA], 
                                x=fill[relevantA], dims=c(nR, nC))
    }
    if (swap && any(relevantB)) { out.mat[(ac1[relevantB] - 1L) * nR + ar2[relevantB]] <- fill[relevantB] }

    # Restoring the original order.
    original.rows <- cumsum(rnd)
    original.rows[ro] <- original.rows
    original.cols <- cumsum(cnd)
    original.cols[co] <- original.cols

    return(ContactMatrix(out.mat[original.rows,original.cols,drop=FALSE], 
                row.chosen[original.rows], col.chosen[original.cols], regions(x)))
})
 
setMethod("inflate", "InteractionSet", function(x, rows, columns, assay=1, sample=1, fill=NULL, swap=TRUE, sparse=FALSE, ...) {
    if (length(fill)==0L) { fill <- assay(x, assay)[,sample] }
    inflate(interactions(x), rows, columns, fill=fill, swap=swap, sparse=sparse, ...)
})

setGeneric("deflate", function(x, ...) standardGeneric("deflate"))

setMethod("deflate", "ContactMatrix", function(x, collapse=TRUE, extract, use.zero, use.na, ...) {
    # Choosing the expansion strategy.
    if (missing(extract)) { 
        is.sparse <- is(as.matrix(x), "sparseMatrix")
        if (missing(use.zero)) {
            use.zero <- !is.sparse
        }
        if (missing(use.na)) { 
            use.na <- is.sparse
        }
        
        if (use.na && use.zero) { 
            is.valid <- seq_along(as.matrix(x))
        } else {
            is.valid <- TRUE
            if (!use.zero) { 
                is.valid <- is.valid & as.matrix(x)!=0
            } 
            if (!use.na) { 
                is.valid <- is.valid & !is.na(as.matrix(x))
            }
            is.valid <- Matrix::which(is.valid)
        }
    } else {
        if (!identical(length(extract), length(x))) { 
            stop("extraction matrix must be of the same length as 'x'")
        }
        is.valid <- Matrix::which(extract)
    }

    valid.coords <- arrayInd(is.valid, dim(x))
    row.index <- anchors(x, type="row", id=TRUE)[valid.coords[,1]]
    col.index <- anchors(x, type="column", id=TRUE)[valid.coords[,2]]

    all.values <- as.matrix(x)[is.valid]
    dim(all.values) <- c(length(all.values), 1L)
    colnames(all.values) <- "1"
        
    final <- InteractionSet(all.values, GInteractions(row.index, col.index, regions(x), 
                mode=ifelse(collapse, "strict", "normal")), ...)
    if (collapse) {
        final <- unique(final)
    }
    return(final)
})

# Convert between GInteractions objects of different strictness.

setAs("GInteractions", "StrictGInteractions", function(from) {
    old_val <- S4Vectors:::disableValidity() # just in case we convert from reverse to strict.
    on.exit(S4Vectors:::disableValidity(old_val))
    S4Vectors:::disableValidity(TRUE)
    new("StrictGInteractions", swapAnchors(from, mode="order"))
})

setAs("GInteractions", "ReverseStrictGInteractions", function(from) {
    old_val <- S4Vectors:::disableValidity()
    on.exit(S4Vectors:::disableValidity(old_val))
    S4Vectors:::disableValidity(TRUE)
    new("ReverseStrictGInteractions", swapAnchors(from, mode="reverse"))
})

# Convert to a Hits object.

setAs("GInteractions", "SelfHits", function(from) {
    SelfHits(from=from@anchor1, to=from@anchor2, nnode=length(regions(from)), sort.by.query=FALSE)
})
