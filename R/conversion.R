# Inflate from InteractionSet to ContactMatrix.

setGeneric("inflate", function(x, ...) { standardGeneric("inflate") })

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

setMethod("inflate", "GInteractions", function(x, rows, columns, fill, ...) {
    row.chosen <- .make_to_indices(regions(x), rows, ...)
    col.chosen <- .make_to_indices(regions(x), columns, ...)
    fill <- rep(fill, length.out=nrow(x))
     
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
    out.mat <- matrix(NA, nR, nC)

    relevantA <- !is.na(ar1) & !is.na(ac2)
    out.mat[(ac2[relevantA] - 1L) * nR + ar1[relevantA]] <- fill[relevantA] 
    relevantB <- !is.na(ar2) & !is.na(ac1)
    out.mat[(ac1[relevantB] - 1L) * nR + ar2[relevantB]] <- fill[relevantB] 

    # Restoring the original order.
    original.rows <- cumsum(rnd)
    original.rows[ro] <- original.rows
    original.cols <- cumsum(cnd)
    original.cols[co] <- original.cols

    return(ContactMatrix(out.mat[original.rows,original.cols,drop=FALSE], 
                row.chosen[original.rows], col.chosen[original.cols], regions(x)))
})
 
setMethod("inflate", "InteractionSet", function(x, rows, columns, assay=1, sample=1, fill=NULL, ...) {
    if (length(fill)==0L) { fill <- assay(x, assay)[,sample] }
    inflate(interactions(x), rows, columns, fill=fill, ...)
})

setGeneric("deflate", function(x, ...) { standardGeneric("deflate") })

setMethod("deflate", "ContactMatrix", function(x, unique=TRUE, ...) {
    is.valid <- !is.na(as.matrix(x))
    valid.coords <- which(is.valid, arr.ind=TRUE)
    row.index <- anchors(x, type="row", id=TRUE)[valid.coords[,1]]
    col.index <- anchors(x, type="column", id=TRUE)[valid.coords[,2]]

    out <- .enforce_order(row.index, col.index)
    all.values <- as.matrix(x)[is.valid]
    dim(all.values) <- c(length(all.values), 1L)
    final <- InteractionSet(all.values, GInteractions(out$anchor1, out$anchor2, regions(x)), ...)

    if (unique) { final <- unique(final) }
    return(final)
})

