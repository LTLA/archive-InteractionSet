# Inflate from InteractionSet to ContactMatrix.

setGeneric("inflate", function(x, ...) { standardGeneric("inflate") })

.make_to_indices <- function(regs, i, ...) {
    nregs <- length(regs)
    if (is.numeric(i)) { 
        i <- as.integer(i)
        if (any(!is.finite(i)) || any(i<=0L) || any(i >= nregs)) { 
            stop("indices must be positive integers") 
        }
        return(i)
    } else if (is.character(i)) { 
        return(which(seqnames(regs) %in% i))
    } else if (is(i, "GRanges")) {
        return(which(overlapsAny(regs, i, ...)))
    }
}


setMethod("inflate", "InteractionSet", function(x, rows, columns, assay=1, sample=1, fill=NULL, ...) {
    row.chosen <- .make_to_indices(regions(x), rows, ...)
    col.chosen <- .make_to_indices(regions(x), columns, ...)
    if (length(fill)==0L) { fill <- assay(x, assay)[,sample] }
    else { fill <- rep(fill, length.out=nrow(x)) }
     
    # Removing duplicated rows and resorting (we'll put them back in later)
    ro <- order(row.chosen)
    co <- order(col.chosen)
    row.chosen <- row.chosen[ro]
    col.chosen <- col.chosen[co]
    rnd <- !duplicated(row.chosen)
    cnd <- !duplicated(col.chosen)
    row.chosen2 <- row.chosen[rnd]   
    col.chosen2 <- col.chosen[cnd]

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

setGeneric("deflate", function(x, ...) { standardGeneric("deflate") })

setMethod("deflate", "ContactMatrix", function(x, unique=TRUE, ...) {
    if (unique) { x <- unique(x) }
    all.values <- as.vector(as.matrix(x))
    row.index <- rep(seq_len(nrow(x)), ncol(x))
    col.index <- rep(seq_len(ncol(x)), each=nrow(x))

    is.valid <- !is.na(all.values)
    row.index <- row.index[is.valid]
    col.index <- col.index[is.valid]
    all.values <- all.values[is.valid]

    out <- .enforce_order(row.index, col.index)
    dim(all.values) <- c(length(all.values), 1L)
    return(InteractionSet(all.values, out$anchor1, out$anchor2, regions(x), ...))
})

