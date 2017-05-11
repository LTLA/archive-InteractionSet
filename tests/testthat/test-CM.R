# Tests the construction and manipulation of ContactMatrix objects.

### Start of loop.
for (type in c("normal", "sparse", "dense")) { 
if (type=="sparse") { 
    mattype <- "dgCMatrix"
    makeMatrix <- function(x, r, c) {
        Matrix::sparseMatrix(rep(seq_len(r), c), rep(seq_len(c), each=r), x=x, dims=c(r, c))
    }
} else if (type=="normal") {
    makeMatrix <- base::matrix
    mattype <- "matrix"
} else {
    makeMatrix <- Matrix::Matrix
    mattype <- "dgeMatrix"
}
###

set.seed(4000)
N <- 30
all.starts <- round(runif(N, 1, 100))
all.ends <- all.starts + round(runif(N, 5, 20))
all.regions <- GRanges(rep(c("chrA", "chrB"), c(N-10, 10)),
                           IRanges(all.starts, all.ends))

Nr <- 10
Nc <- 20
all.anchor1 <- sample(N, Nr)
all.anchor2 <- sample(N, Nc)
counts <- makeMatrix(rpois(Nr*Nc, lambda=10), Nr, Nc)
x <- ContactMatrix(counts, all.anchor1, all.anchor2, all.regions)

expect_output(show(x), sprintf("class: ContactMatrix 
dim: %i %i 
type: %s 
rownames: NULL
colnames: NULL
metadata(0):
regions: %i", Nr, Nc, mattype, N), 
fixed=TRUE)

temp.x <- ContactMatrix(counts, all.anchor1, all.anchor2, all.regions, metadata=list("whee"=1))
expect_output(show(temp.x), sprintf("class: ContactMatrix 
dim: %i %i 
type: %s 
rownames: NULL
colnames: NULL
metadata(1): whee
regions: %i", Nr, Nc, mattype, N),
fixed=TRUE)
metadata(temp.x)$whee <- NULL
expect_identical(temp.x, x)

# Testing all of our new slots:

expect_that(x, is_a("ContactMatrix"))
expect_that(as.matrix(x), equals(counts))
expect_true(!is.unsorted(regions(x)))

o <- order(all.regions)
new.regions <- all.regions[o]
expect_that(regions(x), is_identical_to(new.regions))

new.pos <- integer(length(o))
new.pos[o] <- seq_along(new.pos)
new.anchor1 <- new.pos[all.anchor1]
new.anchor2 <- new.pos[all.anchor2]
expect_that(anchors(x, id=TRUE, type="row"), is_identical_to(new.anchor1))
expect_that(anchors(x, id=TRUE, type="column"), is_identical_to(new.anchor2))
expect_that(anchors(x, id=TRUE), is_identical_to(list(row=new.anchor1, column=new.anchor2)))

expect_that(anchors(x, type="row"), is_identical_to(new.regions[new.anchor1]))
expect_that(anchors(x, type="column"), is_identical_to(new.regions[new.anchor2]))
expect_that(anchors(x), is_identical_to(GRangesList(row=new.regions[new.anchor1], 
        column=new.regions[new.anchor2])))

# Testing alternative construction methods:

x2 <- ContactMatrix(counts, all.regions[all.anchor1], all.regions[all.anchor2])
was.used <- sort(unique(all.regions[union(all.anchor1, all.anchor2)])) # Only includes the regions actually used.
expect_that(regions(x2), is_identical_to(was.used))
expect_that(anchors(x2), is_identical_to(anchors(x)))

x3 <- ContactMatrix(counts, all.regions[all.anchor1], all.regions[all.anchor2], was.used)
expect_that(anchors(x3, id=TRUE), is_identical_to(anchors(x2, id=TRUE)))
expect_that(regions(x3), is_identical_to(regions(x2)))

# Testing with crappy inputs:

expect_that(ContactMatrix(makeMatrix(0, 4, 0), 1:4, integer(0), all.regions), is_a("ContactMatrix")) # No columns.
four.peat <- GRanges("chrA", IRanges(1:4, 1:4))
expect_that(ContactMatrix(makeMatrix(0, 0, 4), integer(0), 1:4, four.peat), is_a("ContactMatrix")) # No rows.
expect_that(ContactMatrix(makeMatrix(0, 0, 4), GRanges(), four.peat), is_a("ContactMatrix")) # Nothing at all

expect_error(ContactMatrix(makeMatrix(0, 3, 1), 1:4, 1, all.regions), "'matrix' nrow must be equal to length of 'anchor1'")
expect_error(ContactMatrix(makeMatrix(0, 4, 0), 1:4, 1, all.regions), "'matrix' ncol must be equal to length of 'anchor2'")
expect_error(ContactMatrix(makeMatrix(0, 4, 0), 0:3, 1:4, all.regions), "all anchor indices must be positive integers")
expect_error(ContactMatrix(makeMatrix(0, 4, 0), c(1,2,3,NA), 1:4, all.regions), "all anchor indices must be finite integers")
expect_error(ContactMatrix(makeMatrix(0, 4, 0), c(1,2,3,-1), 1:4, all.regions), "all anchor indices must be positive integers")
expect_error(ContactMatrix(makeMatrix(0, 4, 0), c(1,2,3,length(all.regions)+1L), 1:4, all.regions), "all anchor indices must refer to entries in 'regions'")

expect_identical(dim(ContactMatrix()), integer(2))
if (type=="normal") { expect_identical(ContactMatrix(), ContactMatrix(makeMatrix(0L,0,0), integer(0), integer(0), GRanges())) }

# Testing setters.

set.seed(4001)
shuffled <- sample(100, N, replace=TRUE)
regions(x)$score <- shuffled
expect_that(regions(x)$score, is_identical_to(shuffled))
expect_false(identical(regions(x), new.regions))
regions(x) <- new.regions # Restoring.
expect_true(identical(regions(x), new.regions))

fresh.anchor1 <- sample(N, Nr)
fresh.anchor2 <- sample(N, Nc)
anchorIds(x) <- list(fresh.anchor1, fresh.anchor2)
expect_that(anchors(x, id=TRUE, type="row"), is_identical_to(fresh.anchor1))
expect_that(anchors(x, id=TRUE, type="column"), is_identical_to(fresh.anchor2))
expect_error(anchorIds(x) <- list(fresh.anchor1, fresh.anchor2, fresh.anchor1), "must be a list of 2 numeric vectors")
expect_error(anchorIds(x) <- list(fresh.anchor2, fresh.anchor2), "nrow must be equal to length of 'anchor1'")
expect_error(anchorIds(x) <- list(fresh.anchor1, fresh.anchor1), "ncol must be equal to length of 'anchor2'")

mod.x <- x
anchorIds(x, type="row") <- new.anchor1 # Restoring; checking that these calls also work.
expect_identical(anchors(x, id=TRUE, type="row"), new.anchor1)
anchorIds(x, type="column") <- new.anchor2
expect_identical(anchors(x, id=TRUE, type="column"), new.anchor2)
anchorIds(mod.x, type="both") <- list(new.anchor1, new.anchor2) # Restoring.
expect_identical(x, mod.x)

x.dump <- x
mod.ranges <- resize(regions(x), fix="center", width=50)
new.ranges <- c(regions(x), mod.ranges) 
expect_error(regions(x.dump) <- new.ranges, "assigned value must be of the same length")
replaceRegions(x.dump) <- new.ranges
expect_identical(anchors(x.dump), anchors(x))
expect_identical(sort(new.ranges), regions(x.dump))
expect_error(replaceRegions(x.dump) <- mod.ranges, "some existing ranges do not exist in replacement GRanges")

x.dump2 <- x
appendRegions(x.dump2) <- mod.ranges
expect_identical(regions(x.dump), regions(x.dump2))
expect_identical(anchors(x.dump), anchors(x.dump2))

x.dump <- reduceRegions(x)
expect_identical(anchors(x), anchors(x.dump))
expect_identical(regions(x)[sort(unique(unlist(anchors(x, id=TRUE))))], regions(x.dump))

x.dump <- x
new.mat <- makeMatrix(sample(N, Nr*Nc, replace=TRUE), Nr, Nc)
as.matrix(x.dump) <- new.mat
expect_identical(new.mat, as.matrix(x.dump))
as.matrix(x.dump) <- 1:Nr
expect_equal(as.matrix(x.dump), makeMatrix(1:Nr, Nr, Nc))

if (type=="sparse") { expect_error(as.matrix(x.dump) <- makeMatrix(1:Nr, Nr, 1), "replacement Matrix must have same dimensions as 'x'") }

# Testing the subsetting.

rchosen <- 1:5
xsub <- x[rchosen,]
expect_output(show(xsub), sprintf("class: ContactMatrix 
dim: 5 20 
type: %s 
rownames: NULL
colnames: NULL
metadata(0):
regions: 30", mattype),
fixed=TRUE)

expect_that(as.matrix(xsub), is_identical_to(as.matrix(x)[rchosen,]))
expect_that(regions(xsub), is_identical_to(regions(x)))
expect_that(anchors(xsub, type="row"), is_identical_to(new.regions[new.anchor1][rchosen]))
expect_that(anchors(xsub, type="column"), is_identical_to(new.regions[new.anchor2]))

cchosen <- 10:20
xsub <- x[,cchosen]
expect_output(show(xsub), sprintf("class: ContactMatrix 
dim: 10 11 
type: %s 
rownames: NULL
colnames: NULL
metadata(0):
regions: 30", mattype),
fixed=TRUE)

expect_that(as.matrix(xsub), is_identical_to(as.matrix(x)[,cchosen]))
expect_that(regions(xsub), is_identical_to(regions(x)))
expect_that(anchors(xsub, type="row"), is_identical_to(new.regions[new.anchor1]))
expect_that(anchors(xsub, type="column"), is_identical_to(new.regions[new.anchor2][cchosen]))

xsub <- subset(x,rchosen,cchosen)
expect_output(show(xsub), sprintf("class: ContactMatrix 
dim: 5 11 
type: %s 
rownames: NULL
colnames: NULL
metadata(0):
regions: 30", mattype),
fixed=TRUE)

expect_that(as.matrix(xsub), is_identical_to(as.matrix(x)[rchosen,cchosen]))
expect_that(regions(xsub), is_identical_to(regions(x)))
expect_that(anchors(xsub, type="row"), is_identical_to(new.regions[new.anchor1][rchosen]))
expect_that(anchors(xsub, type="column"), is_identical_to(new.regions[new.anchor2][cchosen]))

expect_that(nrow(x[0,]), is_identical_to(0L))
expect_that(ncol(x[,0]), is_identical_to(0L))

# Testing subset assignments

temp.x <- x
temp.x[rchosen+5,] <- x[rchosen,]
new.index <- seq_len(nrow(x))
new.index[rchosen+5] <- rchosen
expect_equal(as.matrix(temp.x), as.matrix(x)[new.index,])
expect_identical(anchors(temp.x, type="row"), anchors(x, type="row")[new.index,])
expect_identical(anchors(temp.x, type="column"), anchors(x, type="column"))

temp.x <- x
temp.x[,cchosen-9,] <- x[,cchosen]
new.index <- seq_len(ncol(x))
new.index[cchosen-9] <- cchosen
expect_equal(as.matrix(temp.x), as.matrix(x)[,new.index])
expect_identical(anchors(temp.x, type="row"), anchors(x, type="row"))
expect_identical(anchors(temp.x, type="column"), anchors(x, type="column")[new.index,])

temp.x <- x
temp.x[0,] <- x[0,]
expect_identical(temp.x, x)
temp.x[,0] <- x[,0]
expect_identical(temp.x, x)

temp.x <- x
anchorIds(temp.x[1:5,], type="row") <- 1:5
expect_identical(anchors(temp.x, type="row", id=TRUE)[1:5], 1:5)
anchorIds(temp.x[,1:5], type="column") <- 1:5
expect_identical(anchors(temp.x, type="column", id=TRUE)[1:5], 1:5)
expect_error(anchorIds(temp.x[1:5,1:5], type="row") <- 1:5+1L, "cannot modify row indices for a subset of columns")
expect_error(anchorIds(temp.x[1:5,1:5], type="column") <- 1:5+1L, "cannot modify column indices for a subset of rows")

# Testing the combining.

xsub <- x[1:5,]
xsub2 <- x[6:10,]
expect_that(rbind(xsub, xsub2), equals(x))
expect_error(rbind(xsub, xsub2[,1:2]), "column anchor indices must be identical")
xsub <- x[,1:5]
xsub2 <- x[,6:20]
expect_that(cbind(xsub, xsub2), equals(x))
expect_error(cbind(xsub, xsub2[1:5,]), "row anchor indices must be identical")

xsub.mod <- xsub <- x[1:5,]
chosen.cols <- setdiff(anchors(xsub.mod, type="column", id=TRUE), anchors(xsub.mod, type="row", id=TRUE))
regions(xsub.mod)[chosen.cols] <- resize(regions(xsub.mod)[chosen.cols], width=5)
byrow <- cbind(xsub, xsub.mod)
expect_identical(anchors(byrow, type="row"), anchors(xsub, type="row"))
expect_identical(anchors(byrow, type="column"), c(anchors(xsub, type="column"), anchors(xsub.mod, type="column")))
expect_equal(as.matrix(byrow), cbind(as.matrix(xsub), as.matrix(xsub.mod)))

xsub.mod <- xsub <- x[,1:5]
chosen.cols <- setdiff(anchors(xsub.mod, type="row", id=TRUE), anchors(xsub.mod, type="column", id=TRUE))
regions(xsub.mod)[chosen.cols] <- resize(regions(xsub.mod)[chosen.cols], width=5)
bycol <- rbind(xsub, xsub.mod)
expect_identical(anchors(bycol, type="row"), c(anchors(xsub, type="row"), anchors(xsub.mod, type="row")))
expect_identical(anchors(bycol, type="column"), anchors(xsub, type="column"))
expect_equal(as.matrix(bycol), rbind(as.matrix(xsub), as.matrix(xsub.mod)))

expect_that(nrow(rbind(x[0,], x[0,])), is_identical_to(0L)) # Behaviour with empties.
expect_that(ncol(rbind(x[0,], x[0,])), is_identical_to(ncol(x)))
expect_that(rbind(x, x[0,]), equals(x))
expect_that(nrow(cbind(x[,0], x[,0])), is_identical_to(nrow(x)))
expect_that(ncol(cbind(x[,0], x[,0])), is_identical_to(0L))
expect_that(cbind(x, x[,0]), equals(x))

# Testing the sorting.

o.x <- list(row=order(anchors(x, type="row")), column=order(anchors(x, type="column")))
expect_that(o.x, is_identical_to(order(x)))
expect_that(sort(x), equals(x[o.x$row,o.x$column]))

temp.x <- rbind(x, x)    
temp.x2 <- temp.x
anchorIds(temp.x2) <- list(seq_len(nrow(temp.x2)), anchors(temp.x2, type="column", id=TRUE))
o.x2 <- list(row=order(anchors(temp.x, type="row"), anchors(temp.x2, type="row")),
             column=order(anchors(temp.x, type="column"), anchors(temp.x2, type="column")))
expect_that(o.x2, is_identical_to(order(temp.x, temp.x2)))

is.dup <- list(row=duplicated(anchors(x, type="row")), column=duplicated(anchors(x, type="column")))
expect_that(is.dup, is_identical_to(duplicated(x)))

temp.x <- rbind(x, x)    
is.dup <- list(row=duplicated(anchors(temp.x, type="row")), column=duplicated(anchors(temp.x, type="column")))
expect_that(is.dup, is_identical_to(duplicated(temp.x)))
expect_true(all(tail(is.dup$row, nrow(x)))) # if ordering is stable; only the first occurrence should be true.
expect_equal(x, unique(temp.x))

temp.x <- cbind(x, x)    
is.dup <- list(row=duplicated(anchors(temp.x, type="row")), column=duplicated(anchors(temp.x, type="column")))
expect_that(is.dup, is_identical_to(duplicated(temp.x)))
expect_true(all(tail(is.dup$column, ncol(x)))) 
expect_equal(x, unique(temp.x))

temp.x <- rbind(temp.x, temp.x)
is.dup <- list(row=duplicated(anchors(temp.x, type="row"), fromLast=TRUE), column=duplicated(anchors(temp.x, type="column"), fromLast=TRUE))
expect_that(is.dup, is_identical_to(duplicated(temp.x, fromLast=TRUE)))
expect_true(all(head(is.dup$column, ncol(x)))) # if ordering is stable; only the first occurrence should be true.
expect_true(all(head(is.dup$row, nrow(x)))) # if ordering is stable; only the first occurrence should be true.
expect_equal(x, unique(temp.x, fromLast=TRUE))

dedupped <- duplicated(unique(temp.x))
expect_false(any(dedupped$row))
expect_false(any(dedupped$column))

# Testing the tranposition

tx <- t(x)
expect_identical(anchors(x, type="row"), anchors(tx, type="column"))
expect_identical(anchors(x, type="column"), anchors(tx, type="row"))
expect_identical(t(x[0,]), tx[,0])
expect_identical(t(x[,0]), tx[0,])

# Testing name setting and extraction.

temp.x <- x
rowref <- paste0("X", seq_len(nrow(temp.x)))
rownames(temp.x) <- rowref
colref <- paste0("Y", seq_len(ncol(temp.x)))
colnames(temp.x) <- colref
expect_output(show(temp.x), sprintf("class: ContactMatrix 
dim: %i %i 
type: %s 
rownames(10): X1 X2 ... X9 X10
colnames(20): Y1 Y2 ... Y19 Y20
metadata(0):
regions: %i", nrow(temp.x), ncol(temp.x), mattype, length(regions(temp.x))), 
fixed=TRUE)

expect_identical(rownames(temp.x), rowref)
expect_identical(rownames(temp.x[1:5,]), rowref[1:5])
expect_identical(colnames(temp.x[,1:3]), colref[1:3])

rcombined <- cbind(x, temp.x)
expect_identical(rownames(rcombined), rowref)
expect_identical(colnames(rcombined), c(character(ncol(x)), colref))
ccombined <- rbind(x, temp.x)
expect_identical(rownames(ccombined), c(character(nrow(x)), rowref))
expect_identical(colnames(ccombined), colref)

for (id in c(TRUE, FALSE)) {
    expect_identical(names(anchors(temp.x, id=id)[[1]]), rowref)
    expect_identical(names(anchors(temp.x, id=id)[[2]]), colref)
    expect_identical(names(anchors(temp.x, id=id, type="row")), rowref)
    expect_identical(names(anchors(temp.x, id=id, type="column")), colref)
}

### End of loop.
}
###
