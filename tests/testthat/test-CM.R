# Tests the construction and manipulation of InteractionSet objects.

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
counts <- matrix(rpois(Nr*Nc, lambda=10), Nr, Nc)
x <- ContactMatrix(counts, all.anchor1, all.anchor2, all.regions)

expect_output(show(x), sprintf("class: ContactMatrix 
dim: %i %i 
regions: %i", Nr, Nc, N), 
fixed=TRUE)

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

expect_that(ContactMatrix(matrix(0, 4, 0), 1:4, integer(0), all.regions), is_a("ContactMatrix")) # No columns.
four.peat <- GRanges("chrA", IRanges(1:4, 1:4))
expect_that(ContactMatrix(matrix(0, 0, 4), integer(0), 1:4, four.peat), is_a("ContactMatrix")) # No rows.
expect_that(ContactMatrix(matrix(0, 0, 4), GRanges(), four.peat), is_a("ContactMatrix")) # Nothing at all

expect_error(ContactMatrix(matrix(0, 3, 1), 1:4, 1, all.regions), "'matrix' nrow must be equal to length of 'anchor1'")
expect_error(ContactMatrix(matrix(0, 4, 0), 1:4, 1, all.regions), "'matrix' ncol must be equal to length of 'anchor2'")
expect_error(ContactMatrix(matrix(0, 4, 0), 0:3, 1:4, all.regions), "all anchor indices must be positive integers")
expect_error(ContactMatrix(matrix(0, 4, 0), c(1,2,3,NA), 1:4, all.regions), "all anchor indices must be finite integers")
expect_error(ContactMatrix(matrix(0, 4, 0), c(1,2,3,-1), 1:4, all.regions), "all anchor indices must be positive integers")
expect_error(ContactMatrix(matrix(0, 4, 0), c(1,2,3,length(all.regions)+1L), 1:4, all.regions), "all anchor indices must refer to entries in 'regions'")

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
anchors(x) <- list(fresh.anchor1, fresh.anchor2)
expect_that(anchors(x, id=TRUE, type="row"), is_identical_to(fresh.anchor1))
expect_that(anchors(x, id=TRUE, type="column"), is_identical_to(fresh.anchor2))
expect_error(anchors(x) <- list(fresh.anchor1, fresh.anchor2, fresh.anchor1), "must be a list of 2 numeric vectors")
expect_error(anchors(x) <- list(fresh.anchor2, fresh.anchor2), "nrow must be equal to length of 'anchor1'")
expect_error(anchors(x) <- list(fresh.anchor1, fresh.anchor1), "ncol must be equal to length of 'anchor2'")
anchors(x) <- list(new.anchor1, new.anchor2) # Restoring.

x.dump <- x
mod.ranges <- resize(regions(x), fix="center", width=50)
new.ranges <- c(regions(x), mod.ranges) 
expect_error(regions(x.dump) <- new.ranges, "assigned value must be of the same length")
newRegions(x.dump) <- new.ranges
expect_identical(anchors(x.dump), anchors(x))
expect_error(newRegions(x.dump) <- mod.ranges, "some existing ranges do not exist in replacement GRanges")

# Testing the subsetting.

rchosen <- 1:5
xsub <- x[rchosen,]
expect_output(show(xsub), "class: ContactMatrix 
dim: 5 20 
regions: 30", 
fixed=TRUE)

expect_that(as.matrix(xsub), is_identical_to(as.matrix(x)[rchosen,]))
expect_that(regions(xsub), is_identical_to(regions(x)))
expect_that(anchors(xsub, type="row"), is_identical_to(new.regions[new.anchor1][rchosen]))
expect_that(anchors(xsub, type="column"), is_identical_to(new.regions[new.anchor2]))

cchosen <- 10:20
xsub <- x[,cchosen]
expect_output(show(xsub), "class: ContactMatrix 
dim: 10 11 
regions: 30", 
fixed=TRUE)

expect_that(as.matrix(xsub), is_identical_to(as.matrix(x)[,cchosen]))
expect_that(regions(xsub), is_identical_to(regions(x)))
expect_that(anchors(xsub, type="row"), is_identical_to(new.regions[new.anchor1]))
expect_that(anchors(xsub, type="column"), is_identical_to(new.regions[new.anchor2][cchosen]))

xsub <- subset(x,rchosen,cchosen)
expect_output(show(xsub), "class: ContactMatrix 
dim: 5 11 
regions: 30", 
fixed=TRUE)

expect_that(as.matrix(xsub), is_identical_to(as.matrix(x)[rchosen,cchosen]))
expect_that(regions(xsub), is_identical_to(regions(x)))
expect_that(anchors(xsub, type="row"), is_identical_to(new.regions[new.anchor1][rchosen]))
expect_that(anchors(xsub, type="column"), is_identical_to(new.regions[new.anchor2][cchosen]))

expect_that(nrow(x[0,]), is_identical_to(0L))
expect_that(ncol(x[,0]), is_identical_to(0L))

# Testing the combining.

xsub <- x[1:5,]
xsub2 <- x[6:10,]
expect_that(rbind(xsub, xsub2), equals(x))
expect_error(rbind(xsub, xsub2[,1:2]), "column anchor indices must be identical")
xsub <- x[,1:5]
xsub2 <- x[,6:20]
expect_that(cbind(xsub, xsub2), equals(x))
expect_error(cbind(xsub, xsub2[1:5,]), "row anchor indices must be identical")

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
anchors(temp.x2) <- list(seq_len(nrow(temp.x2)), anchors(temp.x2, type="column", id=TRUE))
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
