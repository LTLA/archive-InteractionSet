# Tests the construction and manipulation of InteractionSet objects.

set.seed(1000)
N <- 30
all.starts <- round(runif(N, 1, 100))
all.ends <- all.starts + round(runif(N, 5, 20))
all.regions <- GRanges(rep(c("chrA", "chrB"), c(N-10, 10)), IRanges(all.starts, all.ends))

Np <- 20
all.anchor1 <- sample(N, Np)
all.anchor2 <- sample(N, Np)
Nlibs <- 4
counts <- matrix(rpois(Np*Nlibs, lambda=10), ncol=Nlibs)
x <- InteractionSet(counts, all.anchor1, all.anchor2, all.regions)

expect_output(show(x), "class: InteractionSet 
dim: 20 4 
metadata(0):
assays(1): ''
rownames: NULL
metadata column names(0):
colnames: NULL
colData names(0):
regions: 30", 
fixed=TRUE)

# Testing all of our new slots:

expect_that(x, is_a("InteractionSet"))
expect_that(assay(x), is_equivalent_to(counts))
expect_true(all(anchors(x, id=TRUE, type="first") >= anchors(x, id=TRUE, type="second")))
expect_true(!is.unsorted(regions(x)))

o <- order(all.regions)
new.regions <- all.regions[o]
expect_that(regions(x), is_identical_to(new.regions))

new.pos <- integer(length(o))
new.pos[o] <- seq_along(new.pos)
new.anchor1 <- pmax(new.pos[all.anchor1], new.pos[all.anchor2])
new.anchor2 <- pmin(new.pos[all.anchor1], new.pos[all.anchor2])
expect_that(anchors(x, id=TRUE, type="first"), is_identical_to(new.anchor1))
expect_that(anchors(x, id=TRUE, type="second"), is_identical_to(new.anchor2))
expect_that(anchors(x, id=TRUE), is_identical_to(list(first=new.anchor1, second=new.anchor2)))

expect_that(anchors(x, type="first"), is_identical_to(new.regions[new.anchor1]))
expect_that(anchors(x, type="second"), is_identical_to(new.regions[new.anchor2]))
expect_that(anchors(x), is_identical_to(GRangesList(first=new.regions[new.anchor1], 
        second=new.regions[new.anchor2])))

# Testing alternative construction methods:

x2 <- InteractionSet(counts, all.regions[all.anchor1], all.regions[all.anchor2])
was.used <- sort(unique(all.regions[union(all.anchor1, all.anchor2)])) # Only includes the regions actually used.
expect_that(regions(x2), is_identical_to(was.used))
expect_that(anchors(x2), is_identical_to(anchors(x)))

x3 <- InteractionSet(counts, all.regions[all.anchor1], all.regions[all.anchor2], was.used)
expect_that(anchors(x3, id=TRUE), is_identical_to(anchors(x2, id=TRUE)))
expect_that(regions(x3), is_identical_to(regions(x2)))

# Testing with crappy inputs:

expect_that(InteractionSet(matrix(0, 4, 0), 1:4, 1:4, all.regions), is_a("InteractionSet")) # No columns.
expect_that(InteractionSet(matrix(0, 0, 4), integer(0), numeric(0), GRanges()), is_a("InteractionSet")) # No rows.
expect_that(InteractionSet(matrix(0, 0, 4), GRanges(), GRanges()), is_a("InteractionSet")) # No rows.

expect_error(InteractionSet(matrix(0, 4, 0), 1:4, 1, all.regions), "first and second anchor vectors have different lengths")
expect_error(InteractionSet(matrix(0, 4, 0), 0:3, 1:4, all.regions), "all anchor indices must be positive integers")
expect_error(InteractionSet(matrix(0, 4, 0), c(1,2,3,NA), 1:4, all.regions), "all anchor indices must be finite integers")
expect_error(InteractionSet(matrix(0, 3, 0), 1:4, 1:4, all.regions), "'assays' nrow differs from 'mcols' nrow")

# Testing setters.

set.seed(1001)
shuffled <- sample(100, N, replace=TRUE)
regions(x)$score <- shuffled
expect_that(regions(x)$score, is_identical_to(shuffled))
expect_false(identical(regions(x), new.regions))
regions(x) <- new.regions # Restoring.
expect_true(identical(regions(x), new.regions))

fresh.anchor1 <- sample(N, Np)
fresh.anchor2 <- sample(N, Np)
anchors(x) <- list(fresh.anchor1, fresh.anchor2)
expect_that(anchors(x, id=TRUE, type="first"), is_identical_to(pmax(fresh.anchor1, fresh.anchor2)))
expect_that(anchors(x, id=TRUE, type="second"), is_identical_to(pmin(fresh.anchor1, fresh.anchor2)))
anchors(x) <- list(new.anchor1, new.anchor2) # Restoring.

lib.sizes <- 1:4*1000L
x$totals <- lib.sizes
expect_that(x$totals, is_identical_to(lib.sizes))
expect_that(colData(x)$totals, is_identical_to(lib.sizes))

x.dump <- x
mod.ranges <- resize(regions(x), fix="center", width=50)
new.ranges <- c(regions(x), mod.ranges) 
expect_error(regions(x.dump) <- new.ranges, "assigned value must be of the same length")
newRegions(x.dump) <- new.ranges
expect_identical(anchors(x.dump), anchors(x))
expect_error(newRegions(x.dump) <- mod.ranges, "some existing ranges do not exist in replacement GRanges")

# Testing the subsetting.

rchosen <- 1:10
xsub <- x[rchosen,]
expect_output(show(xsub), "class: InteractionSet 
dim: 10 4 
metadata(0):
assays(1): ''
rownames: NULL
metadata column names(0):
colnames: NULL
colData names(1): totals
regions: 30", 
fixed=TRUE)

expect_that(assay(xsub), is_identical_to(assay(x)[rchosen,]))
expect_that(regions(xsub), is_identical_to(regions(x)))
expect_that(anchors(xsub, type="first"), is_identical_to(new.regions[new.anchor1][rchosen]))
expect_that(anchors(xsub, type="second"), is_identical_to(new.regions[new.anchor2][rchosen]))

cchosen <- c(2,4)
xsub <- x[,cchosen]
expect_output(show(xsub), "class: InteractionSet 
dim: 20 2 
metadata(0):
assays(1): ''
rownames: NULL
metadata column names(0):
colnames: NULL
colData names(1): totals
regions: 30", 
fixed=TRUE)

expect_that(assay(xsub), is_identical_to(assay(x)[,cchosen]))
expect_that(xsub$totals, is_identical_to(x$totals[cchosen]))
expect_that(regions(xsub), is_identical_to(regions(x)))
expect_that(anchors(xsub, type="first"), is_identical_to(new.regions[new.anchor1]))
expect_that(anchors(xsub, type="second"), is_identical_to(new.regions[new.anchor2]))

xsub <- subset(x, rchosen,cchosen)
expect_output(show(xsub), "class: InteractionSet 
dim: 10 2 
metadata(0):
assays(1): ''
rownames: NULL
metadata column names(0):
colnames: NULL
colData names(1): totals
regions: 30", 
fixed=TRUE)

expect_that(assay(xsub), is_identical_to(assay(x)[rchosen,cchosen]))
expect_that(regions(xsub), is_identical_to(regions(x)))
expect_that(anchors(xsub, type="first"), is_identical_to(new.regions[new.anchor1][rchosen]))
expect_that(anchors(xsub, type="second"), is_identical_to(new.regions[new.anchor2][rchosen]))

# Testing the combining.

xsub <- x[1:5,]
xsub2 <- x[6:20,]
expect_that(rbind(xsub, xsub2), equals(x))
expect_error(rbind(xsub, xsub2[,1:2]), "objects must have the same number of samples")
xsub <- x[,1]
xsub2 <- x[,2:4]
expect_that(cbind(xsub, xsub2), equals(x))
expect_error(cbind(xsub, xsub2[1:10,]), "anchors must be identical")

expect_that(nrow(rbind(x[0,], x[0,])), is_identical_to(0L)) # Behaviour with empties.
expect_that(ncol(rbind(x[0,], x[0,])), is_identical_to(ncol(x)))
expect_that(rbind(x, x[0,]), equals(x))
expect_that(nrow(cbind(x[,0], x[,0])), is_identical_to(nrow(x)))
expect_that(ncol(cbind(x[,0], x[,0])), is_identical_to(0L))
expect_that(cbind(x, x[,0]), equals(x))

set.seed(1002)
N <- 30
next.starts <- round(runif(N, 1, 100))
next.ends <- next.starts + round(runif(N, 5, 20))
next.regions <- GRanges(rep(c("chrA", "chrB"), c(N-10, 10)), IRanges(next.starts, next.ends))

Np <- 20
next.anchor1 <- sample(N, Np)
next.anchor2 <- sample(N, Np)
Nlibs <- 4
counts <- matrix(rpois(Np*Nlibs, lambda=10), ncol=Nlibs)
next.x <- InteractionSet(counts, next.anchor1, next.anchor2, next.regions)

expect_error(rbind(x, next.x), "regions must be identical in 'rbind'") 
c.x <- c(x, next.x)
expect_that(c(anchors(x, type="first"), anchors(next.x, type="first")), is_identical_to(anchors(c.x, type="first")))
expect_that(c(anchors(x, type="second"), anchors(next.x, type="second")), is_identical_to(anchors(c.x, type="second")))
expect_that(unique(sort(c(regions(x), regions(next.x)))), is_identical_to(regions(c.x)))

expect_that(nrow(c(x[0,], next.x[0,])), is_identical_to(0L)) # Behaviour with empties.
expect_that(ncol(c(x[0,], next.x[0,])), is_identical_to(ncol(x)))
expect_that(nrow(c(x, next.x[0,])), is_identical_to(nrow(x))) # Not fully equal, as regions have changed.

# Testing the sorting.

o.x <- order(anchors(x, type="first"), anchors(x, type="second"))
expect_that(o.x, is_identical_to(order(x)))
expect_that(sort(x), equals(x[o.x,]))

o.x2 <- order(anchors(x, type="first"), anchors(x, type="second"), anchors(next.x, type="first"), anchors(next.x, type="second"))
expect_that(o.x2, is_identical_to(order(x, next.x)))

is.dup <- duplicated(paste0(anchors(x, type="first"), ".", anchors(x, type="second")))
expect_that(is.dup, is_identical_to(duplicated(x)))
temp.x <- rbind(x, x)    
is.dup <- duplicated(paste0(anchors(temp.x, type="first"), ".", anchors(temp.x, type="second")))
expect_that(is.dup, is_identical_to(duplicated(temp.x)))
expect_true(all(tail(is.dup, length(x)))) # if ordering is stable; only the first occurrence should be true.
expect_equal(x, unique(temp.x))

is.dup <- duplicated(paste0(anchors(temp.x, type="first"), ".", anchors(temp.x, type="second")), fromLast=TRUE)
expect_that(is.dup, is_identical_to(duplicated(temp.x, fromLast=TRUE)))
expect_true(all(head(is.dup, length(x)))) # if ordering is stable; only the first occurrence should be true.
expect_equal(x, unique(temp.x, fromLast=TRUE))

# Testing the splitting.

flen <- c(5L, 10L, 5L)
f <- rep(1:3, flen)
out <- split(x, f)
expect_that(sapply(out, nrow), is_equivalent_to(flen))

