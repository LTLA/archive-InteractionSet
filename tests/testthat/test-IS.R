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
colnames(counts) <- seq_len(Nlibs)
x <- InteractionSet(counts, GInteractions(all.anchor1, all.anchor2, all.regions))

expect_output(show(x), "class: InteractionSet 
dim: 20 4 
metadata(0):
assays(1): ''
rownames: NULL
rowData names(0):
colnames(4): 1 2 3 4
colData names(0):
type: GInteractions
regions: 30", 
fixed=TRUE)

# Testing all of our new slots:

expect_is(x, "InteractionSet")
expect_is(interactions(x), "GInteractions")
expect_equivalent(assay(x), counts)
ref <- GInteractions(all.anchor1, all.anchor2, all.regions)
expect_identical(interactions(x), ref)

scores <- Nlibs:1 + 50L
x2 <- InteractionSet(counts, ref, colData=DataFrame(score=scores))
expect_identical(x2$score, scores)
expect_identical(colData(x2)$score, scores)
expect_identical(mcols(x2)$score, NULL) # Dollar doesn't assign here.
                     
x3 <- InteractionSet(counts, ref, metadata=list(whee=5L))
expect_identical(metadata(x3)$whee, 5L)

# Testing with crappy inputs:

expect_identical(dim(InteractionSet(matrix(0, 4, 0), GInteractions(1:4, 1:4, all.regions))), c(4L, 0L)) # No columns
expect_identical(dim(InteractionSet(matrix(0, 0, 4, dimnames=list(NULL, seq_len(4))), GInteractions(integer(0), numeric(0), GRanges()))), c(0L, 4L))
expect_error(InteractionSet(matrix(0, 3, 0), GInteractions(1:4, 1:4, all.regions)), "'interactions' length is not equal to the number of rows")

# Testing getters, setters.

set.seed(1001)
shuffled <- sample(100, N, replace=TRUE)
ref <- interactions(x)
expect_identical(regions(x), regions(ref))
expect_identical(anchors(x), anchors(ref))
expect_identical(anchors(x, id=TRUE), anchors(ref, id=TRUE))
expect_identical(anchors(x, type="first"), anchors(ref, type="first"))
expect_identical(anchors(x, type="second"), anchors(ref, type="second"))
expect_identical(first(x), first(ref))
expect_identical(second(x), second(ref))

regions(x)$score <- shuffled
regions(ref)$score <- shuffled
expect_identical(regions(x)$score, regions(ref)$score)
regions(ref)$score <- regions(x)$score <- NULL # Restoring.

orig.x <- x
fresh.anchor1 <- sample(N, Np)
fresh.anchor2 <- sample(N, Np)
anchorIds(x) <- list(fresh.anchor1, fresh.anchor2)
anchorIds(ref) <- list(fresh.anchor1, fresh.anchor2)
expect_identical(anchors(x), anchors(ref))
expect_identical(anchors(x, id=TRUE), anchors(ref, id=TRUE))
expect_identical(anchors(x, type="first"), anchors(ref, type="first"))
expect_identical(anchors(x, type="second"), anchors(ref, type="second"))
expect_identical(first(x), first(ref))
expect_identical(second(x), second(ref))

original.i <- anchors(orig.x, id=TRUE) # Reverting back to original indices, to check individual assignments work.
anchorIds(x, type="first") <- original.i$first
anchorIds(ref, type="first") <- original.i$first
expect_identical(anchors(x, id=TRUE, type="first"), anchors(ref, id=TRUE, type="first"))
anchorIds(x, type="second") <- original.i$second
anchorIds(ref, type="second") <- original.i$second
expect_identical(anchors(x, id=TRUE, type="second"), anchors(ref, id=TRUE, type="second"))
anchorIds(x, type="both") <- original.i
anchorIds(ref, type="both") <- original.i
expect_identical(anchors(x), anchors(ref))

lib.sizes <- 1:4*1000L
x$totals <- lib.sizes
expect_identical(x$totals, lib.sizes)
expect_identical(colData(x)$totals, lib.sizes)

x.dump <- x
ref.dump <- interactions(x)
mod.ranges <- resize(regions(x), fix="center", width=50)
new.ranges <- c(regions(x), mod.ranges) 
expect_error(regions(x.dump) <- new.ranges, "assigned value must be of the same length")
replaceRegions(x.dump) <- new.ranges
replaceRegions(ref.dump) <- new.ranges
expect_identical(anchors(x.dump), anchors(ref.dump))
expect_identical(regions(x.dump), regions(ref.dump))
expect_error(replaceRegions(x.dump) <- mod.ranges, "some existing ranges do not exist in replacement GRanges")

x.dump <- x
ref.dump <- interactions(x)
appendRegions(x.dump) <- mod.ranges
appendRegions(ref.dump) <- mod.ranges
expect_identical(regions(x.dump), regions(ref.dump))

expect_identical(interactions(reduceRegions(x)), reduceRegions(interactions(x)))

x.dump <- x
interactions(x.dump) <- rev(interactions(x))
expect_identical(interactions(x.dump), rev(interactions(x)))
expect_identical(interactions(rev(x.dump)), interactions(x))

new.scores <- 1:Np*10 + 50
mcols(x.dump)$fire <- new.scores
expect_identical(mcols(x.dump)$fire, new.scores)
expect_identical(interactions(x.dump)$fire, new.scores)

new.si <- Seqinfo(seqnames=c("chrA", "chrB"), seqlengths=c(1000, 2000))
new.x <- x
seqinfo(new.x) <- new.si
expect_identical(seqinfo(new.x), seqinfo(interactions(new.x)))

# Testing the subsetting.

rchosen <- 1:10
xsub <- x[rchosen,]
expect_output(show(xsub), "class: InteractionSet 
dim: 10 4 
metadata(0):
assays(1): ''
rownames: NULL
rowData names(0):
colnames(4): 1 2 3 4
colData names(1): totals
type: GInteractions
regions: 30", 
fixed=TRUE)
expect_equal(xsub, x[rchosen])
expect_identical(assay(xsub), assay(x)[rchosen,])
expect_identical(interactions(xsub), interactions(x)[rchosen])

cchosen <- c(2,4)
xsub <- x[,cchosen]
expect_output(show(xsub), "class: InteractionSet 
dim: 20 2 
metadata(0):
assays(1): ''
rownames: NULL
rowData names(0):
colnames(2): 2 4
colData names(1): totals
type: GInteractions
regions: 30", 
fixed=TRUE)
expect_identical(assay(xsub), assay(x)[,cchosen])
expect_identical(xsub$totals, x$totals[cchosen])
expect_identical(interactions(xsub), interactions(x))

xsub <- x[rchosen,cchosen]
expect_output(show(xsub), "class: InteractionSet 
dim: 10 2 
metadata(0):
assays(1): ''
rownames: NULL
rowData names(0):
colnames(2): 2 4
colData names(1): totals
type: GInteractions
regions: 30", 
fixed=TRUE)
expect_equal(xsub, subset(x, rchosen, cchosen))
lrchosen <- logical(nrow(x)); lrchosen[rchosen] <- TRUE
lcchosen <- logical(ncol(x)); lcchosen[cchosen] <- TRUE
expect_equal(xsub, x[lrchosen,lcchosen])

expect_that(assay(xsub), is_identical_to(assay(x)[rchosen,cchosen]))
expect_identical(xsub$totals, x$totals[cchosen])
expect_identical(interactions(xsub), interactions(x)[rchosen])

expect_identical(nrow(x[0,]), 0L)
expect_identical(ncol(x[0,]), as.integer(Nlibs))
expect_identical(ncol(x[,0]), 0L)
expect_identical(nrow(x[,0]), as.integer(Np))

# Testing subset assignment.

temp.x <- x
temp.x[1:5+10,] <- x[1:5,]
new.index <- seq_len(nrow(x))
new.index[1:5+10] <- 1:5
expect_equal(assay(temp.x), assay(x)[new.index,])
expect_identical(interactions(temp.x), interactions(x)[new.index,])

temp.x <- x  
c.from <- 2:3
c.to <- 1:2
tmp <- x[,c.from]
colnames(tmp) <- colnames(x)[c.to] # Avoid errors from duplicated column names.
temp.x[,c.to] <- tmp
new.index <- seq_len(ncol(x))
new.index[c.to] <- c.from
ref.x <- x[,new.index]
colnames(ref.x) <- colnames(x)
expect_equal(assay(temp.x), assay(ref.x))
expect_identical(interactions(temp.x), interactions(x))

temp.x <- x
temp.x[0,] <- x[0,]
expect_equal(temp.x, x)
temp.x[,0] <- x[,0]
expect_equal(temp.x, x)

# Testing the combining.

xsub <- x[1:5,]
xsub2 <- x[6:20,]
expect_equal(rbind(xsub, xsub2), x)
xsub <- x[5:10,]
xsub2 <- x[1:3,]
expect_equal(rbind(xsub, xsub2), x[c(5:10, 1:3),])
expect_error(rbind(xsub, xsub2[,1:2]), "objects must have the same colnames")

xsub <- x[,1]
xsub2 <- x[,2:4]
expect_equal(cbind(xsub, xsub2), x)

first.cols <- 3 
other.cols <- 1:4
xsub <- x[,first.cols]
xsub2 <- x[,other.cols]
xsub.comb <- cbind(xsub, xsub2)
xsub.ref <- x[,c(first.cols, other.cols)]
colnames(xsub.comb) <- colnames(xsub.ref) # Keeping the column names happy again.
expect_equal(xsub.comb, xsub.ref)

expect_error(cbind(xsub, xsub2[1:10,]), "interactions must be identical in 'cbind'")

expect_identical(nrow(rbind(x[0,], x[0,])), 0L) # Behaviour with empties.
expect_identical(ncol(rbind(x[0,], x[0,])), ncol(x))
expect_equal(rbind(x, x[0,]), x)
expect_identical(nrow(cbind(x[,0], x[,0])), nrow(x))
expect_identical(ncol(cbind(x[,0], x[,0])), 0L)
expect_equal(cbind(x, x[,0]), x)

set.seed(1002)
next.starts <- round(runif(N, 1, 100))
next.ends <- next.starts + round(runif(N, 5, 20))
next.regions <- GRanges(rep(c("chrA", "chrB"), c(N-10, 10)), IRanges(next.starts, next.ends))

next.anchor1 <- sample(N, Np)
next.anchor2 <- sample(N, Np)
counts <- matrix(rpois(Np*Nlibs, lambda=10), ncol=Nlibs)
colnames(counts) <- colnames(x)
next.x <- InteractionSet(counts, GInteractions(next.anchor1, next.anchor2, next.regions))

c.x <- rbind(x, next.x)
expect_equivalent(assay(c.x), rbind(assay(x), assay(next.x)))
expect_identical(interactions(c.x), c(interactions(x), interactions(next.x)))

expect_identical(nrow(rbind(x[0,], next.x[0,])), 0L) # Behaviour with empties.
expect_identical(ncol(rbind(x[0,], next.x[0,])), ncol(x))
expect_identical(nrow(rbind(x, next.x[0,])), nrow(x)) # Not fully equal, as regions have changed.

# Testing the sorting.

o.x <- order(x)
expect_identical(o.x, order(interactions(x)))
expect_equal(sort(x), x[o.x,])
o.x2 <- order(x, next.x)
expect_identical(o.x2, order(interactions(x), interactions(next.x)))

expect_identical(duplicated(x), duplicated(interactions(x)))
temp.x <- rbind(x, x)    
expect_identical(duplicated(temp.x), duplicated(interactions(temp.x)))
expect_equal(x, unique(temp.x))

expect_identical(duplicated(x, fromLast=TRUE), duplicated(interactions(x), fromLast=TRUE))
expect_equal(x, unique(temp.x, fromLast=TRUE))

expect_identical(order(x[0,]), integer(0))
expect_identical(duplicated(x[0,]), logical(0))

# Testing the anchor swapping.

expect_identical(interactions(swapAnchors(x)), swapAnchors(interactions(x)))
expect_identical(interactions(swapAnchors(x, mode='reverse')), swapAnchors(interactions(x), mode='reverse'))
expect_identical(interactions(swapAnchors(x, mode='all')), swapAnchors(interactions(x), mode='all'))

# Testing the splitting.

flen <- c(5L, 10L, 5L)
f <- rep(1:3, flen)
out <- split(x, f)
expect_that(sapply(out, nrow), is_equivalent_to(flen))
for (i in seq_along(flen)) {
    expect_equal(out[[i]], x[f==i])
}

# Checking what happens with names.

temp.x <- x
ref.names <- paste0("X", seq_along(temp.x))
names(temp.x) <- ref.names
expect_identical(names(temp.x), ref.names)
expect_identical(names(temp.x), names(interactions(temp.x)))
expect_identical(names(temp.x[2:5]), ref.names[2:5])
expect_identical(names(temp.x[2:5]), names(interactions(temp.x[2:5])))
expect_identical(names(temp.x[2:5]), names(interactions(temp.x)[2:5]))

combined <- rbind(temp.x, temp.x)
expect_identical(names(combined), names(interactions(combined)))
expect_identical(names(combined), names(c(interactions(temp.x), interactions(temp.x))))
combined <- rbind(temp.x, x)
expect_identical(names(combined), c(ref.names, character(length(x))))
expect_identical(names(interactions(combined)), c(ref.names, character(length(x))))

# Testing what happens with strictness.

sx <- InteractionSet(counts, GInteractions(all.anchor1, all.anchor2, all.regions, mode="strict"))

expect_output(show(sx), "class: InteractionSet 
dim: 20 4 
metadata(0):
assays(1): ''
rownames: NULL
rowData names(0):
colnames(4): 1 2 3 4
colData names(0):
type: StrictGInteractions
regions: 30", 
fixed=TRUE)

expect_is(interactions(sx), "StrictGInteractions")
expect_is(interactions(sx[,1:2]), "StrictGInteractions")
expect_is(interactions(sx[1:3,]), "StrictGInteractions")

# Testing pairs.

expect_identical(pairs(x), pairs(interactions(x)))
expect_identical(pairs(x, id=TRUE), pairs(interactions(x), id=TRUE))
expect_identical(pairs(x, as.grlist=TRUE), pairs(interactions(x), as.grlist=TRUE))

