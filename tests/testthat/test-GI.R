# Tests the construction and manipulation of GInteractions objects.

set.seed(7000)
N <- 30
all.starts <- round(runif(N, 1, 100))
all.ends <- all.starts + round(runif(N, 5, 20))
all.regions <- GRanges(rep(c("chrA", "chrB"), c(N-10, 10)), IRanges(all.starts, all.ends))

Np <- 20
all.anchor1 <- sample(N, Np)
all.anchor2 <- sample(N, Np)
x <- GInteractions(all.anchor1, all.anchor2, all.regions)

expect_output(show(x), "GInteractions object with 20 interactions and 0 metadata columns:
       seqnames1   ranges1     seqnames2   ranges2
           <Rle> <IRanges>         <Rle> <IRanges>
   [1]      chrA [82, 100] ---      chrB [67,  84]
   [2]      chrA [76,  95] ---      chrB [64,  78]
   [3]      chrA [87, 104] ---      chrA [ 3,  23]
   [4]      chrB [67,  84] ---      chrA [41,  59]
   [5]      chrA [14,  19] ---      chrA [87, 104]
   ...       ...       ... ...       ...       ...
  [16]      chrA [59,  72] ---      chrA [68,  78]
  [17]      chrA [84, 103] ---      chrB [94, 113]
  [18]      chrB [89,  97] ---      chrA [20,  33]
  [19]      chrA [41,  59] ---      chrA [84, 103]
  [20]      chrA [86, 105] ---      chrA [76,  95]
  -------
  regions: 30 ranges and 0 metadata columns
  seqinfo: 2 sequences from an unspecified genome; no seqlengths", fixed=TRUE)

# Testing all of our new slots:

expect_that(x, is_a("GInteractions"))
expect_true(!is.unsorted(regions(x)))

o <- order(all.regions)
new.regions <- all.regions[o]
expect_identical(regions(x), new.regions)

new.pos <- integer(length(o))
new.pos[o] <- seq_along(new.pos)
new.anchor1 <- new.pos[all.anchor1]
new.anchor2 <- new.pos[all.anchor2]
expect_identical(anchors(x, id=TRUE, type="first"), new.anchor1)
expect_identical(anchors(x, id=TRUE, type="second"), new.anchor2)
expect_identical(anchors(x, id=TRUE), list(first=new.anchor1, second=new.anchor2))

expect_identical(anchors(x, type="first"), new.regions[new.anchor1])
expect_identical(anchors(x, type="second"), new.regions[new.anchor2])
expect_identical(anchors(x), GRangesList(first=new.regions[new.anchor1], second=new.regions[new.anchor2]))

# Testing alternative construction methods for the GInteractions object:

x2 <- GInteractions(all.regions[all.anchor1], all.regions[all.anchor2])
was.used <- sort(unique(all.regions[union(all.anchor1, all.anchor2)])) # Only includes the regions actually used.
expect_identical(regions(x2), was.used)
expect_identical(anchors(x2), anchors(x))

x3 <- GInteractions(all.regions[all.anchor1], all.regions[all.anchor2], was.used)
expect_identical(anchors(x3, id=TRUE), anchors(x2, id=TRUE))
expect_identical(regions(x3), regions(x2))

anno.regions <- all.regions
anno.regions$score <- seq_along(anno.regions) # Seeing what happens with annotation.
anno.regions$revscore <- rev(seq_along(anno.regions))
x4 <- GInteractions(anno.regions[all.anchor1], anno.regions[all.anchor2])
expect_identical(regions(x2), regions(x4))
expect_identical(mcols(x4), DataFrame(anchor1.score=anno.regions$score[all.anchor1], anchor1.revscore=anno.regions$revscore[all.anchor1],
                                      anchor2.score=anno.regions$score[all.anchor2], anchor2.revscore=anno.regions$revscore[all.anchor2]))

# Testing with crappy inputs:

empty <- GInteractions(integer(0), numeric(0), GRanges())
expect_identical(length(empty), 0L)
expect_identical(length(anchors(empty, type="first")), 0L)
expect_identical(empty, GInteractions())

empty <- GInteractions(GRanges(), GRanges())
expect_identical(length(empty), 0L)
expect_identical(length(anchors(empty, type="first")), 0L)
expect_identical(empty, GInteractions())

empty <- GInteractions(GRanges(), GRanges(), all.regions)
expect_identical(length(empty), 0L)
expect_identical(length(anchors(empty, type="first")), 0L)
expect_identical(length(regions(empty)), length(all.regions))

expect_error(GInteractions(1:4, 1, all.regions), "first and second anchor vectors have different lengths")
expect_error(GInteractions(0:3, 1:4, all.regions), "all anchor indices must be positive integers")
expect_error(GInteractions(c(1,2,3,NA), 1:4, all.regions), "all anchor indices must be finite integers")
expect_error(GInteractions(c(1,2,3,length(all.regions)+1L), 1:4, all.regions), "all anchor indices must refer to entries in 'regions'")
missing.value <- GRanges("chrB", IRanges(1000,1000))
expect_error(GInteractions(missing.value, missing.value, all.regions), "anchor regions missing in specified 'regions'")

# Testing metadata inclusion:

m.x <- GInteractions(all.anchor1, all.anchor2, all.regions, score=1L)
expect_identical(m.x$score, rep(1L, length(m.x)))
m.x <- GInteractions(all.anchor1, all.anchor2, all.regions, score=seq_along(all.anchor1))
expect_identical(m.x$score, seq_along(m.x))

m.x <- GInteractions(all.regions[all.anchor1], all.regions[all.anchor2], all.regions, score=1L)
expect_identical(m.x$score, rep(1L, length(m.x)))
m.x <- GInteractions(all.regions[all.anchor1], all.regions[all.anchor2], all.regions, score=seq_along(all.anchor1))
expect_identical(m.x$score, seq_along(m.x))

m.x <- GInteractions(all.regions[all.anchor1], all.regions[all.anchor2], score=1L)
expect_identical(m.x$score, rep(1L, length(m.x)))
m.x <- GInteractions(all.regions[all.anchor1], all.regions[all.anchor2], score=seq_along(all.anchor1))
expect_identical(m.x$score, seq_along(m.x))

m.regions <- all.regions
m.regions$GC <- rev(seq_along(m.regions))
m.x <- GInteractions(m.regions[all.anchor1], m.regions[all.anchor2], score=seq_along(all.anchor1))
expect_identical(m.x$score, seq_along(m.x))
expect_identical(m.x$anchor1.GC, m.regions$GC[all.anchor1])
expect_identical(m.x$anchor2.GC, m.regions$GC[all.anchor2])

m.x <- GInteractions(integer(0), integer(0), all.regions, score=1)
expect_identical(m.x$score, numeric(0))
m.x <- GInteractions(regions=all.regions, score=numeric(0))
expect_identical(m.x$score, numeric(0))
m.x <- GInteractions(regions=all.regions, score=1)
expect_identical(m.x$score, numeric(0))
m.x <- GInteractions(regions=all.regions, score=1:2)
expect_identical(m.x$score, integer(0))

# Testing setters.

set.seed(7001)
shuffled <- sample(100, N, replace=TRUE)
regions(x)$score <- shuffled
expect_that(regions(x)$score, is_identical_to(shuffled))
expect_false(identical(regions(x), new.regions))
regions(x) <- new.regions # Restoring.
expect_true(identical(regions(x), new.regions))

fresh.anchor1 <- sample(N, Np)
fresh.anchor2 <- sample(N, Np)
anchors(x) <- list(fresh.anchor1, fresh.anchor2)
expect_identical(anchors(x, id=TRUE, type="first"), fresh.anchor1)
expect_identical(anchors(x, id=TRUE, type="second"), fresh.anchor2)
expect_error(anchors(x) <- list(new.anchor1, new.anchor2, new.anchor1), "must be a list of 2 numeric vectors")
expect_error(anchors(x) <- list(new.anchor1[1:(Np/2)], new.anchor2), "x@anchor2' is not parallel to 'x'")

mod.x <- x
anchors(mod.x, type="first") <- new.anchor1 # Checking that this also works
expect_identical(anchors(mod.x, id=TRUE, type="first"), new.anchor1)
mod.x <- x
anchors(mod.x, type="second") <- new.anchor2
expect_identical(anchors(mod.x, id=TRUE, type="second"), new.anchor2)
anchors(x, type="both") <- list(new.anchor1, new.anchor2) # Restoring.
expect_identical(anchors(x, id=TRUE, type="first"), new.anchor1)
expect_identical(anchors(x, id=TRUE, type="second"), new.anchor2)

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
expect_identical(anchors(x.dump2), anchors(x))
expect_identical(regions(x.dump), regions(x.dump2))

x.dump <- reduceRegions(x)
expect_identical(anchors(x), anchors(x.dump))
expect_identical(regions(x)[sort(unique(unlist(anchors(x, id=TRUE))))], regions(x.dump))

new.score <- runif(Np)
x$stuff <- new.score
expect_identical(x$stuff, mcols(x)$stuff)
expect_identical(colnames(mcols(x)), "stuff")
expect_output(show(x), "GInteractions object with 20 interactions and 1 metadata column:
       seqnames1   ranges1     seqnames2   ranges2 |              stuff
           <Rle> <IRanges>         <Rle> <IRanges> |          <numeric>
   [1]      chrA [82, 100] ---      chrB [67,  84] |  0.685012309812009
   [2]      chrA [76,  95] ---      chrB [64,  78] |  0.895816484699026
   [3]      chrA [87, 104] ---      chrA [ 3,  23] |  0.618890272220597
   [4]      chrB [67,  84] ---      chrA [41,  59] | 0.0507488863077015
   [5]      chrA [14,  19] ---      chrA [87, 104] |  0.621526781003922
   ...       ...       ... ...       ...       ... .                ...
  [16]      chrA [59,  72] ---      chrA [68,  78] |  0.727361923549324
  [17]      chrA [84, 103] ---      chrB [94, 113] |  0.402884092880413
  [18]      chrB [89,  97] ---      chrA [20,  33] |  0.906575692584738
  [19]      chrA [41,  59] ---      chrA [84, 103] |  0.582026617834345
  [20]      chrA [86, 105] ---      chrA [76,  95] |  0.863604818237945
  -------
  regions: 30 ranges and 0 metadata columns
  seqinfo: 2 sequences from an unspecified genome; no seqlengths", fixed=TRUE)
x$stuff <- NULL

new.si <- Seqinfo(seqnames=c("chrA", "chrB"), seqlengths=c(1000, 2000))
new.x <- x
seqinfo(new.x) <- new.si
expect_identical(seqinfo(new.x), new.si)

# Testing the subsetting.

rchosen <- 1:10
xsub <- x[rchosen,]
expect_output(show(xsub), "GInteractions object with 10 interactions and 0 metadata columns:
       seqnames1   ranges1     seqnames2   ranges2
           <Rle> <IRanges>         <Rle> <IRanges>
   [1]      chrA [82, 100] ---      chrB [67,  84]
   [2]      chrA [76,  95] ---      chrB [64,  78]
   [3]      chrA [87, 104] ---      chrA [ 3,  23]
   [4]      chrB [67,  84] ---      chrA [41,  59]
   [5]      chrA [14,  19] ---      chrA [87, 104]
   [6]      chrB [42,  54] ---      chrB [91,  98]
   [7]      chrB [64,  78] ---      chrA [55,  68]
   [8]      chrA [ 3,  23] ---      chrB [81,  98]
   [9]      chrA [61,  67] ---      chrA [46,  66]
  [10]      chrA [41,  49] ---      chrA [59,  72]
  -------
  regions: 30 ranges and 0 metadata columns
  seqinfo: 2 sequences from an unspecified genome; no seqlengths", fixed=TRUE)
expect_identical(xsub, x[rchosen])

log.chosen <- logical(length(x))
log.chosen[rchosen] <- TRUE
expect_identical(xsub, subset(x, log.chosen))
expect_identical(x[log.chosen], subset(x, log.chosen))

expect_identical(length(xsub), length(rchosen))
expect_identical(regions(xsub), regions(x))
expect_identical(anchors(xsub, type="first"), new.regions[new.anchor1][rchosen])
expect_identical(anchors(xsub, type="second"), new.regions[new.anchor2][rchosen])

temp.x <- x
temp.x$score <- new.score
expect_identical(temp.x[rchosen]$score, new.score[rchosen])
expect_identical(nrow(mcols(temp.x[rchosen])), length(rchosen))

expect_that(length(x[0,]), is_identical_to(0L))
expect_that(length(x[0]), is_identical_to(0L))
expect_error(x[,1], "invalid subsetting")

# Testing subset assignment.

temp.x <- x
temp.x[1:5+10,] <- x[1:5,]
new.index <- seq_along(x)
new.index[1:5+10] <- 1:5
expect_identical(anchors(temp.x, type="first"), anchors(x, type="first")[new.index,])
expect_identical(anchors(temp.x, type="second"), anchors(x, type="second")[new.index,])

temp.x <- x
temp.x[0,] <- x[0,]
expect_identical(temp.x, x)
temp.x[] <- x
expect_identical(temp.x, x)

temp.x <- x
temp.x$score <- new.score
temp.x[1:5]$score <- 5:1
mod.score <- new.score
mod.score[1:5] <- 5:1
expect_identical(temp.x$score, mod.score)
expect_identical(temp.x[rchosen]$score, mod.score[rchosen])
expect_identical(nrow(mcols(temp.x[rchosen])), length(rchosen))

# Testing the combining.

xsub <- x[1:5,]
xsub2 <- x[6:20,]
expect_identical(rbind(xsub, xsub2), x)
expect_identical(c(xsub, xsub2), x)
xsub <- x[1:15]
xsub2 <- x[16:20]
expect_identical(rbind(xsub, xsub2), x)
expect_identical(c(xsub, xsub2), x)

expect_identical(rbind(x[0,], x[0,]), x[0,])
expect_identical(rbind(x, x[0,]), x)

temp.x <- x
temp.x$score <- new.score
double.up <- rbind(temp.x, temp.x)
expect_identical(regions(double.up), regions(x))
expect_identical(anchors(double.up, type="first"), rep(anchors(x, type="first"), 2))
expect_identical(anchors(double.up, type="second"), rep(anchors(x, type="second"), 2))
expect_identical(double.up$score, rep(temp.x$score, 2))

# Combining with different anchor regions.

set.seed(7002)
next.starts <- round(runif(N, 1, 100))
next.ends <- next.starts + round(runif(N, 5, 20))
next.regions <- GRanges(rep(c("chrA", "chrB"), c(N-10, 10)), IRanges(next.starts, next.ends))

next.anchor1 <- sample(N, Np)
next.anchor2 <- sample(N, Np)
next.x <- GInteractions(next.anchor1, next.anchor2, next.regions)

c.x <- rbind(x, next.x)
expect_identical(c(anchors(x, type="first"), anchors(next.x, type="first")), anchors(c.x, type="first"))
expect_identical(c(anchors(x, type="second"), anchors(next.x, type="second")), anchors(c.x, type="second"))
expect_identical(unique(sort(c(regions(x), regions(next.x)))), regions(c.x))

expect_identical(anchors(rbind(x[0,], next.x[0,])), anchors(x[0,])) # Behaviour with empties.
expect_identical(anchors(rbind(x, next.x[0,])), anchors(x)) # Not fully equal, as regions have changed.

next.x2 <- GInteractions(1:10, 1:10, next.regions[1:10]) # What happens with non-equal lengths of the regions?
c.x <- rbind(x, next.x2)
expect_identical(c(anchors(x, type="first"), anchors(next.x2, type="first")), anchors(c.x, type="first"))
expect_identical(c(anchors(x, type="second"), anchors(next.x2, type="second")), anchors(c.x, type="second"))
expect_identical(unique(sort(c(regions(x), regions(next.x2)))), regions(c.x))

# Plus metadata

temp.x <- x
temp.x$score <- new.score
double.up <- rbind(temp.x, temp.x)
expect_identical(regions(double.up), regions(x))
expect_identical(anchors(double.up, type="first"), rep(anchors(x, type="first"), 2))
expect_identical(anchors(double.up, type="second"), rep(anchors(x, type="second"), 2))
expect_identical(double.up$score, rep(temp.x$score, 2))

# Testing the sorting.

o.x <- order(anchors(x, type="first"), anchors(x, type="second"))
expect_identical(o.x, order(x))
expect_identical(sort(x), x[o.x,])

o.x2 <- order(anchors(x, type="first"), anchors(x, type="second"), anchors(next.x, type="first"), anchors(next.x, type="second"))
expect_identical(o.x2, order(x, next.x))

is.dup <- duplicated(paste0(anchors(x, type="first"), ".", anchors(x, type="second")))
expect_identical(is.dup, duplicated(x))
temp.x <- rbind(x, x)    
is.dup <- duplicated(paste0(anchors(temp.x, type="first"), ".", anchors(temp.x, type="second")))
expect_identical(is.dup, duplicated(temp.x))
expect_true(all(tail(is.dup, length(x)))) # if ordering is stable; only the first occurrence should be true.
expect_identical(x, unique(temp.x))

is.dup <- duplicated(paste0(anchors(temp.x, type="first"), ".", anchors(temp.x, type="second")), fromLast=TRUE)
expect_identical(is.dup, duplicated(temp.x, fromLast=TRUE))
expect_true(all(head(is.dup, length(x)))) # if ordering is stable; only the first occurrence should be true.
expect_equal(x, unique(temp.x, fromLast=TRUE))
expect_false(any(duplicated(unique(temp.x))))

expect_identical(order(x[0,]), integer(0))
expect_identical(duplicated(x[0,]), logical(0))

# Testing the anchor swapping.

new.x <- swapAnchors(x)
expect_identical(anchors(new.x, type="first", id=TRUE), pmin(anchors(x, type="first", id=TRUE), anchors(x, type="second", id=TRUE)))
expect_identical(anchors(new.x, type="second", id=TRUE), pmax(anchors(x, type="first", id=TRUE), anchors(x, type="second", id=TRUE)))
expect_identical(regions(x), regions(new.x))

new.x2 <- swapAnchors(x, mode="reverse")
expect_identical(anchors(new.x2, type="first"), anchors(new.x, type="second"))
expect_identical(anchors(new.x, type="first"), anchors(new.x2, type="second"))
expect_identical(regions(x), regions(new.x2))

new.x3 <- swapAnchors(x, mode="all")
expect_identical(anchors(new.x3, type="first"), anchors(x, type="second"))
expect_identical(anchors(x, type="first"), anchors(new.x3, type="second"))
expect_identical(regions(x), regions(new.x3))

# Testing the splitting.

flen <- c(5L, 10L, 5L)
f <- rep(1:3, flen)
out <- split(x, f)
expect_equivalent(lengths(out), flen)
for (i in seq_along(flen)) {
    expect_identical(out[[i]], x[f==i])
}

temp.x <- x
temp.x$score <- new.score
out <- split(temp.x, f)
for (i in seq_along(flen)) {
    expect_identical(out[[i]]$score, temp.x[f==i]$score)
}

# Testing data.frame generation.

out <- as.data.frame(x)
ref1 <- as.data.frame(anchors(x, type="first"))
colnames(ref1) <- paste0(colnames(ref1), "1")
ref2 <- as.data.frame(anchors(x, type="second"))
colnames(ref2) <- paste0(colnames(ref2), "2")
expect_identical(out, data.frame(ref1, ref2))

temp.x <- x
temp.x$stuff <- new.score
out <- as.data.frame(temp.x)
expect_identical(out, data.frame(ref1, ref2, stuff=new.score))

names(temp.x) <- paste0("X", seq_along(temp.x))
out <- as.data.frame(temp.x)
expect_identical(out, data.frame(ref1, ref2, stuff=new.score, row.names=names(temp.x)))

empty <- as.data.frame(x[0,])
expect_identical(colnames(empty), colnames(as.data.frame(x)))
expect_identical(nrow(empty), 0L)

# Testing environment generation via 'with'

expect_identical(anchors(x, id=TRUE, type="first"), with(x, anchor1))
expect_identical(anchors(x, id=TRUE, type="second"), with(x, anchor2))
expect_identical(regions(x), with(x, regions))
expect_identical(anchors(x, type="first"), with(x, regions[anchor1]))
expect_identical(anchors(x, type="second"), with(x, regions[anchor2]))
expect_identical(names(x), with(x, names))

# Testing object flipping to a Pairs/GRangesList.

temp.x <- x
temp.x$score <- new.score
prs <- pairs(temp.x)
expect_identical(anchors(x, type="first"), first(prs))
expect_identical(anchors(x, type="second"), second(prs))
expect_identical(mcols(temp.x), mcols(prs))
expect_identical(names(temp.x), names(prs))
expect_identical(makeGInteractionsFromGRangesPairs(prs), reduceRegions(temp.x))                 

grl <- pairs(temp.x, as.grlist=TRUE)
first <- do.call(c, sapply(grl, function(x) { unname(x[1]) }))
second <- do.call(c, sapply(grl, function(x) { unname(x[2]) }))
expect_identical(anchors(x, type="first"), first)
expect_identical(anchors(x, type="second"), second)
expect_identical(mcols(temp.x), mcols(grl))

out <- pairs(temp.x, id=TRUE)
expect_identical(from(out), anchors(temp.x, id=TRUE, type="first"))
expect_identical(to(out), anchors(temp.x, id=TRUE, type="second"))
expect_identical(nnode(out), length(regions(temp.x)))

# Checking what happens with names.

temp.x <- x
ref.names <- paste0("X", seq_along(temp.x))
names(temp.x) <- ref.names
expect_output(show(temp.x),"GInteractions object with 20 interactions and 0 metadata columns:
      seqnames1   ranges1     seqnames2   ranges2
          <Rle> <IRanges>         <Rle> <IRanges>
   X1      chrA [82, 100] ---      chrB [67,  84]
   X2      chrA [76,  95] ---      chrB [64,  78]
   X3      chrA [87, 104] ---      chrA [ 3,  23]
   X4      chrB [67,  84] ---      chrA [41,  59]
   X5      chrA [14,  19] ---      chrA [87, 104]
  ...       ...       ... ...       ...       ...
  X16      chrA [59,  72] ---      chrA [68,  78]
  X17      chrA [84, 103] ---      chrB [94, 113]
  X18      chrB [89,  97] ---      chrA [20,  33]
  X19      chrA [41,  59] ---      chrA [84, 103]
  X20      chrA [86, 105] ---      chrA [76,  95]
  -------
  regions: 30 ranges and 0 metadata columns
  seqinfo: 2 sequences from an unspecified genome; no seqlengths", fixed=TRUE)

expect_identical(names(temp.x), ref.names)
expect_identical(names(temp.x[2:5]), ref.names[2:5])
expect_identical(names(c(temp.x, temp.x)), c(ref.names, ref.names))
expect_identical(names(c(temp.x, x)), c(ref.names, character(length(x))))

for (id in c(TRUE, FALSE)) {
    expect_identical(names(anchors(temp.x, id=id)[[1]]), ref.names)
    expect_identical(names(anchors(temp.x, id=id)[[2]]), ref.names)
    expect_identical(names(anchors(temp.x, id=id, type="first")), ref.names)
    expect_identical(names(anchors(temp.x, id=id, type="second")), ref.names)
}

# Testing strictness.

sx <- GInteractions(all.anchor1, all.anchor2, all.regions, mode="strict")
expect_that(sx, is_a("StrictGInteractions"))
expect_identical(sx, as(x, "StrictGInteractions"))
expect_identical(anchors(sx, id=TRUE, type="first"), do.call(pmin, anchors(x, id=TRUE)))
expect_identical(anchors(sx, id=TRUE, type="second"), do.call(pmax, anchors(x, id=TRUE)))
expect_identical(regions(sx), regions(x))

temp.sx <- sx
anchors(temp.sx, type="first") <- fresh.anchor1
expect_identical(anchors(temp.sx, id=TRUE, type="first"), pmin(fresh.anchor1, anchors(sx, id=TRUE, type="second")))
expect_identical(anchors(temp.sx, id=TRUE, type="second"), pmax(fresh.anchor1, anchors(sx, id=TRUE, type="second")))

temp.sx <- sx
anchors(temp.sx, type="second") <- fresh.anchor2
expect_identical(anchors(temp.sx, id=TRUE, type="first"), pmin(fresh.anchor2, anchors(sx, id=TRUE, type="first")))
expect_identical(anchors(temp.sx, id=TRUE, type="second"), pmax(fresh.anchor2, anchors(sx, id=TRUE, type="first")))

temp.sx <- sx
anchors(temp.sx) <- list(fresh.anchor1, fresh.anchor2)
expect_identical(anchors(temp.sx, id=TRUE, type="first"), pmin(fresh.anchor1, fresh.anchor2))
expect_identical(anchors(temp.sx, id=TRUE, type="second"), pmax(fresh.anchor1, fresh.anchor2))

temp.sx2 <- sx
anchors(temp.sx2, type="both") <- list(fresh.anchor2, fresh.anchor1)
expect_identical(temp.sx2, temp.sx)

# Testing reverse strictness.

rsx <- GInteractions(all.anchor1, all.anchor2, all.regions, mode="reverse")
expect_that(rsx, is_a("ReverseStrictGInteractions"))
expect_identical(rsx, as(x, "ReverseStrictGInteractions"))
expect_identical(rsx, as(sx, "ReverseStrictGInteractions"))
expect_identical(sx, as(rsx, "StrictGInteractions"))
expect_identical(anchors(rsx, id=TRUE, type="first"), do.call(pmax, anchors(x, id=TRUE)))
expect_identical(anchors(rsx, id=TRUE, type="second"), do.call(pmin, anchors(x, id=TRUE)))
expect_identical(regions(rsx), regions(x))

temp.rsx <- rsx
anchors(temp.rsx, type="first") <- fresh.anchor1
expect_identical(anchors(temp.rsx, id=TRUE, type="first"), pmax(fresh.anchor1, anchors(rsx, id=TRUE, type="second")))
expect_identical(anchors(temp.rsx, id=TRUE, type="second"), pmin(fresh.anchor1, anchors(rsx, id=TRUE, type="second")))

temp.rsx <- rsx
anchors(temp.rsx, type="second") <- fresh.anchor2
expect_identical(anchors(temp.rsx, id=TRUE, type="first"), pmax(fresh.anchor2, anchors(rsx, id=TRUE, type="first")))
expect_identical(anchors(temp.rsx, id=TRUE, type="second"), pmin(fresh.anchor2, anchors(rsx, id=TRUE, type="first")))

temp.rsx <- rsx
anchors(temp.rsx) <- list(fresh.anchor1, fresh.anchor2)
expect_identical(anchors(temp.rsx, id=TRUE, type="first"), pmax(fresh.anchor1, fresh.anchor2))
expect_identical(anchors(temp.rsx, id=TRUE, type="second"), pmin(fresh.anchor1, fresh.anchor2))

temp.rsx2 <- rsx
anchors(temp.rsx2, type="both") <- list(fresh.anchor2, fresh.anchor1)
expect_identical(temp.rsx2, temp.rsx)

# Testing r'binding of objects of different strictness.

expect_identical(c(rsx, sx, x), as(c(x, x, x), "ReverseStrictGInteractions"))
expect_identical(c(x, sx, rsx), as(c(x, swapAnchors(x), swapAnchors(x, mode="reverse")), "GInteractions"))
expect_identical(c(sx, rsx, x), as(c(x, x, x), "StrictGInteractions"))

