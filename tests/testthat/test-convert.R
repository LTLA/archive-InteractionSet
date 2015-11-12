require(InteractionSet); require(testthat)

# This tests the conversion functions of an InteractionSet object.

set.seed(6000)
N <- 30
all.starts <- round(runif(N, 1, 100))
all.ends <- all.starts + round(runif(N, 5, 20))
all.regions <- GRanges(rep(c("chrA", "chrB"), c(N-10, 10)), IRanges(all.starts, all.ends))

Np <- 20
all.anchor1 <- sample(N, Np)
all.anchor2 <- sample(N, Np)
Nlibs <- 4
counts <- matrix(rpois(Np*Nlibs, lambda=10), ncol=Nlibs)
offs <- matrix(rnorm(Np*Nlibs), ncol=Nlibs)
x <- InteractionSet(list(counts, offs), all.anchor1, all.anchor2, all.regions)

##########################################
# Standard construction with integers:

chosen.rows <- 1:10
chosen.cols <- 11:15
out <- inflate(x, chosen.rows, chosen.cols)
expect_identical(regions(out), regions(x))
expect_identical(anchors(out, type="row"), regions(x)[chosen.rows])
expect_identical(anchors(out, type="column"), regions(x)[chosen.cols])
expect_identical(anchors(out, type="row", id=TRUE), chosen.rows)
expect_identical(anchors(out, type="column", id=TRUE), chosen.cols)

ref.fun <- function(x, rows, cols, fill, ass=1, sam=1) { # Slow and steady implementation.
    all.anchors <- anchors(x, id=TRUE)
    if (missing(fill)) { fill <- assay(x, ass)[,sam] }
    ref <- matrix(NA, length(rows), length(cols))
    for (i in seq_len(nrow(x))){ 
        a1 <- all.anchors$first[i]
        a2 <- all.anchors$second[i]
        ref[rows==a1,cols==a2] <- fill[i]
        ref[rows==a2,cols==a1] <- fill[i]
    }
    return(ref)
}

expect_equal(as.matrix(out), ref.fun(x, chosen.rows, chosen.cols))
out <- inflate(x, chosen.rows, chosen.cols, assay=2)
expect_equal(as.matrix(out), ref.fun(x, chosen.rows, chosen.cols, ass=2))
out <- inflate(x, chosen.rows, chosen.cols, sample=4)
expect_equal(as.matrix(out), ref.fun(x, chosen.rows, chosen.cols, sam=4))
blah <- runif(Np)
out <- inflate(x, chosen.rows, chosen.cols, fill=blah)
expect_equal(as.matrix(out), ref.fun(x, chosen.rows, chosen.cols, fill=blah))

# Dealing with duplication and resorting:
chosen.rows <- c(1:10, 1:10)
chosen.cols <- c(11:15, 11:15)
out <- inflate(x, chosen.rows, chosen.cols)
expect_identical(anchors(out, type="row", id=TRUE), chosen.rows)
expect_identical(anchors(out, type="column", id=TRUE), chosen.cols)
expect_equal(as.matrix(out), ref.fun(x, chosen.rows, chosen.cols))

chosen.rows <- as.integer(c(1,3,2,6,7,9,2,2,1))
chosen.cols <- as.integer(c(11,16,2,2,5))
out <- inflate(x, chosen.rows, chosen.cols)
expect_identical(anchors(out, type="row", id=TRUE), chosen.rows)
expect_identical(anchors(out, type="column", id=TRUE), chosen.cols)
expect_equal(as.matrix(out), ref.fun(x, chosen.rows, chosen.cols))

# What happens with silly inputs?
expect_true(nrow(inflate(x, integer(0), 1:10))==0L)
expect_true(ncol(inflate(x, 1:5, integer(0)))==0L)
expect_error(inflate(x, 0, 1:10), "positive integer")
expect_error(inflate(x, as.numeric(NA), 1:10), "positive integer")
expect_error(inflate(x, 10000, 1:10), "positive integer")

##########################################
# Construction with character vectors.

out <- inflate(x, "chrA", "chrA")
chosen.rows <- which(seqnames(regions(out))=="chrA")
chosen.cols <- which(seqnames(regions(out))=="chrA")
expect_identical(anchors(out, type="row", id=TRUE), chosen.rows)
expect_identical(anchors(out, type="column", id=TRUE), chosen.cols)
expect_equal(as.matrix(out), ref.fun(x, chosen.rows, chosen.cols))

out <- inflate(x, "chrA", "chrB")
chosen.rows <- which(seqnames(regions(out))=="chrA")
chosen.cols <- which(seqnames(regions(out))=="chrB")
expect_identical(anchors(out, type="row", id=TRUE), chosen.rows)
expect_identical(anchors(out, type="column", id=TRUE), chosen.cols)
expect_equal(as.matrix(out), ref.fun(x, chosen.rows, chosen.cols))

out <- inflate(x, "chrA", c("chrA", "chrB")) # Multiple chromosomes.
chosen.rows <- which(seqnames(regions(out))=="chrA")
chosen.cols <- which(seqnames(regions(out)) %in%  c("chrA", "chrB"))
expect_identical(anchors(out, type="row", id=TRUE), chosen.rows)
expect_identical(anchors(out, type="column", id=TRUE), chosen.cols)
expect_equal(as.matrix(out), ref.fun(x, chosen.rows, chosen.cols))

expect_true(nrow(inflate(x, "whee", 1:10))==0L)
expect_true(ncol(inflate(x, 1:5, "whee"))==0L)

##########################################
# Construction with GRanges.

of.interest <- GRanges(c("chrA", "chrB"), IRanges(c(1, 10), c(20, 50)))

out <- inflate(x, of.interest, of.interest)
chosen.rows <- chosen.cols <- which(overlapsAny(regions(x), of.interest))
expect_identical(anchors(out, type="row", id=TRUE), chosen.rows)
expect_identical(anchors(out, type="column", id=TRUE), chosen.cols)
expect_equal(as.matrix(out), ref.fun(x, chosen.rows, chosen.cols))

out <- inflate(x, of.interest, of.interest, type="within")
chosen.rows <- chosen.cols <- which(overlapsAny(regions(x), of.interest, type="within"))
expect_identical(anchors(out, type="row", id=TRUE), chosen.rows)
expect_identical(anchors(out, type="column", id=TRUE), chosen.cols)
expect_equal(as.matrix(out), ref.fun(x, chosen.rows, chosen.cols))

out <- inflate(x, of.interest[1], of.interest[2], type="within")
chosen.rows <- which(overlapsAny(regions(x), of.interest[1], type="within"))
chosen.cols <- which(overlapsAny(regions(x), of.interest[2], type="within"))
expect_identical(anchors(out, type="row", id=TRUE), chosen.rows)
expect_identical(anchors(out, type="column", id=TRUE), chosen.cols)
expect_equal(as.matrix(out), ref.fun(x, chosen.rows, chosen.cols))

expect_true(nrow(inflate(x, GRanges(), 1:10))==0L)
expect_true(ncol(inflate(x, 1:5, GRanges()))==0L)
out.of.range <- GRanges("chrC", IRanges(1, 1))
expect_true(nrow(suppressWarnings(inflate(x, out.of.range, 1:10)))==0L)
expect_true(ncol(suppressWarnings(inflate(x, 1:5, out.of.range)))==0L)

all.chr <- range(all.regions)
expect_identical(inflate(x, all.chr[1], all.chr[2]), inflate(x, "chrA", "chrB"))

##########################################
# Deflation tests

y <- inflate(x, "chrA", "chrA")
x2 <- deflate(y)
x2 <- sort(x2)
keep.x <- subsetByOverlaps(x, GRangesList(all.chr[1], all.chr[1])) 
keep.x <- sort(keep.x)
expect_identical(anchors(x2), anchors(keep.x))
expect_identical(assay(x2)[,1], assay(keep.x)[,1])

# What happens when you turn off uniqueness (in this case, we have symmetry):
x2 <- deflate(y, unique=FALSE)
x2 <- sort(x2)
not.diag <- anchors(keep.x, type="first", id=TRUE)!=anchors(keep.x, type="second", id=TRUE)
keep.x <- rbind(keep.x[not.diag], keep.x[not.diag], keep.x[!not.diag])
keep.x <- sort(keep.x)
expect_identical(anchors(x2), anchors(keep.x))
expect_identical(assay(x2)[,1], assay(keep.x)[,1])

# Behaviour for different index sets:
y <- inflate(x, "chrA", "chrB")
x2 <- deflate(y)
x2 <- sort(x2)
keep.x <- subsetByOverlaps(x, GRangesList(all.chr[1], all.chr[2])) 
keep.x <- sort(keep.x)
expect_identical(anchors(x2), anchors(keep.x))
expect_identical(assay(x2)[,1], assay(keep.x)[,1])

expect_true(nrow(deflate(ContactMatrix(matrix(0, 4, 0), 1:4, integer(0), all.regions)))==0L)
expect_true(nrow(deflate(ContactMatrix(matrix(0, 0, 4), integer(0), 1:4, all.regions)))==0L)

