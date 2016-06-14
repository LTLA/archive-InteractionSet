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
colnames(counts) <- seq_len(Nlibs)
offs <- matrix(rnorm(Np*Nlibs), ncol=Nlibs)
x <- InteractionSet(list(counts, offs), GInteractions(all.anchor1, all.anchor2, all.regions))

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

ref.fun <- function(x, rows, cols, fill, ass=1, sam=1, swap=TRUE) { # Slow and steady implementation.
    all.anchors <- anchors(x, id=TRUE)
    if (missing(fill)) { fill <- assay(x, ass)[,sam] }
    ref <- Matrix::Matrix(as(NA, typeof(fill)), length(rows), length(cols))
    for (i in seq_len(nrow(x))){ 
        a1 <- all.anchors$first[i]
        a2 <- all.anchors$second[i]
        ref[rows==a1,cols==a2] <- fill[i]
        if (swap) { ref[rows==a2,cols==a1] <- fill[i] }
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

out2 <- inflate(x, chosen.rows, chosen.cols, fill=blah, sparse=TRUE) # Trying out a sparse matrix.
expect_is(as.matrix(out2), "dgCMatrix")
expect_equal(dim(out), dim(out2))
ref <- as.matrix(out)
not.missing <- Matrix::which(!is.na(ref))
expect_equal(as.matrix(out2)[not.missing], ref[not.missing])
expect_true(all(as.matrix(as.matrix(out2))[!not.missing]==0))

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

out2 <- inflate(x, chosen.rows, chosen.cols, sparse=TRUE) # Trying out a sparse matrix, again.
expect_is(as.matrix(out2), "dgCMatrix")
expect_equal(dim(out), dim(out2))
ref <- as.matrix(out)
not.missing <- Matrix::which(!is.na(ref))
expect_equal(as.matrix(out2)[not.missing], ref[not.missing])
expect_true(all(as.matrix(as.matrix(out2))[!not.missing]==0))

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

out <- inflate(x, chosen.rows, chosen.cols, fill=blah, swap=FALSE) # Symmetric space, so there's guaranteed to swapping.
expect_equal(as.matrix(out), ref.fun(x, chosen.rows, chosen.cols, fill=blah, swap=FALSE))

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
keep.x <- subsetByOverlaps(x, GInteractions(all.chr[1], all.chr[1])) 
keep.x <- sort(swapAnchors(keep.x))
expect_identical(anchors(x2), anchors(keep.x))
expect_equal(assay(x2)[,1], assay(keep.x)[,1]) # Not identical, due to coercion to double.

# What happens when you turn off uniqueness (in this case, we have symmetry):
x2 <- deflate(y, collapse=FALSE)
x2 <- sort(x2)
not.diag <- anchors(keep.x, type="first", id=TRUE)!=anchors(keep.x, type="second", id=TRUE)
keep.x <- rbind(keep.x[not.diag], swapAnchors(keep.x[not.diag], mode='all'), keep.x[!not.diag])
keep.x <- sort(keep.x)
expect_identical(anchors(x2), anchors(keep.x))
expect_equal(assay(x2)[,1], assay(keep.x)[,1])

# Behaviour for different index sets:
y <- inflate(x, "chrA", "chrB")
x2 <- deflate(y)
x2 <- sort(x2)
keep.x <- subsetByOverlaps(x, GInteractions(all.chr[1], all.chr[2])) 
keep.x <- sort(swapAnchors(keep.x))
expect_identical(anchors(x2), anchors(keep.x))
expect_equal(assay(x2)[,1], assay(keep.x)[,1])

# Deflating a sparseMatrix (and other options).
y <- inflate(x, "chrA", "chrA", sparse=TRUE)
xref <- inflate(x, "chrA", "chrA")
ref <- as.matrix(xref)
ref[is.na(ref)] <- 0
expect_identical(as.matrix(y), as(ref, "dgCMatrix"))
expect_equal(deflate(y), deflate(xref, use.zero=FALSE, use.na=FALSE)) # Getting rid of genuine zeros in there.

na.deflated <- deflate(xref, use.na=TRUE)
expect_identical(length(na.deflated), sum(seq_len(nrow(ref))))
xref2 <- xref
as.matrix(xref2) <- ref
expect_equal(deflate(y, use.zero=TRUE), deflate(xref2, use.zero=TRUE))

ex <- !is.na(as.matrix(xref))
alt <- deflate(xref, use.na=FALSE, use.zero=TRUE)
expect_equal(deflate(xref, extract=ex), alt)
expect_false(isTRUE(all.equal(alt, deflate(xref, use.na=TRUE)))) # definitely different
expect_equal(deflate(xref, extract=ex, use.na=TRUE), alt) # ... but 'ex' overrides it.

# Trying out some silliness.

expect_true(nrow(deflate(ContactMatrix(matrix(0, 4, 0), 1:4, integer(0), all.regions)))==0L)
expect_true(nrow(deflate(ContactMatrix(matrix(0, 0, 4), integer(0), 1:4, all.regions)))==0L)

