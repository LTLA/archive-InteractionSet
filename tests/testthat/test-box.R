# Testing the bounding box behaviour.

set.seed(10000)
N <- 30
all.starts <- round(runif(N, 1, 100))
all.ends <- all.starts + round(runif(N, 5, 20))
all.regions <- GRanges(rep(c("chrA", "chrB"), c(N-10, 10)), IRanges(all.starts, all.ends))

Np <- 20
all.anchor1 <- sample(N, Np)
all.anchor2 <- sample(N, Np)
x <- GInteractions(all.anchor1, all.anchor2, all.regions)

# Making up a sensible grouping.

all.chrs <- as.character(seqnames(regions(x)))
for (i in 1:3) { 
    f <- paste0(all.chrs[anchors(x, type="first", id=TRUE)],
                all.chrs[anchors(x, type="second", id=TRUE)],
                sample(i, Np, replace=TRUE))

    ref1 <- unlist(range(split(anchors(x, type="first"), f)))
    ref2 <- unlist(range(split(anchors(x, type="second"), f)))
    ref <- GInteractions(unname(ref1), unname(ref2))
    names(ref) <- names(ref1)
    expect_identical(ref, boundingBox(x, f))
}

# Checking that it works properly when 'f' is not specified.

only.A <- all.chrs[anchors(x, type="first", id=TRUE)] == "chrA" & all.chrs[anchors(x, type="second", id=TRUE)] == "chrA"
x.A <- x[only.A]
ref1 <- unlist(range(anchors(x.A, type="first")))
ref2 <- unlist(range(anchors(x.A, type="second")))
ref <- GInteractions(unname(ref1), unname(ref2))
names(ref) <- 1
expect_identical(boundingBox(x.A), ref)

# Breaking it with silly inputs.

expect_error(boundingBox(x), "multiple chromosomes for group '1'")
f <- rep("whee", Np)
expect_error(boundingBox(x,f), "multiple chromosomes for group 'whee'")
ref <- GInteractions(ref1[0], ref2[0])
names(ref) <- character(0)
expect_identical(boundingBox(x[0]), ref)
