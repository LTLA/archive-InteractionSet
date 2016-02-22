# Tests the construction and manipulation of GInteractions objects.

set.seed(8000)
N <- 30
all.starts <- round(runif(N, 1, 100))
all.ends <- all.starts + round(runif(N, 5, 20))
all.regions <- GRanges(rep(c("chrA", "chrB"), c(N-10, 10)), IRanges(all.starts, all.ends))

Np <- 20
all.anchor1 <- sample(N, Np)
all.anchor2 <- sample(N, Np)
x <- GInteractions(all.anchor1, all.anchor2, all.regions)

another.x <- x[sample(Np)] 

ref.match <- function(x, y) {
    match(do.call(paste, c(anchors(x, id=TRUE), sep=".")), 
          do.call(paste, c(anchors(y, id=TRUE), sep=".")))
}

expect_identical(match(x, another.x), ref.match(x, another.x))
expect_identical(match(x, another.x[1:5]), ref.match(x, another.x[1:5]))
expect_identical(match(x[20:10], another.x), ref.match(x[20:10], another.x))

expect_identical(match(x[0], another.x), integer(0))
expect_identical(match(x, another.x[0]), rep(as.integer(NA), Np))
regions(another.x)$score <- 5
expect_error(match(x, another.x), "'regions' must be identical")

# Trying for ISet objects.

set.seed(8001)
iset <- InteractionSet(matrix(runif(Np), dimnames=list(NULL, 1)), x)
another.x <- x[sample(Np)] 
iset2 <- InteractionSet(matrix(runif(Np), dimnames=list(NULL, 1)), another.x)

expect_identical(match(iset, x), ref.match(iset, x))
expect_identical(match(iset, another.x), ref.match(iset, another.x))
expect_identical(match(iset2, x), ref.match(iset2, x))
expect_identical(match(x, iset), ref.match(x, iset))
expect_identical(match(another.x, iset), ref.match(another.x, iset))
expect_identical(match(x, iset2), ref.match(x, iset2))
expect_identical(match(iset, iset2), ref.match(iset, iset2))

expect_identical(match(iset[10:15], another.x), ref.match(iset[10:15], another.x))
expect_identical(match(another.x, iset[10:15]), ref.match(another.x, iset[10:15]))
expect_identical(match(iset, another.x[20:6]), ref.match(iset, another.x[20:6]))
expect_identical(match(another.x[20:6], iset), ref.match(another.x[20:6], iset))
expect_identical(match(iset, iset2[1:6,]), ref.match(iset, iset2[1:6,]))

# Trying to compare 'GInteractions' objects.

expect_identical(pcompare(x, another.x), ifelse(x@anchor1==another.x@anchor1, x@anchor2-another.x@anchor2, x@anchor1-another.x@anchor1))
sub.x <- x[3:12]
expect_identical(pcompare(sub.x, another.x), ifelse(sub.x@anchor1==another.x@anchor1, sub.x@anchor2-another.x@anchor2, sub.x@anchor1-another.x@anchor1))
expect_identical(pcompare(another.x, sub.x), -pcompare(sub.x, another.x))
expect_identical(pcompare(x[0], another.x[0]), integer(0))

old <- pcompare(x, another.x)
expect_identical(x==x, !logical(length(x)))
expect_identical(x!=x, logical(length(x)))
expect_identical(x==another.x, old==0L)

regions(another.x)$whee <- 1
expect_identical(pcompare(x, another.x), old) # This should be okay, as metadata is ignored.
regions(another.x) <- resize(regions(another.x), width(regions(another.x))*2L)
expect_error(pcompare(x, another.x), "'regions' must be identical")

sx <- as(swapAnchors(x), "StrictGInteractions")
expect_warning(sx==x, "comparison between GInteractions objects of different strictness")
rsx <- as(swapAnchors(x, mode="reverse"), "ReverseStrictGInteractions")
expect_warning(rsx==sx, "comparison between GInteractions objects of different strictness")
expect_warning(rsx==x, "comparison between GInteractions objects of different strictness")
expect_warning(match(rsx, x), "comparison between GInteractions objects of different strictness")

