# Tests the distance calculation methods for an InteractionSet.

set.seed(100)
N <- 50
all.starts <- round(runif(N, 1, 100))
all.ends <- all.starts + round(runif(N, 5, 20))
all.regions <- GRanges(rep(c("chrA", "chrB"), c(N-10, 10)), IRanges(all.starts, all.ends))

Np <- 100
all.anchor1 <- sample(N, Np, replace=TRUE)
all.anchor2 <- sample(N, Np, replace=TRUE)
Nlibs <- 4
counts <- matrix(rpois(Np*Nlibs, lambda=10), ncol=Nlibs)
x <- InteractionSet(counts, all.anchor1, all.anchor2, all.regions)

a1 <- all.regions[all.anchor1]
a2 <- all.regions[all.anchor2]
swap <- a1 < a2
temp <- a2[swap]
a2[swap] <- a1[swap]
a1[swap] <- temp

is.intra <- intrachr(x)
expect_that(is.intra, is_identical_to(as.logical(seqnames(a1)==seqnames(a2))))
expect_that(is.intra, is_identical_to(pairdist(x, type="intra")))
expect_that(pairdist(x), is_identical_to(ifelse(is.intra, as.integer(abs(start(a1)+end(a1)-start(a2)-end(a2))/2L), as.integer(NA)))) # Don't use 'mid', it does its own rounding.
expect_that(pairdist(x,type="gap"), is_identical_to(ifelse(is.intra, pmax(start(a1), start(a2)) - pmin(end(a1), end(a2)) -1L, as.integer(NA)))) 
expect_that(pairdist(x,type="span"), is_identical_to(ifelse(is.intra, pmax(end(a1), end(a2)) - pmin(start(a1), start(a2)) +1L, as.integer(NA)))) 
expect_that(pairdist(x, type="diag"), is_identical_to(ifelse(is.intra, anchors(x, type="first", id=TRUE) - anchors(x, type="second", id=TRUE), as.integer(NA))))

# What happens with empty inputs?

expect_that(pairdist(x[0,]), is_identical_to(integer(0)))
expect_that(pairdist(x[0,], type="intra"), is_identical_to(logical(0)))
expect_that(pairdist(x[!is.intra,]), is_identical_to(rep(as.integer(NA), sum(!is.intra))))

