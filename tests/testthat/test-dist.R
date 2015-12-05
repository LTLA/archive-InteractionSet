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
colnames(counts) <- seq_len(Nlibs)
x <- InteractionSet(counts, GInteractions(all.anchor1, all.anchor2, all.regions))

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

ax1 <- anchors(x, type="first", id=TRUE)
ax2 <- anchors(x, type="second", id=TRUE)
expect_that(pairdist(x, type="diag"), is_identical_to(ifelse(is.intra, pmax(ax1, ax2) - pmin(ax1, ax2), as.integer(NA))))

# What happens with empty inputs?

expect_that(pairdist(x[0,]), is_identical_to(integer(0)))
expect_that(pairdist(x[0,], type="intra"), is_identical_to(logical(0)))
expect_that(pairdist(x[!is.intra,]), is_identical_to(rep(as.integer(NA), sum(!is.intra))))

####################################################
# Testing distance calculation for a ContactMatrix.

set.seed(101)
Nr <- 100
Nc <- 200
all.anchor1 <- sample(N, Nr, replace=TRUE)
all.anchor2 <- sample(N, Nc, replace=TRUE)
counts <- matrix(rpois(Nr*Nc, lambda=10), Nr, Nc)
x <- ContactMatrix(counts, all.anchor1, all.anchor2, all.regions)

a1 <- anchors(x, type="row")
a2 <- anchors(x, type="column")
is.intra <- intrachr(x) 
expect_that(is.intra, is_identical_to(matrix(outer(as.character(seqnames(a1)), as.character(seqnames(a2)), "=="), nrow=Nr, ncol=Nc)))
expect_that(is.intra, is_identical_to(pairdist(x, type="intra")))

ref <- abs(outer(start(a1)+end(a1), start(a2)+end(a2), "-"))
ref <- ref/2L
ref[!is.intra] <- NA
storage.mode(ref) <- "integer"
expect_that(pairdist(x), is_identical_to(ref))
                                                
ref <- abs(outer(anchors(x, type="row", id=TRUE), anchors(x, type="column", id=TRUE), "-"))
ref[!is.intra] <- NA
storage.mode(ref) <- "integer"
expect_that(pairdist(x, type="diag"), is_identical_to(ref))

as1 <- start(a1)
as2 <- start(a2)
ae1 <- end(a1)
ae2 <- end(a2)
ref <- outer(as1, as2, pmax) - outer(ae1, ae2, pmin) - 1L
ref[!is.intra] <- NA                    
expect_that(pairdist(x,type="gap"), is_identical_to(ref))

ref <- outer(ae1, ae2, pmax) - outer(as1, as2, pmin) + 1L
ref[!is.intra] <- NA                    
expect_that(pairdist(x,type="span"), is_identical_to(ref))

# Don't need to worry about NA's in 'x@matrix', they don't affect the calculation.
# Still checking silly inputs, though:

expect_that(pairdist(x[0,]), is_identical_to(matrix(0L, 0, Nc)))
expect_that(pairdist(x[,0], type="intra"), is_identical_to(matrix(FALSE, Nr, 0)))
is.ra <- which(seqnames(anchors(x, type="row"))=="chrA")
is.cb <- which(seqnames(anchors(x, type="column"))=="chrB")
expect_that(pairdist(x[is.ra,is.cb]), is_identical_to(matrix(as.integer(NA), length(is.ra), length(is.cb))))


