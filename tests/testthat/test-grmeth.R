# Tests the trimming of GInteractions objects.

set.seed(10000)

N <- 30
all.starts <- round(runif(N, 1, 100))
all.ends <- all.starts + round(runif(N, 5, 20))
all.regions <- GRanges(rep(c("chrA", "chrB"), c(N-10, 10)), IRanges(all.starts, all.ends))

Np <- 20
all.anchor1 <- sample(N, Np)
all.anchor2 <- sample(N, Np)
x <- GInteractions(all.anchor1, all.anchor2, all.regions)

suppressWarnings(seqlengths(x) <- c(chrA=50, chrB=100))
x2 <- trim(x)
ref <- trim(regions(x))
expect_identical(anchors(x2, type="first"), ref[anchors(x, type="first", id=TRUE)])
expect_identical(anchors(x2, type="second"), ref[anchors(x, type="second", id=TRUE)])

suppressWarnings(regions(x2) <- resize(regions(x2), fix="center", width=1000)) 
x2 <- trim(x2) # reduces down to one range per chromosome.
expect_identical(length(regions(x2)), length(seqlengths(x2)))

iset <- InteractionSet(matrix(0, nrow=Np, ncol=2), x)
iset2 <- trim(iset)
expect_identical(class(iset2), class(iset))
expect_identical(interactions(iset2), trim(x))

# Same again for ContactMatrix.

y <- inflate(x, "chrA", "chrB", fill=0)
y2 <- trim(y)
expect_identical(class(y), class(y2))
expect_identical(anchors(y2, type="row"), ref[anchors(y, type="row", id=TRUE)])
expect_identical(anchors(y2, type="column"), ref[anchors(y, type="column", id=TRUE)])

suppressWarnings(regions(y2) <- resize(regions(y2), fix="center", width=1000)) 
y2 <- trim(y2) 
expect_identical(length(regions(y2)), length(seqlengths(y2)))

# Tests width calculations for everyone.

x2 <- trim(x)
expect_identical(width(x2), DataFrame(anchor1=width(anchors(x2, type="first")), anchor2=width(anchors(x2, type="second"))))
iset2 <- trim(iset)
expect_identical(width(iset2), width(x2))
y2 <- trim(y)
expect_identical(width(y2), list(anchor1=width(anchors(y2, type="row")), anchor2=width(anchors(y2, type="column"))))


