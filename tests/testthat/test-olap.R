# Testing the various overlap methods for InteractionSet objects.
# This implicitly tests the GInteractions methods, because one calls the other.

set.seed(300)
N <- 50
all.starts <- round(runif(N, 1, 100))
all.ends <- all.starts + round(runif(N, 5, 20))
all.regions <- GRanges(rep(c("chrA", "chrB"), c(N-20, 20)), IRanges(all.starts, all.ends))
all.regions <- unique(all.regions)
N <- length(all.regions)

Np <- 100
all.anchor1 <- sample(N, Np, replace=TRUE)
all.anchor2 <- sample(N, Np, replace=TRUE)
Nlibs <- 4
counts <- matrix(rpois(Np*Nlibs, lambda=10), ncol=Nlibs)
colnames(counts) <- seq_len(Nlibs)
x <- InteractionSet(counts, GInteractions(all.anchor1, all.anchor2, all.regions))

#######################################################
# 1D overlaps with Vectors

for (object in c("GRanges", "GRangesList")) {

if (object=="GRanges") {
# GRanges.
set.seed(301)
Nq <- 10
query.starts <- round(runif(Nq, 1, 100))
query.ends <- query.starts + round(runif(Nq, 5, 20))
query.regions <- GRanges(rep(c("chrA", "chrB"), Nq/2), IRanges(query.starts, query.ends))

} else if (object=="GRangesList") {
# Paired regions in a GRangesList (don't be fooled by the GInteractions below!)
set.seed(302)
Nq2 <- 100
query.starts <- round(runif(Nq2, 1, 100))
query.ends <- query.starts + round(runif(Nq2, 5, 20))
query.regions1 <- GRanges(rep(c("chrA", "chrB"), Nq2/2), IRanges(query.starts, query.ends))
query.starts <- round(runif(Nq2, 1, 100))
query.ends <- query.starts + round(runif(Nq2, 5, 20))
query.regions2 <- GRanges(rep(c("chrA", "chrB"), Nq2/2), IRanges(query.starts, query.ends))
pairing <- pairs(GInteractions(query.regions1, query.regions2), as.grlist=TRUE)

}

for (param in seq_len(6)) {
    type <- "any"
    maxgap <- 0L
    minoverlap <- 1L
    use.region <- "both"
    if (param==2L) { 
        maxgap <- 10L
    } else if (param==3L) {
        minoverlap <- 10L
    } else if (param==4L) {
        type <- "within"
    } else if (param==5L) {
        use.region <- "first"
    } else if(param==6L) {
        use.region <- "second"
    }

    # Overlapping with GRanges.
    expected1 <- findOverlaps(anchors(x, type="first"), query.regions, type=type, maxgap=maxgap, minoverlap=minoverlap)
    expected2 <- findOverlaps(anchors(x, type="second"), query.regions, type=type, maxgap=maxgap, minoverlap=minoverlap)
    if (use.region=="both") {
        ref <- Hits(c(queryHits(expected1), queryHits(expected2)), c(subjectHits(expected1), subjectHits(expected2)),
                    nLnode=length(x), nRnode=Nq, sort.by.query=TRUE)
    } else if (use.region=="first") {
        ref <- expected1
    } else if (use.region=="second") { 
        ref <- expected2
    }
    
    olap <- findOverlaps(x, query.regions, type=type, maxgap=maxgap, minoverlap=minoverlap, use.region=use.region)
    expect_that(olap, is_identical_to(unique(sort(ref)))) 

    # Checking 'select' arguments.
    for (select in c("first", "last", "arbitrary")) { 
        selected <- findOverlaps(x, query.regions, type=type, maxgap=maxgap, minoverlap=minoverlap, select=select, use.region=use.region)
        if (select!="arbitrary") {
            expect_that(selected, is_identical_to(selectHits(olap, select)))
        } else {
            is.okay <- !is.na(selected)
            expect_identical(is.okay, !is.na(selectHits(olap, select)))
            arbiter <- Hits(which(is.okay), selected[is.okay], length(x), length(query.regions), sort.by.query=TRUE)
            expect_true(all(!is.na(match(arbiter, olap))))
        }
    }

    # Checking that 'overlapsAny', 'countOverlaps' and 'subsetByOverlaps' works.
    count.lap <- countOverlaps(x, query.regions, type=type, maxgap=maxgap, minoverlap=minoverlap, use.region=use.region)
    expect_identical(count.lap, selectHits(olap, "count"))
    out <- overlapsAny(x, query.regions, type=type, maxgap=maxgap, minoverlap=minoverlap, use.region=use.region)
    expect_identical(out, !is.na(selected))
    expect_equal(subsetByOverlaps(x, query.regions, type=type, maxgap=maxgap, minoverlap=minoverlap, use.region=use.region), x[out,])

    # Flipping it around:
    rolap <- findOverlaps(query.regions, x, type=type, maxgap=maxgap, minoverlap=minoverlap, use.region=use.region)
    rexpected1 <- findOverlaps(query.regions, anchors(x, type="first"), type=type, maxgap=maxgap, minoverlap=minoverlap)
    rexpected2 <- findOverlaps(query.regions, anchors(x, type="second"), type=type, maxgap=maxgap, minoverlap=minoverlap)
    if (use.region=="both") {
        rref <- Hits(c(queryHits(rexpected1), queryHits(rexpected2)), c(subjectHits(rexpected1), subjectHits(rexpected2)),
                     nLnode=Nq, nRnode=length(x), sort.by.query=TRUE)
    } else if (use.region=="first") {
        rref <- rexpected1
    } else if (use.region=="second") {
        rref <- rexpected2
    }
    
    expect_that(rolap, is_identical_to(unique(sort(rref)))) 
    if (type!="within") { expect_that(rolap, is_identical_to(t(olap))) } 

    for (select in c("first", "last", "arbitrary")) { 
        selected <- findOverlaps(query.regions, x, type=type, maxgap=maxgap, minoverlap=minoverlap, select=select, use.region=use.region)
        if (select!="arbitrary") {
            expect_that(selected, is_identical_to(selectHits(rolap, select)))
        } else {
            is.okay <- !is.na(selected)
            expect_identical(is.okay, !is.na(selectHits(rolap, select)))
            arbiter <- Hits(which(is.okay), selected[is.okay], length(query.regions), length(x), sort.by.query=TRUE)
            expect_true(all(!is.na(match(arbiter, rolap))))
        }
    }
    
    count.lap <- countOverlaps(query.regions, x, type=type, maxgap=maxgap, minoverlap=minoverlap, use.region=use.region)
    expect_identical(count.lap, selectHits(rolap, "count"))
    out <- overlapsAny(query.regions, x, type=type, maxgap=maxgap, minoverlap=minoverlap, use.region=use.region)
    expect_identical(out, !is.na(selected))
    expect_identical(subsetByOverlaps(query.regions, x, type=type, maxgap=maxgap, minoverlap=minoverlap, use.region=use.region), query.regions[out])
}
}

# What happens with silly inputs (GRanges)?
    
expect_equal(findOverlaps(x[0,], query.regions), Hits(nLnode=0, nRnode=Nq, sort.by.query=TRUE))
expect_equal(findOverlaps(x, query.regions[0]), Hits(nLnode=length(x), nRnode=0L, sort.by.query=TRUE))
expect_equal(findOverlaps(query.regions, x[0]), Hits(nRnode=0, nLnode=Nq, sort.by.query=TRUE))
expect_equal(findOverlaps(query.regions[0], x), Hits(nRnode=length(x), nLnode=0L, sort.by.query=TRUE))
expect_equal(findOverlaps(x[0,], query.regions[0]), Hits(nLnode=0L, nRnode=0L, sort.by.query=TRUE))
expect_equal(findOverlaps(query.regions[0,], x[0]), Hits(nLnode=0L, nRnode=0L, sort.by.query=TRUE))

expect_equal(findOverlaps(x[0,], query.regions, select="first"), integer(0))
expect_equal(findOverlaps(query.regions[0], x, select="first"), integer(0))
expect_equal(findOverlaps(x, query.regions[0], select="first"), rep(as.integer(NA), nrow(x)))
expect_equal(findOverlaps(query.regions, x[0], select="first"), rep(as.integer(NA), length(query.regions)))
expect_equal(overlapsAny(x[0,], query.regions), logical(0))
expect_equal(overlapsAny(query.regions[0]), logical(0))
expect_equal(overlapsAny(x, query.regions[0]), logical(nrow(x)))
expect_equal(overlapsAny(query.regions, x[0]), logical(length(query.regions)))
    
# What happens with silly inputs (GRangesList)?
   
empty.pairing <- GRangesList()
expect_equal(findOverlaps(x[0,], pairing), Hits(nLnode=0, nRnode=Nq2, sort.by.query=TRUE))
expect_equal(findOverlaps(x, empty.pairing), Hits(nLnode=length(x), nRnode=0L, sort.by.query=TRUE))
expect_equal(findOverlaps(pairing, x[0]), Hits(nRnode=0, nLnode=Nq2, sort.by.query=TRUE))
expect_equal(findOverlaps(empty.pairing, x), Hits(nRnode=length(x), nLnode=0L, sort.by.query=TRUE))
expect_equal(findOverlaps(x[0,], empty.pairing), Hits(nLnode=0L, nRnode=0L, sort.by.query=TRUE))
expect_equal(findOverlaps(empty.pairing, x[0]), Hits(nLnode=0L, nRnode=0L, sort.by.query=TRUE))

expect_equal(findOverlaps(x[0,], pairing, select="first"), integer(0))
expect_equal(findOverlaps(empty.pairing, x, select="first"), integer(0))
expect_equal(findOverlaps(x, empty.pairing, select="first"), rep(as.integer(NA), nrow(x)))
expect_equal(findOverlaps(pairing, x[0], select="first"), rep(as.integer(NA), Nq2))
expect_equal(overlapsAny(x[0,], pairing), logical(0))
expect_equal(overlapsAny(empty.pairing, x), logical(0))
expect_equal(overlapsAny(x, empty.pairing), logical(nrow(x)))
expect_equal(overlapsAny(pairing, x[0]), logical(Nq2))

#######################################################
# Paired overlaps with another InteractionSet (or GInteractions).

set.seed(303)

for (object in c("InteractionSet", "GInteractions")) {

next.starts <- round(runif(N, 1, 100))
next.ends <- next.starts + round(runif(N, 5, 20))
next.regions <- GRanges(rep(c("chrA", "chrB"), c(N-20, 20)), IRanges(next.starts, next.ends))
next.regions <- unique(next.regions)
N2 <- length(next.regions)

next.anchor1 <- sample(N2, Np, replace=TRUE)
next.anchor2 <- sample(N2, Np, replace=TRUE)
counts <- matrix(rpois(Np*Nlibs, lambda=10), ncol=Nlibs)
colnames(counts) <- seq_len(Nlibs)

if (object=="InteractionSet") {
x2 <- InteractionSet(counts, GInteractions(next.anchor1, next.anchor2, next.regions))
} else {
x2 <- GInteractions(next.anchor1, next.anchor2, next.regions)
}
pairing <- anchors(x2)

for (param in seq_len(6)) {
    type <- "any"
    maxgap <- 0L
    minoverlap <- 1L
    use.region <- "both"
    if (param==2L) { 
        maxgap <- 10L
    } else if (param==3L) {
        minoverlap <- 10L
    } else if (param==4L) {
        type <- "within"
    } else if (param==5L) {
        use.region <- "same"
    } else if(param==6L) {
        use.region <- "reverse"
    }

    # Overlapping with the InteractionSet:
    expected1.A <- findOverlaps(anchors(x, type="first"), pairing[[1]], type=type, maxgap=maxgap, minoverlap=minoverlap)
    expected1.B <- findOverlaps(anchors(x, type="first"), pairing[[2]], type=type, maxgap=maxgap, minoverlap=minoverlap)
    expected2.B <- findOverlaps(anchors(x, type="second"), pairing[[1]], type=type, maxgap=maxgap, minoverlap=minoverlap) # Yes, the A/B switch is deliberate.
    expected2.A <- findOverlaps(anchors(x, type="second"), pairing[[2]], type=type, maxgap=maxgap, minoverlap=minoverlap)

    expected1.A <- paste0(queryHits(expected1.A), ".", subjectHits(expected1.A), ".A")
    expected1.B <- paste0(queryHits(expected1.B), ".", subjectHits(expected1.B), ".B")
    expected2.A <- paste0(queryHits(expected2.A), ".", subjectHits(expected2.A), ".A")
    expected2.B <- paste0(queryHits(expected2.B), ".", subjectHits(expected2.B), ".B")
    if (use.region=="both") {
        expected1 <- c(expected1.A, expected1.B)
        expected2 <- c(expected2.A, expected2.B)
    } else if (use.region=="same") {
        expected1 <- expected1.A
        expected2 <- expected2.A
    } else {
        expected1 <- expected1.B
        expected2 <- expected2.B
    }

    expected <- intersect(expected1, expected2)
    harvest <- do.call(rbind, strsplit(expected, "\\."))
    ref <- Hits(as.integer(harvest[,1]), as.integer(harvest[,2]), nLnode=length(x), nRnode=length(x2), sort.by.query=TRUE)
    ref <- sort(unique(ref))
    
    olap <- findOverlaps(x, x2, type=type, maxgap=maxgap, minoverlap=minoverlap, use.region=use.region)
    expect_that(olap, is_identical_to(ref))

    # Checking 'select' arguments.
    for (select in c("first", "last", "arbitrary")) { 
        selected <- findOverlaps(x, x2, type=type, maxgap=maxgap, minoverlap=minoverlap, use.region=use.region, select=select)
        if (select!="arbitrary") {
            expect_that(selected, is_identical_to(selectHits(olap, select)))
        } else {
            is.okay <- !is.na(selected)
            expect_identical(is.okay, !is.na(selectHits(olap, select)))
            arbiter <- Hits(which(is.okay), selected[is.okay], length(x), length(x2), sort.by.query=TRUE)
            expect_true(all(!is.na(match(arbiter, olap))))
        }
    }

    # Checking that 'overlapsAny', 'countOverlaps' and 'subsetByOverlaps' works.
    count.lap <- countOverlaps(x, x2, type=type, maxgap=maxgap, minoverlap=minoverlap, use.region=use.region)
    expect_identical(count.lap, selectHits(olap, "count"))
    out <- overlapsAny(x, x2, type=type, maxgap=maxgap, minoverlap=minoverlap, use.region=use.region)
    expect_identical(out, !is.na(selected))
    expect_equal(subsetByOverlaps(x, x2, type=type, maxgap=maxgap, minoverlap=minoverlap, use.region=use.region), x[out,])

    # Note: No need to flip, it's the same method.
    # However, we do test self-overlaps as well.

    expected1.A <- findOverlaps(anchors(x, type="first"), anchors(x, type="first"), type=type, maxgap=maxgap, minoverlap=minoverlap)
    expected1.B <- findOverlaps(anchors(x, type="first"), anchors(x, type="second"), type=type, maxgap=maxgap, minoverlap=minoverlap)
    expected2.B <- findOverlaps(anchors(x, type="second"), anchors(x, type="first"), type=type, maxgap=maxgap, minoverlap=minoverlap) # Yes, the A/B switch is deliberate.
    expected2.A <- findOverlaps(anchors(x, type="second"), anchors(x, type="second"), type=type, maxgap=maxgap, minoverlap=minoverlap)

    expected1.A <- paste0(queryHits(expected1.A), ".", subjectHits(expected1.A), ".A")
    expected1.B <- paste0(queryHits(expected1.B), ".", subjectHits(expected1.B), ".B")
    expected2.A <- paste0(queryHits(expected2.A), ".", subjectHits(expected2.A), ".A")
    expected2.B <- paste0(queryHits(expected2.B), ".", subjectHits(expected2.B), ".B")
    if (use.region=="both") {
        expected1 <- c(expected1.A, expected1.B)
        expected2 <- c(expected2.A, expected2.B)
    } else if (use.region=="same") {
        expected1 <- expected1.A
        expected2 <- expected2.A
    } else {
        expected1 <- expected1.B
        expected2 <- expected2.B
    }

    expected <- intersect(expected1, expected2)
    harvest <- do.call(rbind, strsplit(expected, "\\."))
    ref <- Hits(as.integer(harvest[,1]), as.integer(harvest[,2]), nLnode=length(x), nRnode=length(x2), sort.by.query=TRUE)
    ref <- sort(unique(ref))
    
    self.olap <- findOverlaps(x, type=type, maxgap=maxgap, minoverlap=minoverlap, use.region=use.region)
    expect_that(self.olap, is_identical_to(ref))

    for (select in c("first", "last", "arbitrary")) { 
        selected <- findOverlaps(x, type=type, maxgap=maxgap, minoverlap=minoverlap, use.region=use.region, select=select)
        if (select!="arbitrary") {
            expect_that(selected, is_identical_to(selectHits(self.olap, select)))
        } else {
            is.okay <- !is.na(selected)
            expect_identical(is.okay, !is.na(selectHits(self.olap, select)))
            arbiter <- Hits(which(is.okay), selected[is.okay], length(x), length(x), sort.by.query=TRUE)
            expect_true(all(!is.na(match(arbiter, self.olap))))
        }
    }

    count.lap <- countOverlaps(x, type=type, maxgap=maxgap, minoverlap=minoverlap, use.region=use.region)
    expect_identical(count.lap, selectHits(self.olap, "count"))
    out <- overlapsAny(x, type=type, maxgap=maxgap, minoverlap=minoverlap, use.region=use.region)
    expect_identical(out, !is.na(selected))
}
}

# What happens with silly inputs?

expect_equal(findOverlaps(x[0], x2), Hits(nLnode=0, nRnode=Np, sort.by.query=TRUE))
expect_equal(findOverlaps(x, x2[0]), Hits(nRnode=0, nLnode=nrow(x), sort.by.query=TRUE))
expect_equal(findOverlaps(x[0], x2[0]), Hits(nLnode=0L, nRnode=0L, sort.by.query=TRUE))

expect_equal(findOverlaps(x[0], x2, select="first"), integer(0))
expect_equal(findOverlaps(x, x2[0], select="first"), rep(as.integer(NA), nrow(x)))
expect_equal(overlapsAny(x[0], x2), logical(0))
expect_equal(overlapsAny(x, x2[0]), logical(nrow(x)))

#######################################################
# overlapsAny for ContactMatrix objects

set.seed(304)
Nr <- 100
Nc <- 200
all.anchor1 <- sample(N, Nr, replace=TRUE)
all.anchor2 <- sample(N, Nc, replace=TRUE)
counts <- matrix(rpois(Nr*Nc, lambda=10), Nr, Nc)
x <- ContactMatrix(counts, all.anchor1, all.anchor2, all.regions)

Nq <- 6
query.starts <- round(runif(Nq, 1, 100))
query.ends <- query.starts + round(runif(Nq, 5, 20))
query.regions <- GRanges(rep(c("chrA", "chrB"), Nq/2), IRanges(query.starts, query.ends))

olap <- overlapsAny(x, query.regions)
expect_identical(olap, list(row=overlapsAny(anchors(x, type="row"), query.regions),
                            column=overlapsAny(anchors(x, type="column"), query.regions)))

olap <- overlapsAny(x, query.regions, type="within")
expect_identical(olap, list(row=overlapsAny(anchors(x, type="row"), query.regions, type="within"),
                            column=overlapsAny(anchors(x, type="column"), query.regions, type="within")))

expect_equal(overlapsAny(x[0,], query.regions), list(row=logical(0), column=overlapsAny(anchors(x, type="column"), query.regions)))
expect_equal(overlapsAny(x[,0], query.regions), list(row=overlapsAny(anchors(x, type="row"), query.regions), column=logical(0)))
expect_equal(overlapsAny(x[0,0], query.regions), list(row=logical(0), column=logical(0)))

# Trying out some 2D overlaps.

Nq3 <- 10
query.starts <- round(runif(Nq3, 1, 100))
query.ends <- query.starts + round(runif(Nq3, 5, 20))
query.regions1 <- GRanges(rep(c("chrA", "chrB"), Nq3/2), IRanges(query.starts, query.ends))
query.starts <- round(runif(Nq3, 1, 100))
query.ends <- query.starts + round(runif(Nq3, 5, 20))
query.regions2 <- GRanges(rep(c("chrA", "chrB"), Nq3/2), IRanges(query.starts, query.ends))
pairing <- GInteractions(query.regions1, query.regions2)

olap <- overlapsAny(x, pairing)
temp.iset <- deflate(x, collapse=TRUE)
ref <- overlapsAny(temp.iset, pairing) # no NAs, everyone's represented here.
ref <- inflate(temp.iset, anchors(x, type="row", id=TRUE), anchors(x, type="column", id=TRUE), fill=ref)
expect_identical(olap, unname(as.matrix(as.matrix(ref))))

olap <- overlapsAny(x, pairing, type="within")
ref <- overlapsAny(temp.iset, pairing, type="within")
ref <- inflate(temp.iset, anchors(x, type="row", id=TRUE), anchors(x, type="column", id=TRUE), fill=ref)
expect_identical(olap, unname(as.matrix(as.matrix(ref))))

olap <- overlapsAny(x, pairing)
new.gi <- GInteractions(query.regions1, query.regions2)
expect_identical(olap, overlapsAny(x, new.gi))

new.iset <- InteractionSet(matrix(runif(10), dimnames=list(NULL, 1)), new.gi)
expect_identical(olap, overlapsAny(x, new.iset))

# Checking with order enforcement.

olap <- overlapsAny(x, pairing, use.region="same")
expect_identical(olap, overlapsAny(x, new.gi, use.region="same"))

olap <- overlapsAny(x, pairing, use.region="reverse")
expect_identical(olap, overlapsAny(x, new.gi, use.region="reverse"))

#######################################################
# End
