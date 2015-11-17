# Testing the various overlap methods for InteractionSet objects.

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
x <- InteractionSet(counts, GInteractions(all.anchor1, all.anchor2, all.regions))

#######################################################
# Linear overlaps with GRanges.

set.seed(301)
Nq <- 10
query.starts <- round(runif(Nq, 1, 100))
query.ends <- query.starts + round(runif(Nq, 5, 20))
query.regions <- GRanges(rep(c("chrA", "chrB"), Nq/2), IRanges(query.starts, query.ends))

for (param in seq_len(4)) {
    type <- "any"
    maxgap <- 0L
    minoverlap <- 1L
    if (param==2L) { 
        maxgap <- 10L
    } else if (param==3L) {
        minoverlap <- 10L
    } else if (param==4L) {
        type <- "within"
    }

    # Overlapping with GRanges.
    expected1 <- findOverlaps(anchors(x, type="first"), query.regions, type=type, maxgap=maxgap, minoverlap=minoverlap)
    expected2 <- findOverlaps(anchors(x, type="second"), query.regions, type=type, maxgap=maxgap, minoverlap=minoverlap)
    ref <- Hits(c(queryHits(expected1), queryHits(expected2)), c(subjectHits(expected1), subjectHits(expected2)),
                queryLength=length(x), subjectLength=Nq)
    
    olap <- findOverlaps(x, query.regions, type=type, maxgap=maxgap, minoverlap=minoverlap)
    expect_that(olap, is_identical_to(unique(sort(ref)))) 

    # Checking 'select' arguments.
    for (select in c("first", "last", "arbitrary")) { 
        selected <- findOverlaps(x, query.regions, type=type, maxgap=maxgap, minoverlap=minoverlap, select=select)
        expect_that(selected, is_identical_to(selectHits(olap, select)))
    }

    # Checking that 'overlapsAny', 'countOverlaps' and 'subsetByOverlaps' works.
    count.lap <- countOverlaps(x, query.regions, type=type, maxgap=maxgap, minoverlap=minoverlap)
    expect_identical(count.lap, selectHits(olap, "count"))
    out <- overlapsAny(x, query.regions, type=type, maxgap=maxgap, minoverlap=minoverlap)
    expect_identical(out, !is.na(selected))
    expect_equal(subsetByOverlaps(x, query.regions, type=type, maxgap=maxgap, minoverlap=minoverlap), x[out,])

    # Flipping it around:
    rolap <- findOverlaps(query.regions, x, type=type, maxgap=maxgap, minoverlap=minoverlap)
    rexpected1 <- findOverlaps(query.regions, anchors(x, type="first"), type=type, maxgap=maxgap, minoverlap=minoverlap)
    rexpected2 <- findOverlaps(query.regions, anchors(x, type="second"), type=type, maxgap=maxgap, minoverlap=minoverlap)
    rref <- Hits(c(queryHits(rexpected1), queryHits(rexpected2)), c(subjectHits(rexpected1), subjectHits(rexpected2)),
                queryLength=Nq, subjectLength=length(x))
    
    expect_that(rolap, is_identical_to(unique(sort(rref)))) 
    if (type!="within") { expect_that(rolap, is_identical_to(t(olap))) } 

    for (select in c("first", "last", "arbitrary")) { 
        selected <- findOverlaps(query.regions, x, type=type, maxgap=maxgap, minoverlap=minoverlap, select=select)
        expect_that(selected, is_identical_to(selectHits(rolap, select)))
    }
    
    count.lap <- countOverlaps(query.regions, x, type=type, maxgap=maxgap, minoverlap=minoverlap)
    expect_identical(count.lap, selectHits(rolap, "count"))
    out <- overlapsAny(query.regions, x, type=type, maxgap=maxgap, minoverlap=minoverlap)
    expect_identical(out, !is.na(selected))
    expect_identical(subsetByOverlaps(query.regions, x, type=type, maxgap=maxgap, minoverlap=minoverlap), query.regions[out])
}

# What happens with silly inputs?
    
expect_equal(findOverlaps(x[0,], query.regions), Hits(queryLength=0, subjectLength=Nq))
expect_equal(findOverlaps(x, query.regions[0]), Hits(queryLength=length(x), subjectLength=0L))
expect_equal(findOverlaps(query.regions, x[0]), Hits(subjectLength=0, queryLength=Nq))
expect_equal(findOverlaps(query.regions[0], x), Hits(subjectLength=length(x), queryLength=0L))
expect_equal(findOverlaps(x[0,], query.regions[0]), Hits(queryLength=0L, subjectLength=0L))
expect_equal(findOverlaps(query.regions[0,], x[0]), Hits(queryLength=0L, subjectLength=0L))

expect_equal(findOverlaps(x[0,], query.regions, select="first"), integer(0))
expect_equal(findOverlaps(query.regions[0], x, select="first"), integer(0))
expect_equal(findOverlaps(x, query.regions[0], select="first"), rep(as.integer(NA), nrow(x)))
expect_equal(findOverlaps(query.regions, x[0], select="first"), rep(as.integer(NA), length(query.regions)))
expect_equal(overlapsAny(x[0,], query.regions), logical(0))
expect_equal(overlapsAny(query.regions[0]), logical(0))
expect_equal(overlapsAny(x, query.regions[0]), logical(nrow(x)))
expect_equal(overlapsAny(query.regions, x[0]), logical(length(query.regions)))

#######################################################
# Paired overlaps with a GRangesList.

set.seed(302)
Nq2 <- 100
query.starts <- round(runif(Nq2, 1, 100))
query.ends <- query.starts + round(runif(Nq2, 5, 20))
query.regions1 <- GRanges(rep(c("chrA", "chrB"), Nq2/2), IRanges(query.starts, query.ends))
query.starts <- round(runif(Nq2, 1, 100))
query.ends <- query.starts + round(runif(Nq2, 5, 20))
query.regions2 <- GRanges(rep(c("chrA", "chrB"), Nq2/2), IRanges(query.starts, query.ends))
pairing <- GRangesList(first=query.regions1, second=query.regions2)

for (param in seq_len(4)) {
    type <- "any"
    maxgap <- 0L
    minoverlap <- 1L
    if (param==2L) { 
        maxgap <- 10L
    } else if (param==3L) {
        minoverlap <- 10L
    } else if (param==4L) {
        type <- "within"
    }

    # Overlapping with the GRangesList:
    expected1.A <- findOverlaps(anchors(x, type="first"), pairing[[1]], type=type, maxgap=maxgap, minoverlap=minoverlap)
    expected1.B <- findOverlaps(anchors(x, type="first"), pairing[[2]], type=type, maxgap=maxgap, minoverlap=minoverlap)
    expected2.B <- findOverlaps(anchors(x, type="second"), pairing[[1]], type=type, maxgap=maxgap, minoverlap=minoverlap) # Yes, the A/B switch is deliberate.
    expected2.A <- findOverlaps(anchors(x, type="second"), pairing[[2]], type=type, maxgap=maxgap, minoverlap=minoverlap)

    expected1 <- c(paste0(queryHits(expected1.A), ".", subjectHits(expected1.A), ".A"), 
                   paste0(queryHits(expected1.B), ".", subjectHits(expected1.B), ".B"))
    expected2 <- c(paste0(queryHits(expected2.A), ".", subjectHits(expected2.A), ".A"), 
                   paste0(queryHits(expected2.B), ".", subjectHits(expected2.B), ".B"))
    expected <- intersect(expected1, expected2)
    harvest <- do.call(rbind, strsplit(expected, "\\."))
    ref <- Hits(as.integer(harvest[,1]), as.integer(harvest[,2]), queryLength=length(x), subjectLength=Nq2)
    ref <- sort(unique(ref))
    
    olap <- findOverlaps(x, pairing, type=type, maxgap=maxgap, minoverlap=minoverlap)
    expect_that(olap, is_identical_to(ref))

    # Checking 'select' arguments.
    for (select in c("first", "last", "arbitrary")) { 
        selected <- findOverlaps(x, pairing, type=type, maxgap=maxgap, minoverlap=minoverlap, select=select)
        expect_that(selected, is_identical_to(selectHits(olap, select)))
    }

    # Checking that 'overlapsAny', 'countOverlaps' and 'subsetByOverlaps' works.
    count.lap <- countOverlaps(x, pairing, type=type, maxgap=maxgap, minoverlap=minoverlap)
    expect_identical(count.lap, selectHits(olap, "count"))
    out <- overlapsAny(x, pairing, type=type, maxgap=maxgap, minoverlap=minoverlap)
    expect_identical(out, !is.na(selected))
    expect_equal(subsetByOverlaps(x, pairing, type=type, maxgap=maxgap, minoverlap=minoverlap), x[out,])

    # Flipping it around:
    rexpected1.A <- findOverlaps(pairing[[1]], anchors(x, type="first"), type=type, maxgap=maxgap, minoverlap=minoverlap)
    rexpected1.B <- findOverlaps(pairing[[2]], anchors(x, type="first"), type=type, maxgap=maxgap, minoverlap=minoverlap)
    rexpected2.B <- findOverlaps(pairing[[1]], anchors(x, type="second"), type=type, maxgap=maxgap, minoverlap=minoverlap) 
    rexpected2.A <- findOverlaps(pairing[[2]], anchors(x, type="second"), type=type, maxgap=maxgap, minoverlap=minoverlap)

    rexpected1 <- c(paste0(queryHits(rexpected1.A), ".", subjectHits(rexpected1.A), ".A"), 
                   paste0(queryHits(rexpected1.B), ".", subjectHits(rexpected1.B), ".B"))
    rexpected2 <- c(paste0(queryHits(rexpected2.A), ".", subjectHits(rexpected2.A), ".A"), 
                   paste0(queryHits(rexpected2.B), ".", subjectHits(rexpected2.B), ".B"))
    rexpected <- intersect(rexpected1, rexpected2)
    rharvest <- do.call(rbind, strsplit(rexpected, "\\."))
    rref <- Hits(as.integer(rharvest[,1]), as.integer(rharvest[,2]), queryLength=length(x), subjectLength=Nq2)
    rref <- sort(unique(rref))

    rolap <- findOverlaps(pairing, x, type=type, maxgap=maxgap, minoverlap=minoverlap)
    expect_that(rolap, is_identical_to(unique(sort(rref)))) 
    if (type!="within") { expect_that(rolap, is_identical_to(t(olap))) } 
    
    for (select in c("first", "last", "arbitrary")) { 
        selected <- findOverlaps(pairing, x, type=type, maxgap=maxgap, minoverlap=minoverlap, select=select)
        expect_that(selected, is_identical_to(selectHits(rolap, select)))
    }

    count.lap <- countOverlaps(pairing, x, type=type, maxgap=maxgap, minoverlap=minoverlap)
    expect_identical(count.lap, selectHits(rolap, "count"))
    out <- overlapsAny(pairing, x, type=type, maxgap=maxgap, minoverlap=minoverlap)
    expect_identical(out, !is.na(selected))
    expect_identical(subsetByOverlaps(pairing, x, type=type, maxgap=maxgap, minoverlap=minoverlap), 
           GRangesList(first=query.regions1[out], second=query.regions2[out]))
}
    
# What happens with silly inputs?
   
expect_error(findOverlaps(x, GRangesList(GRanges(), query.regions1)), "component GRanges in the GRangesList must be of the same length", fixed=TRUE)
expect_error(findOverlaps(x, GRangesList(query.regions1, query.regions2, query.regions1)), "input GRangesList must be of length 2", fixed=TRUE)

empty.pairing <- GRangesList(GRanges(), GRanges())
expect_equal(findOverlaps(x[0,], pairing), Hits(queryLength=0, subjectLength=Nq2))
expect_equal(findOverlaps(x, empty.pairing), Hits(queryLength=length(x), subjectLength=0L))
expect_equal(findOverlaps(pairing, x[0]), Hits(subjectLength=0, queryLength=Nq2))
expect_equal(findOverlaps(empty.pairing, x), Hits(subjectLength=length(x), queryLength=0L))
expect_equal(findOverlaps(x[0,], empty.pairing), Hits(queryLength=0L, subjectLength=0L))
expect_equal(findOverlaps(empty.pairing, x[0]), Hits(queryLength=0L, subjectLength=0L))

expect_equal(findOverlaps(x[0,], pairing, select="first"), integer(0))
expect_equal(findOverlaps(empty.pairing, x, select="first"), integer(0))
expect_equal(findOverlaps(x, empty.pairing, select="first"), rep(as.integer(NA), nrow(x)))
expect_equal(findOverlaps(pairing, x[0], select="first"), rep(as.integer(NA), Nq2))
expect_equal(overlapsAny(x[0,], pairing), logical(0))
expect_equal(overlapsAny(empty.pairing, x), logical(0))
expect_equal(overlapsAny(x, empty.pairing), logical(nrow(x)))
expect_equal(overlapsAny(pairing, x[0]), logical(Nq2))

#######################################################
# Paired overlaps with another InteractionSet.

set.seed(303)
next.starts <- round(runif(N, 1, 100))
next.ends <- next.starts + round(runif(N, 5, 20))
next.regions <- GRanges(rep(c("chrA", "chrB"), c(N-20, 20)), IRanges(next.starts, next.ends))
next.regions <- unique(next.regions)
N2 <- length(next.regions)

next.anchor1 <- sample(N2, Np, replace=TRUE)
next.anchor2 <- sample(N2, Np, replace=TRUE)
counts <- matrix(rpois(Np*Nlibs, lambda=10), ncol=Nlibs)
x2 <- InteractionSet(counts, GInteractions(next.anchor1, next.anchor2, next.regions))
pairing <- anchors(x2)

for (param in seq_len(4)) {
    type <- "any"
    maxgap <- 0L
    minoverlap <- 1L
    if (param==2L) { 
        maxgap <- 10L
    } else if (param==3L) {
        minoverlap <- 10L
    } else if (param==4L) {
        type <- "within"
    }

    # Overlapping with the GRangesList:
    expected1.A <- findOverlaps(anchors(x, type="first"), pairing[[1]], type=type, maxgap=maxgap, minoverlap=minoverlap)
    expected1.B <- findOverlaps(anchors(x, type="first"), pairing[[2]], type=type, maxgap=maxgap, minoverlap=minoverlap)
    expected2.B <- findOverlaps(anchors(x, type="second"), pairing[[1]], type=type, maxgap=maxgap, minoverlap=minoverlap) # Yes, the A/B switch is deliberate.
    expected2.A <- findOverlaps(anchors(x, type="second"), pairing[[2]], type=type, maxgap=maxgap, minoverlap=minoverlap)

    expected1 <- c(paste0(queryHits(expected1.A), ".", subjectHits(expected1.A), ".A"), 
                   paste0(queryHits(expected1.B), ".", subjectHits(expected1.B), ".B"))
    expected2 <- c(paste0(queryHits(expected2.A), ".", subjectHits(expected2.A), ".A"), 
                   paste0(queryHits(expected2.B), ".", subjectHits(expected2.B), ".B"))
    expected <- intersect(expected1, expected2)
    harvest <- do.call(rbind, strsplit(expected, "\\."))
    ref <- Hits(as.integer(harvest[,1]), as.integer(harvest[,2]), queryLength=length(x), subjectLength=length(x2))
    ref <- sort(unique(ref))
    
    olap <- findOverlaps(x, x2, type=type, maxgap=maxgap, minoverlap=minoverlap)
    expect_that(olap, is_identical_to(ref))

    # Checking 'select' arguments.
    for (select in c("first", "last", "arbitrary")) { 
        selected <- findOverlaps(x, x2, type=type, maxgap=maxgap, minoverlap=minoverlap, select=select)
        expect_that(selected, is_identical_to(selectHits(olap, select)))
    }

    # Checking that 'overlapsAny', 'countOverlaps' and 'subsetByOverlaps' works.
    count.lap <- countOverlaps(x, x2, type=type, maxgap=maxgap, minoverlap=minoverlap)
    expect_identical(count.lap, selectHits(olap, "count"))
    out <- overlapsAny(x, x2, type=type, maxgap=maxgap, minoverlap=minoverlap)
    expect_identical(out, !is.na(selected))
    expect_equal(subsetByOverlaps(x, x2, type=type, maxgap=maxgap, minoverlap=minoverlap), x[out,])

    # Note: No need to flip, it's the same method.
}

# What happens with silly inputs?

expect_equal(findOverlaps(x[0], x2), Hits(queryLength=0, subjectLength=nrow(x2)))
expect_equal(findOverlaps(x, x2[0]), Hits(subjectLength=0, queryLength=nrow(x)))
expect_equal(findOverlaps(x[0], x2[0]), Hits(queryLength=0L, subjectLength=0L))

expect_equal(findOverlaps(x[0], x2, select="first"), integer(0))
expect_equal(findOverlaps(x, x2[0], select="first"), rep(as.integer(NA), nrow(x)))
expect_equal(overlapsAny(x[0], x2), logical(0))
expect_equal(overlapsAny(x, x2[0]), logical(Nq2))

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
pairing <- GRangesList(first=query.regions1, second=query.regions2)

olap <- overlapsAny(x, pairing)
temp.iset <- deflate(x, unique=TRUE)
ref <- overlapsAny(temp.iset, pairing) # no NAs, everyone's represented here.
ref <- inflate(temp.iset, anchors(x, type="row", id=TRUE), anchors(x, type="column", id=TRUE), fill=ref)
expect_identical(olap, as.matrix(ref))

olap <- overlapsAny(x, pairing, type="within")
ref <- overlapsAny(temp.iset, pairing, type="within")
ref <- inflate(temp.iset, anchors(x, type="row", id=TRUE), anchors(x, type="column", id=TRUE), fill=ref)
expect_identical(olap, as.matrix(ref))

olap <- overlapsAny(x, pairing)
new.gi <- GInteractions(query.regions1, query.regions2)
expect_identical(olap, overlapsAny(x, new.gi))

new.iset <- InteractionSet(as.matrix(runif(10)), new.gi)
expect_identical(olap, overlapsAny(x, new.iset))

#######################################################
# End
