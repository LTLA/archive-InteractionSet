# Tests the linkOverlaps method

set.seed(9000)
N <- 30
all.starts <- round(runif(N, 1, 100))
all.ends <- all.starts + round(runif(N, 5, 20))
all.regions <- GRanges(rep(c("chrA", "chrB"), c(N-10, 10)), IRanges(all.starts, all.ends))

Np <- 20
all.anchor1 <- sample(N, Np)
all.anchor2 <- sample(N, Np)
x <- GInteractions(all.anchor1, all.anchor2, all.regions)

# Generating some random regions to test against.

Ngenes <- 10
gene.starts <- round(runif(Ngenes, 1, 100))
gene.ends <- gene.starts + round(runif(Ngenes, 5, 20))
gene.regions <- GRanges(rep(c("chrA", "chrB"), Ngenes/2), IRanges(gene.starts, gene.ends))

Nenh <- 4
enh.starts <- round(runif(Nenh, 1, 100))
enh.ends <- enh.starts + round(runif(Nenh, 5, 20))
enh.regions <- GRanges(rep(c("chrA", "chrB"), Nenh/2), IRanges(enh.starts, enh.ends))

# Testing for ISets and GIs, in a range of scenarios.

for (cls in 1:2) { 
    if (cls==1L) {
        obj <- x
    } else {
        obj <- InteractionSet(matrix(0, Np, 4), x)
    }
    
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

        # Getting the reference results by doing it manually via 'merge'.
        olap1.g <- findOverlaps(anchors(obj, type="first"), gene.regions, maxgap=maxgap, minoverlap=minoverlap, type=type)
        olap2.g <- findOverlaps(anchors(obj, type="second"), gene.regions, maxgap=maxgap, minoverlap=minoverlap, type=type)
        olap1.e <- findOverlaps(anchors(obj, type="first"), enh.regions, maxgap=maxgap, minoverlap=minoverlap, type=type)
        olap2.e <- findOverlaps(anchors(obj, type="second"), enh.regions, maxgap=maxgap, minoverlap=minoverlap, type=type)

        combo <- rbind(merge(olap1.g, olap2.e, by.x=1, by.y=1), merge(olap2.g, olap1.e, by.x=1, by.y=1))
        colnames(combo) <- c("query", "subject1", "subject2")
        is.dup <- duplicated(paste0(combo$query, ".", combo$subject1, ".", combo$subject2))
        combo <- combo[!is.dup,]
        o <- order(combo$query, combo$subject1, combo$subject2)
        combo <- combo[o,]
        rownames(combo) <- NULL

        expect_identical(combo, linkOverlaps(obj, gene.regions, enh.regions, type=type, maxgap=maxgap, minoverlap=minoverlap))

        # Testing against self-interactions.
        combo2 <- merge(olap1.g, olap2.g, by.x=1, by.y=1)
        colnames(combo2) <- c("query", "subject1", "subject2")
        new.s1 <- pmax(combo2$subject1, combo2$subject2)
        new.s2 <- pmin(combo2$subject1, combo2$subject2)
        combo2$subject1 <- new.s1
        combo2$subject2 <- new.s2

        is.dup <- duplicated(paste0(combo2$query, ".", combo2$subject1, ".", combo2$subject2)) | is.na(combo2$query)
        combo2 <- combo2[!is.dup,]
        o <- order(combo2$query, combo2$subject1, combo2$subject2)
        combo2 <- combo2[o,]
        rownames(combo2) <- NULL

        expect_identical(combo2, linkOverlaps(obj, gene.regions, type=type, maxgap=maxgap, minoverlap=minoverlap))
    }

    # Testing against empty slots.
    expect_identical(data.frame(query=integer(0), subject1=integer(0), subject2=integer(0)), linkOverlaps(obj, gene.regions[0], gene.regions))
    expect_identical(data.frame(query=integer(0), subject1=integer(0), subject2=integer(0)), linkOverlaps(obj, gene.regions, gene.regions[0]))
    expect_identical(data.frame(query=integer(0), subject1=integer(0), subject2=integer(0)), linkOverlaps(obj, gene.regions[0], gene.regions[0]))
}


