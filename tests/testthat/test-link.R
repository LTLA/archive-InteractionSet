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
        obj <- InteractionSet(matrix(0, Np, 4, dimnames=list(NULL, seq_len(4))), x)
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
            use.region <- "same"
        } else if (param==6L) {
            use.region <- "reverse"
        }

        # Getting the reference results by doing it manually via 'merge'.
        olap1.g <- findOverlaps(anchors(obj, type="first"), gene.regions, maxgap=maxgap, minoverlap=minoverlap, type=type)
        olap2.g <- findOverlaps(anchors(obj, type="second"), gene.regions, maxgap=maxgap, minoverlap=minoverlap, type=type)
        olap1.e <- findOverlaps(anchors(obj, type="first"), enh.regions, maxgap=maxgap, minoverlap=minoverlap, type=type)
        olap2.e <- findOverlaps(anchors(obj, type="second"), enh.regions, maxgap=maxgap, minoverlap=minoverlap, type=type)

        combo1 <- base::merge(olap1.g, olap2.e, by.x=1, by.y=1)
        combo2 <- base::merge(olap2.g, olap1.e, by.x=1, by.y=1)
        if (use.region=="both") { 
            combo <- rbind(combo1, combo2)
        } else if (use.region=="same") {
            combo <- combo1
        } else {
            combo <- combo2
        }
        colnames(combo) <- c("query", "subject1", "subject2")
        is.dup <- duplicated(paste0(combo$query, ".", combo$subject1, ".", combo$subject2))
        combo <- combo[!is.dup,]
        o <- order(combo$query, combo$subject1, combo$subject2)
        combo <- combo[o,]
        rownames(combo) <- NULL

        expect_identical(combo, linkOverlaps(obj, gene.regions, enh.regions, type=type, maxgap=maxgap, minoverlap=minoverlap, use.region=use.region))

        # Testing against self-interactions.
        combo.S <- base::merge(olap1.g, olap2.g, by.x=1, by.y=1)
        colnames(combo.S) <- c("query", "subject1", "subject2")
        new.s1 <- pmax(combo.S$subject1, combo.S$subject2)
        new.s2 <- pmin(combo.S$subject1, combo.S$subject2)
        combo.S$subject1 <- new.s1
        combo.S$subject2 <- new.s2

        is.dup <- duplicated(paste0(combo.S$query, ".", combo.S$subject1, ".", combo.S$subject2)) | is.na(combo.S$query)
        combo.S <- combo.S[!is.dup,]
        o <- order(combo.S$query, combo.S$subject1, combo.S$subject2)
        combo.S <- combo.S[o,]
        rownames(combo.S) <- NULL

        expect_identical(combo.S, linkOverlaps(obj, gene.regions, type=type, maxgap=maxgap, minoverlap=minoverlap, use.region=use.region))
    }

    # Testing against empty slots.
    expect_identical(data.frame(query=integer(0), subject1=integer(0), subject2=integer(0)), linkOverlaps(obj, gene.regions[0], gene.regions))
    expect_identical(data.frame(query=integer(0), subject1=integer(0), subject2=integer(0)), linkOverlaps(obj, gene.regions, gene.regions[0]))
    expect_identical(data.frame(query=integer(0), subject1=integer(0), subject2=integer(0)), linkOverlaps(obj, gene.regions[0], gene.regions[0]))
}


