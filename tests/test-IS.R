suppressPackageStartupMessages(require(InteractionSet))

############################
# InteractionSet testing:

out <- InteractionSet(matrix(0, 4, 2), 1:4, 4:1, GRanges("chr1", IRanges(1:4, 1:4)), 
            colData=DataFrame(whee=1:2), metadata=list())

# Getter test
regions(out)
anchors(out)
anchors(out, type="first")
anchors(out, type="second")

anchors(out, id=TRUE)
anchors(out, type="first", id=TRUE)
anchors(out, type="second", id=TRUE)

# Setter test
regions(out)$x <- "whee"
regions(out)
anchors(out) <- list(1:4, 1:4)
anchors(out, id=TRUE)
anchors(out) <- list(rep(4L, 4), 1:4)
anchors(out, id=TRUE)

# Subset test
out[,2]
out[1,]
anchors(out[4:1,])

# Combining.
cbind(out, out)
rbind(out, out)

elementMetadata(regions(out)) <- NULL
out2 <- InteractionSet(matrix(0, 4, 2), 1:4, 4:1, GRanges("chr1", IRanges(2:5, 2:5)))
try(rbind(out, out2))
out3 <- c(out, out2)
regions(out3)
anchors(out3, id=TRUE)

# Other methods test
order(out)
order(out[4:1,])
sort(out)

sout <- split(out, c(1,1,1,2))
sout[[1]]
anchors(sout[[1]])

# Alternative construction
out <- InteractionSet(matrix(0, 2, 2), GRanges("chr1", IRanges(4:3, 4:3)),
           GRanges("chr1", IRanges(2:1, 2:1))) 
anchors(out)
out <- InteractionSet(matrix(0, 2, 2), GRanges("chr1", IRanges(4:3, 4:3)),
           GRanges("chr1", IRanges(2:1, 2:1)), regions=GRanges("chr1", IRanges(1:5, 1:5))) 
anchors(out, id=TRUE)
                                              
############################
# Testing findOverlaps:

out <- InteractionSet(matrix(0, 4, 2), 1:4, 4:1, GRanges("chr1", IRanges(1:4, 1:4)), 
            colData=DataFrame(whee=1:2), metadata=list())

findOverlaps(out, GRanges("chr1", IRanges(1, 1)))
overlapsAny(out, GRanges("chr1", IRanges(1, 1)))
findOverlaps(out, GRanges("chr1", IRanges(1, 2)))
findOverlaps(out, GRanges("chr1", IRanges(1:2, 1:2)))

subsetByOverlaps(out, GRanges("chr1", IRanges(2, 2)))
subsetByOverlaps(GRanges("chr1", IRanges(2, 2)), out)

findOverlaps(out, GRangesList(first=GRanges("chr1", IRanges(1:2, 1:2)),
                              second=GRanges("chr1", IRanges(1:2, 1:2))))
overlapsAny(out, GRangesList(first=GRanges("chr1", IRanges(1:2, 1:2)),
                              second=GRanges("chr1", IRanges(1:2, 1:2))))
findOverlaps(out, GRangesList(first=GRanges("chr1", IRanges(2, 2)),
                              second=GRanges("chr1", IRanges(3, 3))))
overlapsAny(out, GRangesList(first=GRanges("chr1", IRanges(2, 2)),
                   second=GRanges("chr1", IRanges(3, 3))))
overlapsAny(GRangesList(first=GRanges("chr1", IRanges(1:4, 1:4)),
        second=GRanges("chr1", IRanges(1:4, 1:4))), out)

out2 <- InteractionSet(matrix(0, 4, 2), 1:4, 1:4, GRanges("chr1", IRanges(1:4, 1:4))) 
findOverlaps(out, out2)

out2 <- InteractionSet(matrix(0, 4, 2), c(1,2,2,4), 1:4, GRanges("chr1", IRanges(1:4, 1:4))) 
findOverlaps(out, out2)
overlapsAny(out, out2)

subsetByOverlaps(out, out2)

##############################
# Testing distance functions:

out <- InteractionSet(matrix(0, 4, 2), 1:4, 4:1, GRanges("chr1", IRanges(1:4, 1:4)), 
            colData=DataFrame(whee=1:2), metadata=list())
pairdist(out)
pairdist(out, type="gap")
pairdist(out, type="span")
pairdist(out, type="diag")
intrachr(out)

