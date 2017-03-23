for (siglist in c("GInteractions", "ContactMatrix")) {
    setMethod("trim", siglist, function(x, use.names=TRUE) {
        out <- trim(regions(x), use.names=use.names)
        collated <- .collate_GRanges(out)
        unchecked_regions(x) <- collated$ranges
        unchecked_anchor1(x) <- collated$indices[[1]][anchor1(x)]
        unchecked_anchor2(x) <- collated$indices[[1]][anchor2(x)]
        return(x)
    })
}

setMethod("width", "GInteractions", function(x) {
    w <- width(regions(x))          
    DataFrame(anchor1=w[anchor1(x)], anchor2=w[anchor2(x)])
})

setMethod("width", "ContactMatrix", function(x) {
    w <- width(regions(x))          
    list(anchor1=w[anchor1(x)], anchor2=w[anchor2(x)])
})

setMethod("trim", "InteractionSet", function(x, use.names=TRUE) { 
    interactions(x) <- trim(interactions(x), use.names=use.names) 
    return(x)
})
setMethod("width", "InteractionSet", function(x) { width(interactions(x)) })
