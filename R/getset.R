###############################################################
# Getters:

setGeneric("anchors", function(x, ...) { standardGeneric("anchors") })
setGeneric("regions", function(x, ...) { standardGeneric("regions") })

# A generating function, to capture differences in 'type' for 'anchors' call.
GI.args <- c("both", "first", "second") 
CM.args <- c("both", "row", "column") 
anchor.fun.gen <- function(is.GI) { 
    if (is.GI) { 
        type.arg <- GI.args
    } else {
        type.arg <- CM.args
    }
    out.names <- type.arg[2:3]
    type1 <- out.names[1]

    function(x, type="both", id=FALSE) {
        type <- match.arg(type, type.arg)
        if (id) {
            if (type=="both") {
                out <- list(x@anchor1, x@anchor2)
                names(out) <- out.names
                return(out)
            } else if (type==type1) {
                return(x@anchor1)
            } else {
                return(x@anchor2)
            }
        } else {
            if (type=="both") {
                out <- GRangesList(x@regions[x@anchor1], x@regions[x@anchor2])
                names(out) <- out.names
                return(out)
            } else if (type==type1) {
                return(x@regions[x@anchor1])
            } else {
                return(x@regions[x@anchor2])
            }
        }
    }
}

# Defining the methods:
setMethod("anchors", "GInteractions", anchor.fun.gen(TRUE))
setMethod("anchors", "ContactMatrix", anchor.fun.gen(FALSE))

for (siglist in list("GInteractions", "ContactMatrix")) {
    setMethod("regions", siglist, function(x) { x@regions })
}

###############################################################
# Setters:

setGeneric("regions<-", function(x, value) { standardGeneric("regions<-") })
for (siglist in c("GInteractions", "ContactMatrix")) { 
    setReplaceMethod("regions", siglist, function(x, value) {
        if (length(value)!=length(regions(x))) { 
            stop("assigned value must be of the same length as 'regions(x)'")
        }
        out <- .resort_regions(x@anchor1, x@anchor2, value)
        x@anchor1 <- out$anchor1
        x@anchor2 <- out$anchor2
        x@regions <- out$regions
        validObject(x)
        return(x)
    })
}

# Also allow setting of regions of different length.

setGeneric("replaceRegions<-", function(x, value) { standardGeneric("replaceRegions<-") })
for (siglist in c("GInteractions", "ContactMatrix")) { 
    setReplaceMethod("replaceRegions", siglist, function(x, value) {
        converter <- match(regions(x), value)
        new.a1 <- converter[x@anchor1]
        new.a2 <- converter[x@anchor2]
        if (any(is.na(new.a1)) || any(is.na(new.a2))) { 
            stop("some existing ranges do not exist in replacement GRanges") 
        }

        out <- .resort_regions(new.a1, new.a2, value)
        x@anchor1 <- out$anchor1
        x@anchor2 <- out$anchor2
        x@regions <- out$regions
        return(x)
    })
}

# Append regions.

setGeneric("appendRegions<-", function(x, value) { standardGeneric("appendRegions<-") })
for (siglist in c("GInteractions", "ContactMatrix")) {
    setReplaceMethod("appendRegions", siglist, function(x, value) {
        out <- .resort_regions(x@anchor1, x@anchor2, c(x@regions, value)) 
        x@anchor1 <- out$anchor1
        x@anchor2 <- out$anchor2
        x@regions <- out$regions
        return(x)
    })
}

# Reduce regions.

setGeneric("reduceRegions", function(x) { standardGeneric("reduceRegions") })
for (siglist in c("GInteractions", "ContactMatrix")) { 
    setMethod("reduceRegions", siglist, function(x) {
        used <- logical(length(x@regions))
        used[x@anchor1] <- TRUE
        used[x@anchor2] <- TRUE
        new.dex <- integer(length(used))
        new.dex[used] <- seq_len(sum(used))
        x@anchor1 <- new.dex[x@anchor1]
        x@anchor2 <- new.dex[x@anchor2]
        x@regions <- x@regions[used]
        return(x)
    })
}

###############################################################

setGeneric("anchors<-", function(x, ..., value) { standardGeneric("anchors<-") })

anchor.repfun.gen <- function(is.GI) { 
    if (is.GI) { 
        type.arg <- GI.args
    } else {
        type.arg <- CM.args
    }
    type1 <- type.arg[2]

    function(x, type="both", ..., value) {
        type <- match.arg(type, type.arg)
        if (type=="both") { 
            if (length(value)!=2L) { 
                stop("'value' must be a list of 2 numeric vectors")
            }
            x@anchor1 <- as.integer(value[[1]])
            x@anchor2 <- as.integer(value[[2]])
        } else if (type==type1) {
            x@anchor1 <- as.integer(value)
        } else {
            x@anchor2 <- as.integer(value)
        }

        validObject(x)
        return(x)
    }
}

setReplaceMethod("anchors", "GInteractions", anchor.repfun.gen(TRUE))
setReplaceMethod("anchors", "ContactMatrix", anchor.repfun.gen(FALSE))

###############################################################
# Methods on InteractionSet that operate on GInteractions.

setGeneric("interactions", function(x, ...) { standardGeneric("interactions") })
setMethod("interactions", "InteractionSet", function(x) { return(x@interactions) })

setGeneric("interactions<-", function(x, value) { standardGeneric("interactions<-") })
setReplaceMethod("interactions", "InteractionSet", function(x, value) { 
    x@interactions <- value
    return(x)
}) 

setMethod("anchors", "InteractionSet", function(x, type="both", id=FALSE) { 
    anchors(x@interactions, type=type, id=id) 
})

setMethod("regions", "InteractionSet", function(x) { regions(x@interactions) })

setReplaceMethod("anchors", "InteractionSet", function(x, type="both", ..., value) { 
    anchors(x@interactions, type=type, ...) <- value 
    return(x)
})

setReplaceMethod("regions", "InteractionSet", function(x, value) { 
    regions(x@interactions) <- value
    return(x)
})

setReplaceMethod("replaceRegions", "InteractionSet", function(x, value) { 
    replaceRegions(x@interactions) <- value
    return(x)
})

setReplaceMethod("appendRegions", "InteractionSet", function(x, value) { 
    appendRegions(x@interactions) <- value
    return(x)
})

setMethod("reduceRegions", "InteractionSet", function(x) {
    x@interactions <- reduceRegions(x@interactions)
    return(x)
})

###############################################################
# Defining some other getters and setters.

setMethod("nrow", signature("GInteractions"), function(x) { 
    length(x) 
})

setMethod("$", "GInteractions", function(x, name) {
    return(x@elementMetadata[[name]])
})

setReplaceMethod("$", "GInteractions", function(x, name, value) {
    x@elementMetadata[[name]] <- value
    return(x)
})

setMethod("mcols", "InteractionSet", function(x, use.names=FALSE) {
    mcols(interactions(x), use.names=use.names)
})

setReplaceMethod("mcols", "InteractionSet", function(x, ..., value) {
    mcols(interactions(x), ...) <- value
    return(x)
})

setMethod("seqinfo", "GInteractions", function(x) {
     seqinfo(x@regions)
})

setReplaceMethod("seqinfo", "GInteractions", function(x, value) {
    seqinfo(x@regions) <- value
    validObject(x)
    return(x)
})

setMethod("seqinfo", "InteractionSet", function(x) {
     seqinfo(interactions(x))
})

setReplaceMethod("seqinfo", "InteractionSet", function(x, value) {
    seqinfo(interactions(x)) <- value
    return(x)
})

###############################################################
# End
