###############################################################
# Getters:

setGeneric("anchors", function(x, ...) standardGeneric("anchors"))
setGeneric("regions", function(x, ...) standardGeneric("regions"))

# A generating function, to capture differences in 'type' for 'anchors' call.
GI.args <- c("both", "first", "second") 
CM.args <- c("both", "row", "column") 
anchor.fun.gen <- function(is.GI) { 
    if (is.GI) { 
        type.arg <- GI.args
        n1fun <- n2fun <- names
    } else {
        type.arg <- CM.args
        n1fun <- rownames
        n2fun <- colnames
    }
    out.names <- type.arg[2:3]
    type1 <- out.names[1]

    function(x, type="both", id=FALSE) {
        type <- match.arg(type, type.arg)
        if (id) {
            if (type=="both") {
                out <- list(x@anchor1, x@anchor2)
                names(out[[1]]) <- n1fun(x)
                names(out[[2]]) <- n2fun(x)
                names(out) <- out.names
            } else if (type==type1) {
                out <- x@anchor1
                names(out) <- n1fun(x)
            } else {
                out <- x@anchor2
                names(out) <- n2fun(x)
            }
        } else {
            if (type=="both") {
                out <- GRangesList(x@regions[x@anchor1], x@regions[x@anchor2])
                names(out[[1]]) <- n1fun(x)
                names(out[[2]]) <- n2fun(x)
                names(out) <- out.names
            } else if (type==type1) {
                out <- x@regions[x@anchor1]
                names(out) <- n1fun(x)
            } else {
                out <-x@regions[x@anchor2]
                names(out) <- n2fun(x)
            }
        }
        return(out)
    }
}

# Defining the methods:
setMethod("anchors", "GInteractions", anchor.fun.gen(TRUE))
setMethod("anchors", "ContactMatrix", anchor.fun.gen(FALSE))

for (siglist in list("GInteractions", "ContactMatrix")) {
    setMethod("regions", siglist, function(x) { x@regions })
}

# Also defining some internal getters, for environment uses: 
setGeneric("anchor1", function(x) standardGeneric("anchor1"))
setGeneric("anchor2", function(x) standardGeneric("anchor2"))
setMethod("anchor1", "GInteractions", function(x) { x@anchor1 })
setMethod("anchor2", "GInteractions", function(x) { x@anchor2 })

###############################################################
# Setters:

setGeneric("regions<-", function(x, value) standardGeneric("regions<-"))
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

setGeneric("replaceRegions<-", function(x, value) standardGeneric("replaceRegions<-"))
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

setGeneric("appendRegions<-", function(x, value) standardGeneric("appendRegions<-"))
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

setGeneric("reduceRegions", function(x) standardGeneric("reduceRegions"))
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

setGeneric("anchors<-", function(x, ..., value) standardGeneric("anchors<-"))

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

setReplaceMethod("anchors", "StrictGInteractions", function(x, type="both", ..., value) {
    x <- as(x, "GInteractions")
    anchors(x, type=type, ...) <- value
    x <- swapAnchors(x)
    as(x, "StrictGInteractions")
})

setReplaceMethod("anchors", "ReverseStrictGInteractions", function(x, type="both", ..., value) {
    x <- as(x, "GInteractions")
    anchors(x, type=type, ...) <- value
    x <- swapAnchors(x, mode="reverse")
    as(x, "ReverseStrictGInteractions")
})

###############################################################
# Methods on InteractionSet that operate on GInteractions.

setGeneric("interactions", function(x, ...) standardGeneric("interactions"))
setMethod("interactions", "InteractionSet", function(x) { return(x@interactions) })

setGeneric("interactions<-", function(x, value) standardGeneric("interactions<-"))
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

###############################################################

setMethod("names", "GInteractions", function(x) { 
    x@NAMES 
})

setReplaceMethod("names", "GInteractions", function(x, value) {
    if (!is.null(value) && !is.character(value)) { value <- as.character(value) }                
    x@NAMES <- value
    validObject(x)
    return(x)
})

setMethod("names", "InteractionSet", function(x) {
    names(interactions(x))
})

setReplaceMethod("names", "InteractionSet", function(x, value) {
    names(interactions(x)) <- value
    return(x)
})

setMethod("dimnames", "ContactMatrix", function(x) {
    dimnames(x@matrix)
})

setReplaceMethod("dimnames", "ContactMatrix", function(x, value) {
    dimnames(x@matrix) <- value
    return(x)
})

###############################################################

for (siglist in c("GInteractions", "ContactMatrix")) { 
    setMethod("seqinfo", siglist, function(x) {
        seqinfo(x@regions)
    })
    
    setReplaceMethod("seqinfo", siglist, function(x, value) {
        seqinfo(x@regions) <- value
        validObject(x)
        return(x)
    })
}

setMethod("seqinfo", "InteractionSet", function(x) {
     seqinfo(interactions(x))
})

setReplaceMethod("seqinfo", "InteractionSet", function(x, value) {
    seqinfo(interactions(x)) <- value
    return(x)
})

##############################################
# Matrix dimensions

setMethod("dim", "ContactMatrix", function(x) { 
    dim(x@matrix)
})

setMethod("length", "ContactMatrix", function(x) { 
    length(x@matrix)
})


setMethod("as.matrix", "ContactMatrix", function(x) {
    return(x@matrix)
}) 

setGeneric("as.matrix<-", function(x, ..., value) standardGeneric("as.matrix<-"))
setReplaceMethod("as.matrix", "ContactMatrix", function(x, value) {
    if (is(value, "Matrix")) {
        if (!identical(dim(x), dim(value))) { 
            stop("replacement Matrix must have same dimensions as 'x'")
        }
        x@matrix <- value
    } else {
        x@matrix[] <- value
    }
    return(x)
}) 

###############################################################
# End
