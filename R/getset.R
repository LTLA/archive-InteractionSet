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
                out <- list(anchor1(x), anchor2(x))
                names(out[[1]]) <- n1fun(x) 
                names(out[[2]]) <- n2fun(x)
                names(out) <- out.names
            } else if (type==type1) {
                out <- anchor1(x)
                names(out) <- n1fun(x)
            } else {
                out <- anchor2(x)
                names(out) <- n2fun(x)
            }
        } else {
            if (type=="both") {
                out <- GRangesList(regions(x)[anchor1(x)], regions(x)[anchor2(x)])
                names(out[[1]]) <- n1fun(x)
                names(out[[2]]) <- n2fun(x)
                names(out) <- out.names
            } else if (type==type1) {
                out <- regions(x)[anchor1(x)]
                names(out) <- n1fun(x)
            } else {
                out <- regions(x)[anchor2(x)]
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

# Defining some convenience methods for GInteractions.
setMethod("first", "GInteractions", function(x) { anchors(x, type="first") })
setMethod("second", "GInteractions", function(x) { anchors(x, type="second") })

# Also defining some internal getters, for environment uses: 
setGeneric("anchor1", function(x) standardGeneric("anchor1"))
setGeneric("anchor2", function(x) standardGeneric("anchor2"))
for (siglist in list("GInteractions", "ContactMatrix")) { 
    setMethod("anchor1", siglist, function(x) { x@anchor1 })
    setMethod("anchor2", siglist, function(x) { x@anchor2 })
}

###############################################################
# Setters:

setGeneric("regions<-", function(x, value) standardGeneric("regions<-"))
for (siglist in c("GInteractions", "ContactMatrix")) { 
    setReplaceMethod("regions", siglist, function(x, value) {
        if (length(value)!=length(regions(x))) { 
            stop("assigned value must be of the same length as 'regions(x)'")
        }
        out <- .resort_regions(anchor1(x), anchor2(x), value)
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
        new.a1 <- converter[anchor1(x)]
        new.a2 <- converter[anchor2(x)]
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
        out <- .resort_regions(anchor1(x), anchor2(x), c(regions(x), value)) 
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
        used <- logical(length(regions(x)))
        used[anchor1(x)] <- TRUE
        used[anchor2(x)] <- TRUE
        new.dex <- integer(length(used))
        new.dex[used] <- seq_len(sum(used))
        x@anchor1 <- new.dex[anchor1(x)]
        x@anchor2 <- new.dex[anchor2(x)]
        x@regions <- regions(x)[used]
        return(x)
    })
}

###############################################################

setGeneric("anchorIds<-", function(x, ..., value) standardGeneric("anchorIds<-"))

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

setReplaceMethod("anchorIds", "GInteractions", anchor.repfun.gen(TRUE))
setReplaceMethod("anchorIds", "ContactMatrix", anchor.repfun.gen(FALSE))

setReplaceMethod("anchorIds", "StrictGInteractions", function(x, type="both", ..., value) {
    x <- as(x, "GInteractions")
    anchorIds(x, type=type, ...) <- value
    x <- swapAnchors(x)
    as(x, "StrictGInteractions")
})

setReplaceMethod("anchorIds", "ReverseStrictGInteractions", function(x, type="both", ..., value) {
    x <- as(x, "GInteractions")
    anchorIds(x, type=type, ...) <- value
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
    anchors(interactions(x), type=type, id=id) 
})

setMethod("first", "InteractionSet", function(x) { anchors(x, type="first") })

setMethod("second", "InteractionSet", function(x) { anchors(x, type="second") })

setMethod("regions", "InteractionSet", function(x) { regions(interactions(x)) })

setReplaceMethod("anchorIds", "InteractionSet", function(x, type="both", ..., value) { 
    anchorIds(interactions(x), type=type, ...) <- value 
    return(x)
})

setReplaceMethod("regions", "InteractionSet", function(x, value) { 
    regions(interactions(x)) <- value
    return(x)
})

setReplaceMethod("replaceRegions", "InteractionSet", function(x, value) { 
    replaceRegions(interactions(x)) <- value
    return(x)
})

setReplaceMethod("appendRegions", "InteractionSet", function(x, value) { 
    appendRegions(interactions(x)) <- value
    return(x)
})

setMethod("reduceRegions", "InteractionSet", function(x) {
    interactions(x) <- reduceRegions(interactions(x))
    return(x)
})

###############################################################
# Defining some other getters and setters.

setMethod("$", "GInteractions", function(x, name) {
    return(mcols(x)[[name]])
})

setReplaceMethod("$", "GInteractions", function(x, name, value) {
    mcols(x)[[name]] <- value
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
# Name getting and setting.

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
    dimnames(as.matrix(x))
})

setReplaceMethod("dimnames", "ContactMatrix", function(x, value) {
    dimnames(as.matrix(x)) <- value
    return(x)
})

###############################################################
# Seqinfo getting and setting.

for (siglist in c("GInteractions", "ContactMatrix")) { 
    setMethod("seqinfo", siglist, function(x) {
        seqinfo(regions(x))
    })
    
    setReplaceMethod("seqinfo", siglist, function(x, value) {
        seqinfo(regions(x)) <- value
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

###############################################################
# Deprecation of anchors<- setter, use anchorIds<- instead
# (to avoid confusion with anchors(x), which returns GRanges by default).

setGeneric("anchors<-", function(x, ..., value) standardGeneric("anchors<-"))
for (siglist in c("GInteractions", "InteractionSet", "ContactMatrix")) { 
    setReplaceMethod("anchors", siglist, function(x, type="both", ..., value) {
        .Deprecated("anchorIds<-", old="anchors<-")    
        anchorIds(x, type=type, ...) <- value
        return(x)
    })
}

##############################################
# Matrix dimensions

setMethod("dim", "ContactMatrix", function(x) { 
    dim(as.matrix(x))
})

setMethod("length", "ContactMatrix", function(x) { 
    length(as.matrix(x))
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
