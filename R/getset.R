###############################################################
# Simple getters and setters (some of these are exported):

for (siglist in list("GInteractions", "ContactMatrix")) {
    setMethod("regions", siglist, function(x) { x@regions })
}

for (siglist in list("GInteractions", "ContactMatrix")) { 
    # internal, for use with environments and with().
    setMethod("anchor1", siglist, function(x) { x@anchor1 })
    setMethod("anchor2", siglist, function(x) { x@anchor2 })

    setReplaceMethod("unchecked_regions", siglist, function(x, value) {
        x@regions <- value
        return(x)        
    })

    setReplaceMethod("unchecked_anchor1", siglist, function(x, value) {
        x@anchor1 <- value 
        return(x)        
    })

    setReplaceMethod("unchecked_anchor2", siglist, function(x, value) {
        x@anchor2 <- value 
        return(x)        
    })
}

setMethod("interactions", "InteractionSet", function(x) { return(x@interactions) })

setReplaceMethod("unchecked_interactions", "InteractionSet", function(x, value) { 
    x@interactions <- value
    return(x)
}) 

setMethod("as.matrix", "ContactMatrix", function(x) {
    return(x@matrix)
}) 

setReplaceMethod("unchecked_matrix", "ContactMatrix", function(x, value) {
    x@matrix <- value
    return(x)
}) 

###############################################################
# Exported getters for anchors. Using a generating function, 
# to capture differences in 'type' for 'anchors' call.

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

# Defining some convenience methods for GInteractions.
setMethod("first", "GInteractions", function(x) { anchors(x, type="first") })
setMethod("second", "GInteractions", function(x) { anchors(x, type="second") })

###############################################################
# Setters for regions:

for (siglist in c("GInteractions", "ContactMatrix")) { 
    setReplaceMethod("regions", siglist, function(x, value) {
        if (length(value)!=length(regions(x))) { 
            stop("assigned value must be of the same length as 'regions(x)'")
        }
        out <- .resort_regions(anchor1(x), anchor2(x), value)
        unchecked_anchor1(x) <- out$anchor1
        unchecked_anchor2(x) <- out$anchor2
        unchecked_regions(x) <- out$regions
        validObject(x)
        return(x)
    })
}

# Also allow setting of regions of different length.

for (siglist in c("GInteractions", "ContactMatrix")) { 
    setReplaceMethod("replaceRegions", siglist, function(x, value) {
        converter <- match(regions(x), value)
        new.a1 <- converter[anchor1(x)]
        new.a2 <- converter[anchor2(x)]
        if (any(is.na(new.a1)) || any(is.na(new.a2))) { 
            stop("some existing ranges do not exist in replacement GRanges") 
        }

        out <- .resort_regions(new.a1, new.a2, value)
        unchecked_anchor1(x) <- out$anchor1
        unchecked_anchor2(x) <- out$anchor2
        unchecked_regions(x) <- out$regions
        return(x)
    })
}

# Append regions.

for (siglist in c("GInteractions", "ContactMatrix")) {
    setReplaceMethod("appendRegions", siglist, function(x, value) {
        out <- .resort_regions(anchor1(x), anchor2(x), c(regions(x), value)) 
        unchecked_anchor1(x) <- out$anchor1
        unchecked_anchor2(x) <- out$anchor2
        unchecked_regions(x) <- out$regions
        return(x)
    })
}

# Reduce regions.

for (siglist in c("GInteractions", "ContactMatrix")) { 
    setMethod("reduceRegions", siglist, function(x) {
        used <- logical(length(regions(x)))
        used[anchor1(x)] <- TRUE
        used[anchor2(x)] <- TRUE
        new.dex <- integer(length(used))
        new.dex[used] <- seq_len(sum(used))
        unchecked_anchor1(x) <- new.dex[anchor1(x)]
        unchecked_anchor2(x) <- new.dex[anchor2(x)]
        unchecked_regions(x) <- regions(x)[used]
        return(x)
    })
}

###############################################################
# Setters for anchors.

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
            unchecked_anchor1(x) <- as.integer(value[[1]])
            unchecked_anchor2(x) <- as.integer(value[[2]])
        } else if (type==type1) {
            unchecked_anchor1(x) <- as.integer(value)
        } else {
            unchecked_anchor2(x) <- as.integer(value)
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

setReplaceMethod("interactions", "InteractionSet", function(x, value) { 
    unchecked_interactions(x) <- value
    validObject(x)
    return(x)
}) 

setMethod("anchors", "InteractionSet", function(x, type="both", id=FALSE) { 
    anchors(interactions(x), type=type, id=id) 
})

setMethod("first", "InteractionSet", function(x) { anchors(x, type="first") })

setMethod("second", "InteractionSet", function(x) { anchors(x, type="second") })

setMethod("regions", "InteractionSet", function(x) { regions(interactions(x)) })

# Modification of regions in GInteractions doesn't affect validity of InteractionSet.

setReplaceMethod("anchorIds", "InteractionSet", function(x, type="both", ..., value) { 
    i <- interactions(x)
    anchorIds(i, type=type, ...) <- value 
    unchecked_interactions(x) <- i
    return(x)
})

setReplaceMethod("regions", "InteractionSet", function(x, value) { 
    i <- interactions(x)
    regions(i) <- value
    unchecked_interactions(x) <- i
    return(x)
})

setReplaceMethod("replaceRegions", "InteractionSet", function(x, value) { 
    i <- interactions(x)
    replaceRegions(i) <- value
    unchecked_interactions(x) <- i
    return(x)
})

setReplaceMethod("appendRegions", "InteractionSet", function(x, value) { 
    i <- interactions(x)
    appendRegions(i) <- value
    unchecked_interactions(x) <- i
    return(x)
})

setMethod("reduceRegions", "InteractionSet", function(x) {
    unchecked_interactions(x) <- reduceRegions(interactions(x))
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
    i <- interactions(x)                 
    mcols(i, ...) <- value
    unchecked_interactions(x) <- i
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
    i <- interactions(x)                 
    names(i) <- value
    unchecked_interactions(x) <- i
    return(x)
})

setMethod("dimnames", "ContactMatrix", function(x) {
    dimnames(as.matrix(x))
})

setReplaceMethod("dimnames", "ContactMatrix", function(x, value) {
    m <- as.matrix(x)
    dimnames(m) <- value
    unchecked_matrix(x) <- m
    return(x)
})

###############################################################
# Seqinfo getting and setting.

for (siglist in c("GInteractions", "ContactMatrix")) { 
    setMethod("seqinfo", siglist, function(x) {
        seqinfo(regions(x))
    })
    
    setReplaceMethod("seqinfo", siglist, function(x, value) {
        r <- regions(x)
        seqinfo(r) <- value
        unchecked_regions(x) <- r
        return(x)
    })
}

setMethod("seqinfo", "InteractionSet", function(x) {
     seqinfo(interactions(x))
})

setReplaceMethod("seqinfo", "InteractionSet", function(x, value) {
    i <- interactions(x)                 
    seqinfo(i) <- value
    unchecked_interactions(x) <- i
    return(x)
})

###############################################################
# Deprecation of anchors<- setter, use anchorIds<- instead
# (to avoid confusion with anchors(x), which returns GRanges by default).

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

setReplaceMethod("as.matrix", "ContactMatrix", function(x, value) {
    if (is(value, "Matrix")) {
        if (!identical(dim(x), dim(value))) { 
            stop("replacement Matrix must have same dimensions as 'x'")
        }
        unchecked_matrix(x) <- value
    } else {
        x@matrix[] <- value
    }
    return(x)
}) 

###############################################################
# End
