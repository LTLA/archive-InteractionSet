###############################################################
# Getters:

setGeneric("anchors", function(x, ...) { standardGeneric("anchors") })
setGeneric("regions", function(x, ...) { standardGeneric("regions") })

# A generating function, to capture differences in 'type' for 'anchors' call.
anchor.fun.gen <- function(is.IS) { 
    if (is.IS) { 
        type.arg <- c("both", "first", "second") 
    } else {
        type.arg <- c("both", "row", "column") 
    }
    out.names <- type.arg[2:3]

    function(x, type="both", id=FALSE) {
        type <- match.arg(type, type.arg)
        if (id) {
            if (type=="both") {
                out <- list(x@anchor1, x@anchor2)
                names(out) <- out.names
                return(out)
            } else if (type=="first") {
                return(x@anchor1)
            } else {
                return(x@anchor2)
            }
        } else {
            if (type=="both") {
                out <- GRangesList(x@regions[x@anchor1], x@regions[x@anchor2])
                names(out) <- out.names
                return(out)
            } else if (type=="first") {
                return(x@regions[x@anchor1])
            } else {
                return(x@regions[x@anchor2])
            }
        }
    }
}

# Defining the methods:
setMethod("anchors", "InteractionSet", anchor.fun.gen(TRUE))
setMethod("anchors", "ContactMatrix", anchor.fun.gen(FALSE))

for (siglist in list("InteractionSet", "ContactMatrix")) {
    setMethod("regions", siglist, function(x) { x@regions })
}

###############################################################
# Setters:

setGeneric("regions<-", function(x, value) { standardGeneric("regions<-") })

region.fun.gen <- function(enforce.order) {
    function(x, value) {
        if (length(value)!=length(regions(x))) { 
            stop("assigned value must be of the same length as 'regions(x)'")
        }
        out <- .resort_regions(x@anchor1, x@anchor2, value, enforce.order=enforce.order)
        x@anchor1 <- out$anchor1
        x@anchor2 <- out$anchor2
        x@regions <- out$regions
        validObject(x)
        return(x)
    }
}

setReplaceMethod("regions", "InteractionSet", region.fun.gen(TRUE))
setReplaceMethod("regions", "ContactMatrix", region.fun.gen(FALSE))

# Also allow setting of regions of different length.

setGeneric("newRegions<-", function(x, value) { standardGeneric("newRegions<-") })

newRegion.fun.gen <- function(enforce.order) {
    function(x, value) {
        converter <- match(regions(x), value)
        new.a1 <- converter[x@anchor1]
        new.a2 <- converter[x@anchor2]
        if (any(is.na(new.a1)) || any(is.na(new.a2))) { 
            stop("some existing ranges do not exist in replacement GRanges") 
        }

        out <- .resort_regions(new.a1, new.a2, value, enforce.order=enforce.order)
        x@anchor1 <- out$anchor1
        x@anchor2 <- out$anchor2
        x@regions <- out$regions
        return(x)
    }
}

setReplaceMethod("newRegions", "InteractionSet", newRegion.fun.gen(TRUE))
setReplaceMethod("newRegions", "ContactMatrix", newRegion.fun.gen(FALSE))

###############################################################

# 'anchors<-' is necessarily different between classes, as there is no requirement for equal length 'anchor1' and 'anchor2'
# nor is there any need to enforce 'anchor1 >= anchor2' in "ContactMatrix" objects.

setGeneric("anchors<-", function(x, ..., value) { standardGeneric("anchors<-") })

setReplaceMethod("anchors", "InteractionSet", function(x, value) {
    if (length(value)!=2L) { 
        stop("'value' must be a list of 2 numeric vectors")
    }
    first <- as.integer(value[[1]])
    second <- as.integer(value[[2]])
    msg <- .check_inputs(first, second, regions(x)) # Need check here, otherwise .enforce_order might fail.
    if (is.character(msg)) { stop(msg) }

    out <- .enforce_order(first, second)
    x@anchor1 <- out$anchor1
    x@anchor2 <- out$anchor2
    validObject(x)
    return(x)
})

setReplaceMethod("anchors", "ContactMatrix", function(x, value) {
    if (length(value)!=2L) { 
        stop("'value' must be a list of 2 numeric vectors")
    }
    x@anchor1 <- as.integer(value[[1]])
    x@anchor2 <- as.integer(value[[2]])
    validObject(x)
    return(x)
})

###############################################################
# End
