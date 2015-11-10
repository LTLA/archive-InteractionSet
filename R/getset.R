###############################################################
# Getters and setters

setGeneric("anchors", function(x, ...) { standardGeneric("anchors") })
setGeneric("regions", function(x, ...) { standardGeneric("regions") })
setGeneric("regions<-", function(x, value) { standardGeneric("regions<-") })
setGeneric("anchors<-", function(x, ..., value) { standardGeneric("anchors<-") })

# A generating function, to capture differences in 'type' for 'anchors' call.
anchorfun.gen <- function(is.IS) { 
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

# Defining the methods
for (siglist in list(
        "InteractionSet",
        "ContactMatrix"
    )) {

    setMethod("anchors", siglist, anchorfun.gen(siglist=="InteractionSet"))
    setMethod("regions", siglist, function(x) { x@regions })

    setReplaceMethod("regions", siglist, function(x, value) {
        stopifnot(length(value)==length(regions(x)))
        out <- .resort_regions(x@anchor1, x@anchor2, value)
        x@anchor1 <- out$anchor1
        x@anchor2 <- out$anchor2
        x@regions <- out$regions
        validObject(x)
        return(x)
    })

    setReplaceMethod("anchors", siglist, function(x, value) {
        if (length(value)!=2L) { stop("'value' must be a list of 2 numeric vectors") }
        if (length(value[[1]])!=length(value[[2]])) { stop("vectors in 'value' must be of the same length") }
        first <- as.integer(value[[1]])
        second <- as.integer(value[[2]])
        if (!all(is.finite(first)) || !all(is.finite(second))) { stop("all anchor indices should be finite") }

        out <- .enforce_order(first, second)
        x@anchor1 <- out$anchor1
        x@anchor2 <- out$anchor2
        validObject(x)
        return(x)
    })
}

