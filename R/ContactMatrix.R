##############################################
# Defines the ContactMatrix class.

setClass("ContactMatrix", 
    slots=list(
        matrix="matrix", 
        anchor1="integer",
        anchor2="integer",
        regions="GRanges"
    )		
)

setValidity("ContactMatrix", function(object) {
    if (is.unsorted(object@regions)) {
        return("'regions' should be sorted")
    }
    msg <- .check_inputs(object@anchor1, object@anchor2, object@regions, same.length=FALSE)
    if (is.character(msg)) { return(msg) }

    if (nrow(object@matrix)!=length(object@anchor1)) { 
        return("'matrix' nrow must be equal to length of 'anchor1'")
    }
    if (ncol(object@matrix)!=length(object@anchor2)) {
        return("'matrix' ncol must be equal to length of 'anchor2'")
    }
    return(TRUE)
}) 

setMethod("show", signature("ContactMatrix"), function(object) {
    cat("class:", class(object), "\n")
    cat("dim:", dim(object@matrix), "\n")
    cat(sprintf("regions: %i\n", length(object@regions)))
})

##############################################
# Constructor:

.new_ContactMatrix <- function(matrix, anchor1, anchor2, regions) {
    anchor1 <- as.integer(anchor1)
    anchor2 <- as.integer(anchor2)
    
    msg <- .check_inputs(anchor1, anchor2, regions, same.length=FALSE)
    if (is.character(msg)) { stop(msg) }
    out <- .resort_regions(anchor1, anchor2, regions, enforce.order=FALSE)

    new("ContactMatrix", matrix=matrix, anchor1=out$anchor1, anchor2=out$anchor2, regions=out$regions)
}

setGeneric("ContactMatrix", function(matrix, anchor1, anchor2, ...) { standardGeneric("ContactMatrix") })
setMethod("ContactMatrix", c("ANY", "numeric", "numeric"), .new_ContactMatrix)

setMethod("ContactMatrix", c("ANY", "GRanges", "GRanges"), 
    function(matrix, anchor1, anchor2, regions) { 

        if (missing(regions)) { 
            collated <- .collate_GRanges(anchor1, anchor2)
            regions <- collated$ranges
            anchor1 <- collated$indices[[1]]
            anchor2 <- collated$indices[[2]]
        } else {
            anchor1 <- match(anchor1, regions)
            anchor2 <- match(anchor2, regions)
            if (any(is.na(anchor1)) || any(is.na(anchor2))) {
                 stop("anchor regions missing in specified 'regions'")
            }
        }
        
        .new_ContactMatrix(matrix, anchor1, anchor2, regions)
    }
)

##############################################
# Subsetting

setMethod("[", c("ContactMatrix", "ANY", "ANY"), function(x, i, j, ..., drop=TRUE) {
    if (!missing(i)) { 
        x@anchor1 <- x@anchor1
    }
    if (!missing(j)) {
        x@anchor2 <- x@anchor2[j]
    }
    x@matrix <- x@matrix[i,j,drop=FALSE]
    return(x)
}) 

setMethod("subset", "ContactMatrix", function(x, i, j) {
    x[i,j]
})

##############################################
# Other getters

setMethod("as.matrix", "ContactMatrix", function(x) {
    return(x@matrix)
}) 

##############################################
# End


# Testing:
# cm <- InteractionSet:::ContactMatrix(matrix(0, 4, 4), 1:4, 1:4, GRanges("chr1", IRanges(1:4, 1:4)))

