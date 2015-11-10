##############################################
# Defines the ContactMatrix class.

setClass("ContactMatrix", 
    contains="matrix",
    slots=list(
        anchor1="integer",
        anchor2="integer",
        regions="GRanges"
    )		
)

setValidity("ContactMatrix", function(object) {
    if (nrow(object)!=length(object@anchor1)) { 
        return("'assays' nrow must be equal to length of 'anchor1'")
    }
    if (ncol(object)!=length(object@anchor2)) {
        return("'assays' ncol must be equal to length of 'anchor2'")
    }

    if (is.unsorted(object@regions)) { 
        return("'regions' should be sorted")
    }
    msg <- .check_inputs(object@anchor1, object@anchor2, object@regions)
    if (is.character(msg)) { return(msg) }
   
    return(TRUE)
}) 

setMethod("show", signature("ContactMatrix"), function(object) {
    cat("class:", class(object), "\n")
    cat("dim:", dim(object), "\n")
    cat(sprintf("regions: %i\n", length(object@regions)))
})

##############################################
# Constructor:

.new_ContactMatrix <- function(mat, anchor1, anchor2, regions) {
    anchor1 <- as.integer(anchor1)
    anchor2 <- as.integer(anchor2)
    
    msg <- .check_inputs(anchor1, anchor2, regions)
    if (is.character(msg)) { stop(msg) }
    out <- .resort_regions(anchor1, anchor2, regions)

    new("ContactMatrix", mat, anchor1=out$anchor1, anchor2=out$anchor2, regions=out$regions)
}

setGeneric("ContactMatrix", function(mat, anchor1, anchor2, ...) { standardGeneric("ContactMatrix") })
setMethod("ContactMatrix", c("ANY", "numeric", "numeric"), .new_ContactMatrix)

setMethod("ContactMatrix", c("ANY", "GRanges", "GRanges"), 
    function(mat, anchor1, anchor2, regions) { 

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
        
        .new_ContactMatrix(mat, anchor1, anchor2, regions)
    }
)

##############################################
# End
