# getset.R

setGeneric("anchors", function(x, ...) standardGeneric("anchors"))
setGeneric("anchor1", function(x) standardGeneric("anchor1"))
setGeneric("anchor2", function(x) standardGeneric("anchor2"))
setGeneric("unchecked_anchor1<-", function(x, value) standardGeneric("unchecked_anchor1<-"))
setGeneric("unchecked_anchor2<-", function(x, value) standardGeneric("unchecked_anchor2<-"))

setGeneric("regions", function(x, ...) standardGeneric("regions"))
setGeneric("unchecked_regions<-", function(x, value) standardGeneric("unchecked_regions<-"))
setGeneric("regions<-", function(x, value) standardGeneric("regions<-"))
setGeneric("replaceRegions<-", function(x, value) standardGeneric("replaceRegions<-"))
setGeneric("appendRegions<-", function(x, value) standardGeneric("appendRegions<-"))
setGeneric("reduceRegions", function(x) standardGeneric("reduceRegions"))

setGeneric("anchorIds<-", function(x, ..., value) standardGeneric("anchorIds<-"))
setGeneric("anchors<-", function(x, ..., value) standardGeneric("anchors<-"))

setGeneric("interactions", function(x, ...) standardGeneric("interactions"))
setGeneric("unchecked_interactions<-", function(x, value) standardGeneric("unchecked_interactions<-"))
setGeneric("interactions<-", function(x, value) standardGeneric("interactions<-"))

setGeneric("unchecked_matrix<-", function(x, value) standardGeneric("unchecked_matrix<-"))
setGeneric("as.matrix<-", function(x, ..., value) standardGeneric("as.matrix<-"))

# boundingBox.R

setGeneric("boundingBox", function(x, f) standardGeneric("boundingBox"))

# ContactMatrix-class.R

setGeneric("ContactMatrix", function(matrix, anchor1, anchor2, regions, ...) standardGeneric("ContactMatrix"))

# conversion.R

setGeneric("inflate", function(x, ...) standardGeneric("inflate"))
setGeneric("deflate", function(x, ...) standardGeneric("deflate"))

# distances.R

setGeneric("pairdist", function(x, ...) standardGeneric("pairdist"))
setGeneric("intrachr", function(x) standardGeneric("intrachr"))

# GInteractions-class.R

setGeneric("GInteractions", function(anchor1, anchor2, regions, ...) standardGeneric("GInteractions"))

# InteractionSet-class.R

setGeneric("InteractionSet", function(assays, interactions, ...) standardGeneric("InteractionSet"))

# linearize.R

setGeneric("linearize", function(x, ref, ...) standardGeneric("linearize"))

# linkOverlaps.R

setGeneric("linkOverlaps", function(query, subject1, subject2, ...) standardGeneric('linkOverlaps'))

# pairs.R

setGeneric("pairs", function(x, ...) standardGeneric("pairs"))

# swapAnchors.R

setGeneric("swapAnchors", function(x, ...) { standardGeneric("swapAnchors") })

