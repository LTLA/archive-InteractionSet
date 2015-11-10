\name{InteractionSet accessors}
\alias{InteractionSet-access}

\alias{anchors}
\alias{anchors,InteractionSet-method}
\alias{anchors<-}
\alias{anchors<-,InteractionSet-method}
\alias{regions}
\alias{regions,InteractionSet-method}
\alias{regions<-}
\alias{regions<-,InteractionSet-method}

\alias{show}
\alias{show,InteractionSet-method}

\title{InteractionSet accessors}

\description{
Methods to get and set fields in an InteractionSet object.
} 

\usage{
\S4method{anchors}{InteractionSet}(x, type="both", id=FALSE)
\S4method{anchors}{InteractionSet}(x) <- value

\S4method{regions}{InteractionSet}(x)
\S4method{regions}{InteractionSet}(x) <- value

\S4method{show}{InteractionSet}(x)
}

\arguments{
\item{x}{An InteractionSet object.}

\item{type}{a string specifying which anchors are to be extracted or replaced.}

\item{id}{a scalar logical indicating whether indices or GRanges should be returned.}

\item{value}{
    For \code{anchors<-,InteractionSet-method}, a list of two integer vectors if \code{type="both"}, and an integer vector otherwise.
    All vectors must be of length equal to the number of rows in \code{x}.

    For \code{regions<-,InteractionSet-method}, a GRanges object of length equal to that of \code{regions(x)}.
}
}

\details{
The return value of \code{anchors,InteractionSet-method} varies depending on \code{type} and \code{id}.
If \code{id=FALSE}, a GRanges object will be returned containing the regions specified by the \code{anchor1} or \code{anchor2} slots in \code{x}, 
    for \code{type=="first"} or \code{"second"}, respectively.
If \code{type="both"}, a GRangesList will be returned with two entries - \code{first} and \code{second}, containing regions specified by \code{anchor1} and \code{anchor2} respectively.
    
However, if \code{id=TRUE}, the integer vectors in the \code{anchor1} or \code{anchor2} slots of \code{x} are returned directly 
    (a list of length two is returned for \code{type="both"}, containing both of the vectors).
Replacement in \code{anchors<-,InteractionSet-method} can only be performed using a list of two integer vectors in \code{value}.
Swapping is performed automatically to ensure that the entries in the \code{anchor1} slot are not less than those in \code{anchor2} in the modified object.
   
For \code{regions,InteractionSet-method}, a GRanges is returned equal to the \code{regions} slot in \code{x}.
For \code{regions<-,InteractionSet-method}, resorting of the input GRanges in \code{value} is performed automatically, along with re-indexing of the anchors.    
In addition, the input GRanges must be of the same length as the existing object in the \code{regions} slot.

The \code{show} method will print information equivalent to that done for a SummarizedExperiment object.
An additional line is added indicating the number of regions in the \code{regions} slot of the object.

For access and setting of all other slots, see \code{\link{SummarizedExperiment0-class}} for details on the appropriate methods.
This includes \code{assays}, \code{colData}, \code{mcol} or \code{metadata}, which can be applied to InteractionSet objects.
}

\examples{
example(InteractionSet, echo=FALSE) # Generate a nice object.
show(x)

anchors(x)
anchors(x, type="first")
anchors(x, id=TRUE)

anchors(x, id=TRUE, type="first")
anchors(x, id=TRUE, type="second")
anchors(x) <- list(sample(length(regions(x)), nrow(x), replace=TRUE),
    sample(length(regions(x)), nrow(x), replace=TRUE))
anchors(x, id=TRUE, type="first")
anchors(x, id=TRUE, type="second")

regions(x)
regions(x)$score <- runif(length(regions(x)))
}

\author{
Aaron Lun
}

\seealso{
\code{\link{InteractionSet-class}},
\code{\link{SummarizedExperiment0-class}}
}
