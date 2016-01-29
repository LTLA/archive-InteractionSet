\name{ContactMatrix accessors}
\alias{ContactMatrix-access}

\alias{anchors,ContactMatrix-method}
\alias{anchors<-,ContactMatrix-method}

\alias{regions,ContactMatrix-method}
\alias{regions<-,ContactMatrix-method}

\alias{replaceRegions<-,ContactMatrix-method}
\alias{appendRegions<-,ContactMatrix-method}
\alias{reduceRegions,ContactMatrix-method}

\alias{show,ContactMatrix-method}
\alias{as.matrix,ContactMatrix-method}
\alias{t,ContactMatrix-method}

\alias{as.matrix<-}
\alias{as.matrix<-,ContactMatrix-method}

\alias{dim,ContactMatrix-method}
\alias{dimnames,ContactMatrix-method}
\alias{dimnames<-,ContactMatrix-method}
\alias{dimnames<-,ContactMatrix,ANY-method}
\alias{length,ContactMatrix-method}

\title{ContactMatrix accessors}

\description{
Methods to get and set fields in an ContactMatrix object.
} 

\usage{
\S4method{anchors}{ContactMatrix}(x, type="both", id=FALSE)
\S4method{anchors}{ContactMatrix}(x) <- value

\S4method{regions}{ContactMatrix}(x)
\S4method{regions}{ContactMatrix}(x) <- value
\S4method{replaceRegions}{ContactMatrix}(x) <- value
\S4method{appendRegions}{ContactMatrix}(x) <- value
\S4method{reduceRegions}{ContactMatrix}(x)

\S4method{show}{ContactMatrix}(x)
\S4method{as.matrix}{ContactMatrix}(x)
\S4method{as.matrix}{ContactMatrix}(x) <- value
\S4method{t}{ContactMatrix}(x)

\S4method{dim}{ContactMatrix}(x)
\S4method{dimnames}{ContactMatrix}(x)
\S4method{dimnames}{ContactMatrix}(x) <- value
\S4method{length}{ContactMatrix}(x)
}

\arguments{
\item{x}{A ContactMatrix object.}

\item{type}{a string specifying which anchors are to be extracted or replaced.}

\item{id}{a scalar logical indicating whether indices or GRanges should be returned.}

\item{value}{
    For \code{anchors<-} and \code{type="row"} or \code{"column"}, an integer vector of length equal to the number of rows or columns in \code{x}, respectively.
    For \code{type="both"}, a list of two such vectors must be supplied (first and second for rows and columns, respectively).

    For \code{anchors<-}, a list of two integer vectors when \code{type="both"}.
    The first and second vectors must have length equal to the number of rows and columns of \code{x}, respectively.
    For \code{type="row"} or \code{"column"}, only one vector needs to be supplied (rows or columns, respectively.

    For \code{regions<-}, a GRanges object of length equal to that of \code{regions(x)}.
    For \code{newRegions<-}, a GRanges object that is a superset of all entries in \code{regions(x)} involved in interactions.
    For \code{appendRegions<-}, a GRanges of any length containing additional regions.

    For \code{as.matrix<-}, a Matrix object of the same dimensions as that in the \code{matrix} slot.
    
    For \code{dimnames<-}, a list of two character vectors corresponding to the row and column names, respectively.
    These can also be passed separately via \code{rownames<-} and \code{colnames<-}.
}
}

\details{
The return value of \code{anchors} varies depending on \code{type} and \code{id}.
If \code{id=FALSE}, a GRanges object will be returned containing the regions specified by the \code{anchor1} or \code{anchor2} slots in \code{x}, 
    for \code{type=="row"} or \code{"column"}, respectively.
If \code{type="both"}, a GRangesList will be returned with two entries - \code{row} and \code{column}, containing regions specified by \code{anchor1} and \code{anchor2} respectively.
However, if \code{id=TRUE}, the integer vectors in the \code{anchor1} or \code{anchor2} slots of \code{x} are returned directly 
    (a list of length two is returned for \code{type="both"}, containing both of the vectors).

Replacement in \code{anchors<-} can only be performed using anchor indices.
If \code{type="both"}, a list of two integer vectors is required in \code{value}, specifying the indices of the row- and column-wise interacting regions in \code{x}.
If \code{type="row"} or \code{"column"}, an integer vector is required to replace the existing row- or column-wise indices in the \code{anchor1} or \code{anchor2} slot, respectively.
   
For \code{regions}, a GRanges is returned equal to the \code{regions} slot in \code{x}.
For \code{regions<-}, the GRanges in \code{value} is used to replace the \code{regions} slot.
Resorting of the replacement GRanges is performed automatically, along with re-indexing of the anchors.    
In addition, the input GRanges must be of the same length as the existing object in the \code{regions} slot.
The \code{newRegions} replacement method can take variable length GRanges, but requires that the replacement contain (at least) all ranges contained in \code{anchors(x)}.
The \code{appendRegions} replacement method appends extra intervals to the existing \code{regions} slot of \code{x}.
The \code{reduceRegions} method removes unused entries in the \code{regions} slot, to save memory -- see \code{\link{reduceRegions,GInteractions-method}} for more details.

The \code{show} method will print out various details of the object, such as the dimensions of the \code{matrix} slot and the length of the \code{regions} slot.
The \code{as.matrix} method will return the value of the \code{matrix} slot, containing the Matrix of interaction data.
Replacement with a Matrix of the same dimensions can be performed using the \code{as.matrix<-} function.
The \code{t} method will transpose the matrix, i.e., switch the rows and columns (and switch the vectors in the \code{anchor1} and \code{anchor2} slots).

The \code{dim} method will return a vector of length 2, containing the dimensions of the \code{matrix} slot.
The \code{dimnames} method will return a list of length 2, containing the row and column names of \code{matrix} (these can be modified with the \code{dimnames<-} method).
The \code{length} method will return an integer scalar corresponding to the total number of entries in the \code{matrix} slot.

As the ContactMatrix class inherits from the \code{\link{Annotated}} class, additional metadata can be stored in the \code{metadata} slot.
This can be accessed or modified with \code{\link{metadata,Annotated-method}}.
}

\value{
For the getters, values in various slots of \code{x} are returned, while for the setters, the slots of \code{x} are modified accordingly -- see Details.
}

\examples{
example(ContactMatrix, echo=FALSE) # Generate a nice object.
show(x)

# Various matrix methods:
as.matrix(x)
t(x)

nrow(x)
ncol(x)
length(x)

# Accessing anchor ranges or indices:
anchors(x)
anchors(x, type="row")
anchors(x, id=TRUE)

anchors(x, id=TRUE, type="row")
anchors(x, id=TRUE, type="column")

# Modifying anchor indices:
nregs <- length(regions(x))
anchors(x) <- list(sample(nregs, nrow(x), replace=TRUE),
                   sample(nregs, ncol(x), replace=TRUE))
anchors(x, id=TRUE, type="row")
anchors(x, id=TRUE, type="column")

# Accessing or modifying regions:
regions(x)
regions(x)$score <- runif(length(regions(x)))

new.ranges <- c(regions(x), resize(regions(x), fix="center", width=50))
try(regions(x) <- new.ranges) # Fails
replaceRegions(x) <- new.ranges # Succeeds

length(regions(x))
appendRegions(x) <- GRanges("chrA", IRanges(5:10+1000, 1100+5:10), score=runif(6))
length(regions(x))

reduceRegions(x)

# Setting metadata
metadata(x)$name <- "I am a contact matrix"
metadata(x)
}

\author{
Aaron Lun
}

\seealso{
\code{\link{ContactMatrix-class}}
}
