\name{Convert classes}

\alias{inflate}
\alias{inflate,GInteractions-method}
\alias{inflate,InteractionSet-method}
\alias{deflate}
\alias{deflate,ContactMatrix-method}

\title{Convert between classes}

\description{
Inflate a GInteractions or InteractionSet into a ContactMatrix, or deflate a ContactMatrix to an InteractionSet.
}

\usage{
\S4method{inflate}{GInteractions}(x, rows, columns, fill, swap=TRUE, sparse=FALSE, ...)

\S4method{inflate}{InteractionSet}(x, rows, columns, assay=1L, sample=1L, fill=NULL, swap=TRUE, sparse=FALSE, ...)

\S4method{deflate}{ContactMatrix}(x, collapse=TRUE, extract, use.zero, use.na, ...)
}

\arguments{
\item{x}{
    A GInteractions or InteractionSet object for \code{inflate}, or a ContactMatrix object for \code{deflate}.
}

\item{rows, columns}{
    An integer, logical or character vector, a GRanges object or \code{NULL}, indicating the regions of interest to be used as the rows or columns of the ContactMatrix.
}

\item{assay}{
    A numeric scalar indicating the assay of the InteractionSet object, from which values are extracted to fill the ContactMatrix.
}

\item{sample}{
    A numeric scalar indicating the sample (i.e., column) of the assay to extract values to fill the ContactMatrix.
}

\item{fill}{
    A numeric vector of length equal to \code{nrow(x)}, containing values with which to fill the ContactMatrix.
    If specified, this overrides extraction of assay values for \code{inflate,InteractionSet-method}.
}

\item{swap}{
    A logical scalar indicating whether filling should also be performed after swapping anchor indices.
}

\item{sparse}{ 
    A logical scalar indicating whether the inflated matrix should use a sparseMatrix representation.
}

\item{collapse}{
    A logical scalar indicating whether duplicated interactions should be removed from \code{x} prior to deflation.
}

\item{extract}{
    A logical vector or matrix indicating which entries of \code{x} to convert into pairwise interactions.
}

\item{use.zero, use.na}{ 
    A logical scalar indicating whether to convert zero or \code{NA} entries into pairwise interactions.
}

\item{...}{
    For \code{inflate}, additional arguments to pass to \code{\link{overlapsAny}} when \code{rows} or \code{columns} is a GRanges.

    For \code{deflate}, additional arguments to pass to the InteractionSet constructor.
}
}

\value{
For \code{inflate}, a ContactMatrix is returned.

For \code{deflate}, an InteractionSet object is returned.
}

\section{Inflating to a ContactMatrix}{
The \code{inflate} method will return a ContactMatrix where the rows and columns correspond to specified regions of interest in \code{rows} and \code{columns}.
Regions can be specified by supplying an object of various types: 
\itemize{
\item{If it is an integer vector, it is assumed to refer to intervals in the \code{regions} slot of the input object \code{x}.
Values of the vector need not be sorted or unique, but must lie within \code{[1, regions(x)]}.}
\item{If it is a logical vector, it will subset to retain intervals in \code{regions(x)} that are \code{TRUE}.}
\item{If it is a character vector, it is assumed to contain the names of the reference sequences of interest (i.e., chromosome names).}
\item{If it is a GRanges object, \code{\link{overlapsAny}} will be called to identify the overlapping intervals of \code{regions(x)}.}
\item{If it is \code{NULL}, all regions in \code{regions(x)} will be used to construct that dimension of the ContactMatrix.}
}

For the GInteractions method, values in the matrix are filled based on user-supplied values in \code{fill}.
Each entry of \code{fill} corresponds to an interaction in \code{x} is used to set the matrix value at the matching row/column.
Some entries of the matrix will correspond to pairwise interactions that are not present in \code{x} - these are filled with \code{NA} values.

By default, filling is repeated after swapping the anchor indices.
This means that the value of the matrix at (1, 2) will be the same as that at (2, 1), i.e., the matrix is symmetric around the diagonal of the interaction space.
However, if \code{swap=FALSE}, filling is performed so that the first and second anchor indices correspond strictly to rows and columns, respectively.
This may be preferable if the order of the anchors contains some relevant information. 
In all cases, if duplicated interactions are present in \code{x} (and redundant permutations, when \code{swap=TRUE}), one will be arbitrarily chosen to fill the matrix.
 
For the InteractionSet \code{inflate} method, entries in the matrix are filled in based on the values in the first sample of the first assay by default.
For more complex \code{x}, values from different assays and samples can be extracted using the \code{assay} and \code{sample} arguments.
Note that if \code{fill} is specified, it will override any extraction of values from the assays.

If \code{sparse=TRUE}, \code{inflate} will return a ContactMatrix containing a sparseMatrix in the \code{matrix} slot.
Here, entries without a corresponding interaction in \code{x} are set to zero, not \code{NA}.
}

\section{Deflating from a ContactMatrix}{
The \code{deflate} method will return an InteractionSet where each relevant entry in the ContactMatrix is converted into a pairwise interaction.
Relevant entries are defined as those that are non-zero, if \code{use.zero} is \code{FALSE}; and non-\code{NA}, if \code{use.na} is \code{FALSE}.
If \code{x} contains a sparseMatrix representation, the former is set to \code{FALSE} while the latter is set to \code{TRUE}, if either are not specified.
For all other matrices, \code{use.zero=TRUE} and \code{use.na=FALSE} by default.

If \code{extract} is specified, this overrides all values of \code{use.zero} and \code{use.na}.
A typical application would be to \code{deflate} a number of ContactMatrix objects with the same \code{extract} matrix.
This ensures that the resulting InteractionSet objects can be easily combined with \code{cbind}, as the interactions are guaranteed to be the same.
Otherwise, different interactions may be extracted depending on the presence of zero or \code{NA} values.

The values of all matrix entries are stored as a one-sample assay, with each value corresponding to its pairwise interaction after conversion.
Duplicate interactions are removed by default, along with redundant permutations of the anchor indices.
These can be included in the returned object by setting \code{collapse=FALSE}.
This setting will also store the pairs as a GInteractions object, rather than using the default StrictGInteractions object where duplicates are not stored.

Additional arguments can be used to specify the \code{colData} and \code{metadata}, which are stored in the ContactMatrix itself.
}

\section{Interpreting zeroes in a sparse matrix}{
Storing data as a sparseMatrix may be helpful as it is more memory-efficient for sparse areas of the interaction space.
However, users should keep in mind that the zero values in the sparseMatrix may not represent zeroes in \code{fill}.
The majority of these values are likely to be zero just because there was no corresponding interaction in \code{x} to set it to a non-zero value.

Whether or not this is a problem depends on the application.
For example, if \code{fill} represents count data and only interactions with non-zero counts are stored in \code{x}, then setting all other entries to zero is sensible.
However, in other cases, it is not appropriate to fill entries corresponding to missing interactions with zero.
If \code{fill} represents, e.g., log-fold changes, then setting missing entries to a value of zero will be misleading.
One could simply ignore zeroes altogether, though this will also discard entries that are genuinely zero. 

These problems are largely avoided with the default denseMatrix classes, where missing entries are simply set to \code{NA}.
}

\examples{
example(InteractionSet, echo=FALSE)

inflate(iset, 1:10, 1:10)
inflate(iset, 1:10, 1:10, sparse=TRUE)
inflate(iset, 1:10, 1:5+10)
inflate(iset, "chrA", 1:5+10)
inflate(iset, "chrA", "chrB")
inflate(iset, "chrA", GRanges("chrB", IRanges(1, 10)))

y <- inflate(iset, 1:10, 1:10)
iset2 <- deflate(y)
iset2
assay(iset2)

y <- inflate(iset, 1:10, 1:10, swap=FALSE)
iset2 <- deflate(y)
iset2
assay(iset2)

# Testing with different fillings:
y <- inflate(iset, 1:10, 1:10, sample=2)
iset2 <- deflate(y)
assay(iset2)

y <- inflate(iset, 1:10, 1:10, fill=rowSums(assay(iset)))
iset2 <- deflate(y)
assay(iset2)

y2 <- inflate(interactions(iset), 1:10, 1:10, rowSums(assay(iset)))
identical(y, y2) # should be TRUE

# Effect of 'collapse'
y <- inflate(iset, c(8, 1:10), 1:10)
deflate(y)
deflate(y, collapse=FALSE)
}

\author{
Aaron Lun
}

\seealso{
\code{\link{InteractionSet-class}},
\code{\link{GInteractions-class}},
\code{\link{ContactMatrix-class}},
\code{\link{sparseMatrix-class}}
}
