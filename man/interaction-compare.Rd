\name{Interaction compare}
\alias{interaction-compare}

\alias{order,InteractionSet-method}
\alias{sort,InteractionSet-method}
\alias{duplicated,InteractionSet-method}
\alias{unique,InteractionSet-method}

\alias{order,GInteractions-method}
\alias{sort,GInteractions-method}
\alias{duplicated,GInteractions-method}
\alias{unique,GInteractions-method}

\alias{match,GInteractions,GInteractions-method}
\alias{match,InteractionSet,GInteractions-method}
\alias{match,GInteractions,InteractionSet-method}
\alias{match,InteractionSet,InteractionSet-method}

\alias{swapAnchors}
\alias{swapAnchors,GInteractions-method}
\alias{swapAnchors,InteractionSet-method}

\alias{pcompare,GInteractions,GInteractions-method}

\title{Interaction comparisons}

\description{
Methods to order, compare and de-duplicate GInteractions or InteractionSet objects, based on the anchor indices.
}

\usage{
\S4method{order}{GInteractions}(..., na.last=TRUE, decreasing=FALSE)

\S4method{sort}{GInteractions}(x, decreasing=FALSE, ...)

\S4method{duplicated}{GInteractions}(x, incomparables=FALSE, fromLast=FALSE, ...)

\S4method{unique}{GInteractions}(x, incomparables=FALSE, fromLast=FALSE, ...)

\S4method{swapAnchors}{GInteractions}(x, mode=c("order", "reverse", "all"))

## Each of the above methods has an identical equivalent for
## InteractionSet objects (not shown for brevity).

\S4method{match}{GInteractions,GInteractions}(x, table, nomatch=NA_integer_, incomparables=NULL, ...)

# The same call for 'match' is used for all possible combinations of
# GInteractions and InteractionSet objects in 'x' or 'table'. 
# For brevity, these will not be listed here.

\S4method{pcompare}{GInteractions,GInteractions}(x, y)
}

\arguments{
\item{...}{
    For \code{order}, one or more InteractionSet or GInteractions objects with the same number of rows.
    Otherwise, ignored in all other methods.
}

\item{x}{
    An InteractionSet or GInteractions object.
    For \code{pcompare}, a GInteractions object only.
}

\item{na.last}{
    A logical scalar indicating whether \code{NA} values should be ordered last.
    This should not be relevant as anchor indices should be finite.
}

\item{decreasing}{
    A logical scalar indicating whether rows should be sorted in decreasing order.
}

\item{incomparables}{
    A logical scalar.
    See \code{?\link[base]{match}} for a description of this in \code{match}.
    Otherwise, it is ignored.
}

\item{fromLast}{
    A logical scalar indicating whether the last entry of a repeated set in \code{x} should be considered as a non-duplicate.
}

\item{mode}{
    A string indicating what type of swapping should be performed in \code{swapAnchors}.
}

\item{table}{A GInteractions or InteractionSet object.}

\item{nomatch}{An integer scalar indicating the value to use for unmatched entries.}

\item{y}{A GInteractions object, of the same length as \code{x}.}
}

\value{
For \code{sort} and \code{unique}, a GInteractions or InteractionSet object is returned, depending on the class of \code{x}.

For \code{order} and \code{duplicated}, an integer vector of permutations or a logical vector of duplicate specifications is returned, respectively.
}

\section{Sorting and ordering}{
Sorting is performed based on the anchor indices of the GInteraction object.
Rows are ordered for increasing values of the \code{anchor1} slot - if these are equal, ordering is performed with values of the \code{anchor2} slot.
This equates to ordering by the coordinates directly, as the GRanges in the \code{regions} slot is always sorted.
Based on this, \code{sort} will return a sorted copy of \code{x} with permuted rows in increasing order.

The \code{order} method returns an integer vector indicating the permutation to rearrange \code{x} in increasing order.
If multiple objects are supplied to \code{order}, ordering will be computed using anchor indices from successive objects.
For example, ordering is first performed using anchor indices from the first object; 
    any rows with the same \code{anchor1} and \code{anchor2} will be split using the second object; and so on.

For both of these methods, the sorting can be reversed by setting \code{decreasing=TRUE}.
This will sort or order for decreasing values of \code{anchor1} and \code{anchor2}, rather than for increasing values.
}

\section{Removing duplicates}{
The \code{duplicated} method returns a logical vector indicating whether the rows of \code{x} are duplicated.
Duplicated rows are identified on the basis of identical entries in the \code{anchor1} and \code{anchor2} slots.
Values in other slots (e.g., in \code{mcols}) are ignored.
For a set of duplicated rows, the first occurrence in \code{x} is marked as the non-duplicate if \code{fromLast=FALSE}, and the last entry otherwise.

\code{unique} returns a GInteractions object where all duplicate rows have been removed from \code{x}.
This is equivalent to \code{x[!duplicated(x),]}, with any additional arguments passed to \code{duplicated} as specified.
}

\section{Matching and comparing}{
The \code{match} function will return an integer vector of length equal to that of \code{x}.
Each entry of the vector corresponds to an interaction in \code{x} and contains the index of the first interaction \code{table} with the same anchor regions.
Interactions in \code{x} without any matches are assigned \code{NA} values by default.
Note that the \code{regions} slot must be identical between \code{x} and \code{table}.

The \code{pcompare} function will return an integer vector of length equal to \code{x} and \code{y}.
This performs parallel comparisons between corresponding entries in two GInteractions objects, based on the values of the anchor indices 
(\code{anchor1} first, and then \code{anchor2} if \code{anchor1} is tied).
Negative, zero and positive values indicate that the corresponding interaction in \code{x} is 'lesser', equal or 'greater' than the corresponding interaction in \code{y}.
Again, the \code{regions} must be identical between the two objects.
}

\section{Swapping anchors}{
For GInteractions objects, \code{swapAnchors} returns a GInteractions object where \code{anchor1} and \code{anchor2} values are swapped.
If \code{mode="order"}, this is done so that all values in the \code{anchor2} slot are not less than values in \code{anchor1} in the returned object.
If \code{mode="reverse"}, all values in \code{anchor1} are not less than all values in \code{anchor2}.
If \code{mode="all"}, the anchor indices in \code{x} are directly swapped without consideration of ordering.

It is recommended to apply this method before running methods like \code{order} or \code{duplicated}.
This ensures that redundant permutations are coerced into the same format for a valid comparison.
In many applications, permutations of pairwise interactions are not of interest, i.e., an interaction between regions 1 and 2 is the same as that between 2 and 1.
Application of \code{swapAnchors} with \code{mode="order"} ensures that all indices are arranged in a comparable manner.
Alternatively, users can use a (Reverse)StrictGInteractions object which enforces a consistent arrangement of indices across interactions.
}

\section{Methods for InteractionSet objects}{
For all InteractionSet methods, the corresponding method is called on the GInteractions object in the \code{interactions} slot of the InteractionSet object.
Return values for each InteractionSet method is the same as those for the corresponding GInteractions method 
    - except for \code{sort} and \code{unique}, which return a row-permuted or row-subsetted InteractionSet, respectively, instead of a GInteractions object;
      and \code{swapAnchors}, which returns an InteractionSet object where the internal GInteractions has its anchor indices swapped around.      

Note that no additional information from the InteractionSet (beyond that in \code{interactions}) is used for sorting or duplicate marking,
    i.e., the assay or metadata values for each interaction are \emph{not} used in distinguishing rows with the same interaction.
For this reason, the \code{pcompare} method is not implemented for InteractionSet objects.
It makes little sense to do a parallel comparison in an InteractionSet without examining the data.
}

\examples{
example(GInteractions, echo=FALSE)

anchors(gi, id=TRUE)
anchors(swapAnchors(gi, mode="all"), id=TRUE)
gi <- swapAnchors(gi)

anchors(gi)
gi2 <- sort(gi)
gi2
anchors(gi2)

# Can also order them.
o <- order(gi)
o
stopifnot(all.equal(gi[o], gi2))

# Checking duplication.
summary(duplicated(gi))
summary(duplicated(rbind(gi, gi)))
stopifnot(all.equal(gi, unique(rbind(gi, gi))))

# Matching and comparing.
another.gi <- gi[sample(length(gi))]
match(gi, another.gi)
match(gi, another.gi[1:5])

pcompare(gi, another.gi)

#################
# Same can be done for an InteractionSet object:

example(InteractionSet, echo=FALSE)
iset <- swapAnchors(iset)

anchors(iset)
anchors(sort(iset))
order(iset)
summary(duplicated(iset))
unique(iset)
}

\author{
Aaron Lun
}

\seealso{
\code{\link{GInteractions-class}},
\code{\link{InteractionSet-class}},
\code{\link[base]{match}},
\code{\link[S4Vectors]{pcompare}}
}

