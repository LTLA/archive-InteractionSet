\name{boundingBox}

\alias{boundingBox}
\alias{boundingBox,GInteractions-method}
\alias{boundingBox,InteractionSet-method}

\title{Get the bounding box}

\description{Computing a minimum bounding box for groups of pairwise interactions.}

\usage{
\S4method{boundingBox}{GInteractions}(x, f)

# Note, the same call is used for the InteractionSet method.
}

\arguments{
\item{x}{A GInteractions or InteractionSet object.}

\item{f}{
A factor or vector of length equal to that of \code{x}, indicating the group to which each pairwise interaction belongs.
}
}

\value{
A GInteractions object containing the coordinates of each bounding box.
}

\details{
For any group of pairwise interactions, the minimum bounding box is the smallest rectangle in the interaction space that contains all interactions in the group.
Each side of the box has coordinates spanning the most extreme anchor regions on the corresponding chromosome.
This is often useful for summarizing clusters of interactions.

Grouping of interactions is specified using \code{f}, where interactions in \code{x} with the same level of \code{f} are considered to be in the same group.
If \code{f} is not specified, all interactions in \code{x} are assumed to be in a single group (named as ``1'').
An error will be raised if a group spans multiple chromosomes for either the first or second anchor regions.

The function returns a GInteractions object containing the coordinates of the bounding boxes for all groups.
Each interaction represents a bounding box for a group, where the anchor regions represent the sides of the box.
Entries are named according to the levels of \code{f}, in order to specify which bounding box corresponds to which group.

It is recommended to run \code{\link{swapAnchors}} prior to computing the bounding box for intra-chromosomal groups.
If all \code{anchor1 >= anchor2} or all \code{anchor1 <= anchor2}, all interactions will lie on one side of the diagonal of the intra-chromosomal interaction space.
This results in the smallest possible minimum bounding box, which will only increase in size if interactions are placed on the other side of the diagonal.
Alternatively, users can specify a StrictGInteractions object as an input into \code{x}, in which \code{anchor1 <= anchor2} is enforced automatically.

% Bit of a pain to prove, but basically, if you flip a point to be above the diagonal, the Chebyshev distance to a point below the diagonal will always increase.
% This means that you must increase the size of one of your sides of your bounding box.
}

\seealso{
\code{\link{GInteractions-class}},
\code{\link{InteractionSet-class}}
}

\author{
Aaron Lun
}

\examples{
example(GInteractions, echo=FALSE)

# Making up a sensible grouping.
gi <- sort(gi)
all.chrs <- as.character(seqnames(regions(gi)))
f <- paste0(all.chrs[anchors(gi, type="first", id=TRUE)], ".",
            all.chrs[anchors(gi, type="second", id=TRUE)])

boundingBox(gi, f)
boundingBox(swapAnchors(gi), f)

# Fails for multiple chromosomes
try(out <- boundingBox(gi))
in.A <- f=="chrA.chrA"
out <- boundingBox(gi[in.A])
}


