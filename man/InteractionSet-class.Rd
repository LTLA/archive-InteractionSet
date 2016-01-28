\name{InteractionSet class}
\docType{class}

\alias{InteractionSet-class}
\alias{InteractionSet}
\alias{InteractionSet,ANY,GInteractions-method}
\alias{InteractionSet,missing,missing-method}

\title{InteractionSet class and constructors}

\description{
The InteractionSet class stores information about pairwise genomic interactions, and is intended for use in data analysis from Hi-C or ChIA-PET experiments.
Each row of the InteractionSet corresponds to a pairwise interaction between two loci, as defined in the GInteractions object.
Each column corresponds to a library or sample.
Each InteractionSet also contains one or more assays, intended to hold experimental data about interaction frequencies for each interaction in each sample.
}

\usage{
\S4method{InteractionSet}{ANY,GInteractions}(assays, interactions, ...)

\S4method{InteractionSet}{missing,missing}(assays, interactions, ...)
}

\arguments{
\item{assays}{
    A numeric matrix or a list or \code{SimpleList} of matrices, containing data for each interaction.
}

\item{interactions}{
    A GInteractions object of length equal to the number of rows in \code{assays}.
}

\item{...}{Other arguments to be passed to \code{\link{SummarizedExperiment}}.}
}

\value{
For the constructors, an InteractionSet object is returned.
}

\details{
The InteractionSet class inherits from the SummarizedExperiment class and has access to all of its data members and methods 
    (see \code{?\link{SummarizedExperiment-class}} for more details). 
It also contains an additional \code{interactions} slot which holds a GInteractions object (or an object from any derived classes, e.g., StrictGInteractions).
Each row of the InteractionSet object corresponds to a pairwise interaction between two genomic loci in \code{interactions}.

The constructor will return an InteractionSet object containing all of the specified information 
    - for \code{InteractionSet,missing,missing-method}, an empty InteractionSet object is returned.
Note that any \code{metadata} arguments will be placed in the metadata of the internal SummarizedExperiment object, 
    \emph{not} the metadata of the internal GInteractions object.
This is consistent with the behaviour of similar classes like RangedSummarizedExperiment.
}

\examples{
example(GInteractions, echo=FALSE)
Nlibs <- 4
counts <- matrix(rpois(Np*Nlibs, lambda=10), ncol=Nlibs)
colnames(counts) <- seq_len(Nlibs)

iset <- InteractionSet(counts, gi)
iset <- InteractionSet(counts, gi, colData=DataFrame(lib.size=1:Nlibs*1000))
iset <- InteractionSet(counts, gi, metadata=list(name="My Hi-C data"))

# Note differences in metadata storage:
metadata(iset)
metadata(interactions(iset))
}

\author{
Aaron Lun
}

\seealso{
\code{\link{SummarizedExperiment}},
\code{\link{interaction-access}},
\code{\link{interaction-subset}},
\code{\link{interaction-compare}},
\code{\link{SummarizedExperiment-class}}
}

