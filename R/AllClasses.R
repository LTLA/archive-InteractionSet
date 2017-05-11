# Defines the ContactMatrix class.

setClass("ContactMatrix",
    contains="Annotated", 
    slots=list(
        matrix="ANY", 
        anchor1="integer",
        anchor2="integer",
        regions="GRanges"
    )       
)

# Defines the GInteractions class, to hold interacting coordinates. 
# (Could inherit from Hits with an extra 'regions'; but I don't want to
# deal with 'queryHits' and 'subjectHits', which would get very confusing
# when you're dealing with queries and subjects in the overlap section.)


setClass("GInteractions", 
    contains="Vector",
    representation(
        anchor1="integer",
        anchor2="integer",
        regions="GRanges",
        NAMES="character_OR_NULL"
    )
)

setClass("StrictGInteractions", contains="GInteractions")
setClass("ReverseStrictGInteractions", contains="GInteractions")

# Defines the InteractionSet class, based on the SummarizedExperiment base class.
# This allows us to avoid re-defining various standard functions.

setClass("InteractionSet", 
    contains="SummarizedExperiment",
    representation(
        interactions="GInteractions"
    )
)

