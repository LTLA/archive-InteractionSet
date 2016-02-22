# Defines the match function.

setMethod("match", c("GInteractions", "GInteractions"), 
    function(x, table, nomatch = NA_integer_, incomparables = NULL, ...) {

        .strict_check(x, table)
        if (!identical(regions(x), regions(table))) { 
            stop("'regions' must be identical for arguments to 'match'")
        }

        # Using the Hits method, for convenience (this automatically resorts).
        nR <- length(regions(x))
        a1 <- Hits(anchors(x, id=TRUE, type="first"),
                   anchors(x, id=TRUE, type="second"), nR, nR, order=seq_along(x))
        a2 <- Hits(anchors(table, id=TRUE, type="first"),
                   anchors(table, id=TRUE, type="second"), nR, nR, order=seq_along(table))
        out <- match(a1, a2, nomatch=nomatch, incomparables=incomparables, ...)

        # Unscrambling:
        out[mcols(a1)$order] <- mcols(a2)$order[out]
        return(out)
    }
)

setMethod("match", c("GInteractions", "InteractionSet"), 
    function(x, table, nomatch = NA_integer_, incomparables = NULL, ...) {
        match(x, interactions(table), nomatch=nomatch, incomparables=incomparables, ...)
    }
)

setMethod("match", c("InteractionSet", "GInteractions"), 
    function(x, table, nomatch = NA_integer_, incomparables = NULL, ...) {
        match(interactions(x), table, nomatch=nomatch, incomparables=incomparables, ...)
    }
)

setMethod("match", c("InteractionSet", "InteractionSet"), 
    function(x, table, nomatch = NA_integer_, incomparables = NULL, ...) {
        match(interactions(x), table, nomatch=nomatch, incomparables=incomparables, ...)
    }
)

# Defining the compare function (only for two GInteractions, it wouldn't make sense otherwise).

setMethod("pcompare",  c("GInteractions", "GInteractions"), function(x, y) { 
    .strict_check(x, y)
    if (length(regions(x))!=length(regions(y)) || any(regions(x)!=regions(y))) { 
        stop("'regions' must be identical for arguments to 'pcompare'")
    }
    output <- x@anchor1 - y@anchor1
    tied1 <- output==0L
    output[tied1] <- x@anchor2[tied1] - y@anchor2[tied1]
    return(output)
})

.strict_check <- function(x, y) {
    if (is(x, "StrictGInteractions")!=is(y, "StrictGInteractions") 
        || is(x, "ReverseStrictGInteractions")!=is(y, "ReverseStrictGInteractions")) { 
        warning("comparison between GInteractions objects of different strictness")
    }
}
