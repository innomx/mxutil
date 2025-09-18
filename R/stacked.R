.myseqalong <- match.fun("seq_along")
.mymatchfun <- match.fun


#' @export
difftime_default <- function(time1, time2) {
    if (length(time1) != length(time2)) {
        stop("Length of time1 and time2 should match.")
    }
    if (inherits(time1, "POSIXt") && inherits(time2, "POSIXt")) {
        as.double(difftime(time1, time2), units="hours")
    } else if (inherits(time1, "Date") && inherits(time2, "Date")) {
        as.numeric(time1 - time2) * 24
    } else if (is.numeric(time1) && is.numeric(time2)) {
        time1 - time2
    } else {
        stop("Cannot subtract times that are not either numeric or of class POSIXt or Date.")
    }
}


#' @export
difftime_hours <- function(time1, time2) {
    if (length(time1) != length(time2)) {
        stop("Length of time1 and time2 should match.")
    }
    if (!(inherits(time1, "POSIXt") && inherits(time2, "POSIXt"))) {
        stop("Both time1 and time2 should be POSIXt objects")
    }
    as.double(difftime(time1, time2), units="hours")
}


#' @export
difftime_days <- function(time1, time2) {
    if (length(time1) != length(time2)) {
        stop("Length of time1 and time2 should match.")
    }
    if (!((inherits(time1, "POSIXt") || (inherits(time1, "Date"))) && (inherits(time2, "POSIXt")) || inherits(time2, "Date"))) {
        stop("Both time1 and time2 should be Date or POSIXt objects")
    }
    if (inherits(time1, "POSIXt")) {
        time1 <- as.Date(strftime(time1, "%Y-%m-%d"))  # Drop the time
    }
    if (inherits(time2, "POSIXt")) {
        time2 <- as.Date(strftime(time2, "%Y-%m-%d"))  # Drop the time
    }
    as.integer(as.double(difftime(time1, time2), units="days"))
}


# ordering to restore the original order of the data after unlist(tapply(...))
#.restoreOrder <- function(id) {
#    order(unlist(tapply(.myseqalong(id), id, I, simplify=FALSE)))
#}

#' Restore Order
#'
#' @param o a permunation of indices from 1 to n (where n is `length(o)`)
#' @return a permunation of indices from 1 to n
#' @examples
#' x <- c(2, 3, 9, 1, 8, 4, 6, 7, 5)
#' x <- x[o <- order(x)] # Now x is sorted
#' x[restoreOrder(o)]    # But we can restore the original order
#' @export
restoreOrder <- function(o) {
    invord <- numeric(length(o))
    invord[o] <- 1:length(o)
    invord
}

#' Make a Valid Block Format ID From One or More Vectors
#'
#' A [stacked] ID is constructed from one or more atomic vectors. If more than
#' one atomic vector is supplied, then a compound ID is constructed by
#' concatenating the values together with a separator.
#' @param id An atomic vector.
#' @param ... More atomic vectors.
#' @param sep A character to use as a separator when forming a compound ID.
#' @return A [stacked] ID, i.e. a factor `x` such that [is.unsorted]`(x)` is `FALSE`.
#' @author Benjamin Rich <mail@benjaminrich.net>
#' @seealso [stacked]
#' @examples
#' require(nlme)
#' data(Phenobarb)
#' dat <- Phenobarb[1:56,]  # First 4 subjects
#' 
#' attach(dat)
#' 
#' is.unsorted(Subject)           # TRUE  - the original column is NOT sorted
#' is.unsorted(asID(Subject))     # FALSE - the ID is sorted
#' 
#' all(as.character(Subject) == as.character(asID(Subject)))  # TRUE - textually they are the same
#' 
#' cbind(ID=asID(Subject), dat)
#' 
#' detach(dat)
#' @export
asID <- function(id, ..., sep="|") {
    args <- list(...)
    if (length(args) > 0) {
        if (!all(sapply(args, length) == length(id))) {
            stop("Elements of a compound ID must have the same length")
        }
        if (any(is.na(id)) || !all(sapply(args, nmissing) == 0)) {
            warning("One or more arguments contains missing values. 'NA' will be treated as a distict value when forming the compound ID.")
        }
        id2 <- do.call(function(...) paste(..., sep=sep), args)
        id <- paste(id, id2, sep=sep)
    } else {
        if (any(is.na(id))) {
            warning("Argument contains missing values.")
        }
    }
    id <- factor(id, levels=unique(id))
    if (is.unsorted(id)) {
        # Don't issue a warning here.  The caller's intension may be to sort later.
        #warning("Not a valid block-format ID")
    }
    id
}


noID <- function(x) {
    asID(rep("", length(x)))
}


# How many missing values?
nmissing <- function(x) {
    if (is.null(x)) { warning("x is NULL", call.=F) }
    sum(is.na(x))
}


# How many unique values?
nunique <- function(x, na.rm=FALSE) {
    if (na.rm) {
        length(unique(x[!is.na(x)]))
    } else {
        length(unique(x))
    }
}


#' @export
is_sorted <- function(...) {
    args <- list(...)
    if ((length(args) == 1) && is.atomic(args[[1]])) {
        # Handle the most common cases (numeric and factor) more efficiently
        !is.unsorted(args[[1]])
    } else {
        !is.unsorted(order(...))
    }
}


.checkID <- function(id) {
    if (!is.factor(id)) {
        stop("Not a valid ID: A block-format ID must be a factor.")
    } else if (any(is.na(id))) {
        stop("Not a valid ID: A block-format ID must not contain missing values.")
    } else if (is.unsorted(unique(id))) {
        stop("Not a valid ID: A block-format ID must be sorted.  Try using asID(...) to create a valid ID.")
    } else if (is.unsorted(id)) {
        first.prob1 <- which(base::diff(order(id)) > 1)[1]
        first.prob2 <- which(base::diff(as.integer(id)) < 0)[1] + 1
        id2 <- gsub("\r", "|", id)
        temp <- c(        "    position    id",
                  sprintf("        %4s    %-s", first.prob1, id[first.prob1]),
                          "         ...    ...",
                  sprintf("        %4s    %-s", first.prob2-1, id[first.prob2-1]),
                  sprintf("   ---> %4s    %-s", first.prob2, id[first.prob2]))

        temp2 <- paste(temp, collapse="\n")
        stop("Not a valid ID: A block-format ID must consist of contiguous blocks.\n",
             "There can only be one block per ID.\n",
             "The first recurring ID was identified here:\n\n\n",
             temp2,
             "\n\n")
    }
}


#' Extract a Single Value for Each Block in a Block Data Set
#'
#' Extract a single value for each block in a [stacked] data set.
#'
#' @param x A vector in \code{\link{block-format}} with respect to `id`.
#' @param id A valid \code{\link{block-format}} ID for the input `x`.
#' @param id2 A valid \code{\link{block-format}} ID for the return value.
#' @param ind A vector of logicals used to filter the values of `x`.
#' @param select When more than one value exists, which one to select.  See details.
#' @param fill A value to use when no other value is appropriate. See details.
#'
#' @details
#' These functions allow one to extract a single value per block in a
#' \code{\link{block-format}} data set.  This can be useful in many
#' contexts.  If the values of \code{x} in each block are not unique,
#' then a specific value needs to be determined.  The indicator vector
#' \code{ind} can be used to filter out specific rows that contain
#' values of interest.  If there is more than one value, then
#' \code{select} is used to choose either the first or the last
#' (according to the ordering of \code{x}).  If, for a given block,
#' there is no value at all, either because none of the rows matched
#' the \code{ind} criteria or because the block is of size zero (see
#' \code{\link{block-format}}), then the value of \code{fill} is used
#' for that block.
#' 
#' Both \code{blk_singleValue} and \code{blk_repeatValue} determine a
#' unique value for each block.  The difference between them is how
#' many times each value is repeated in the result.  The vector
#' returned by \code{blk_singleValue} has length equal to the number of
#' blocks in \code{id}, i.e. each value appears exactly once.  For
#' \code{blk_repeatValue}, each value is repeated the appropriate
#' number of times so that the result is a vector in block format with
#' respect to \code{id2}.  Thus, \code{blk_repeatValue} can effectively
#' be used to perform a simple ``left outer join'' (or ``merge'')
#' operation on a single variable (see example).
#' 
#' Neither of \code{levels(id)} and \code{levels(id2)} need be a proper
#' subset of the other.  For levels of \code{id2} (resp. \code{id}) that
#' are not levels of \code{id} (resp. \code{id2}), the corresponding blocks
#' are assumed to be of size zero in \code{id} (resp. \code{id2}).
#' 
#' It is also possible to use \code{blk_repeatValue} in a different
#' way, without specifying \code{id}.  In that case, \code{x} must
#' have length equal to \code{nlevels(id2)}, i.e. it must contain a
#' unique value for each possible value of \code{id2}, and the
#' correspondence of values to ID's is taken from their respective
#' ordering.  Then, each value in \code{x} is repeated the number of
#' times that the corresponding ID appears in \code{id2}.
#' This is equivalent to \code{blk_repeatValue(x, id=asID(levels(id2)), id2, ...)}.
#'
#' @return
#' `blk_singleValue` returns a vector containing one value for each level of
#' `id`.
#' 
#' `blk_repeatValue` returns a vector in `\link{block-format}` with respect to
#' `id2`.
#'
#' @author Benjamin Rich <mail@benjaminrich.net>
#' 
#' @seealso
#' - [block-format]
#' - [tapply]
#' - [merge]
#' - [duplicated]
#'
#' @examples
#' 
#' # EXAMPLE 1
#' require(nlme)
#' data(Phenobarb)
#' dat <- Phenobarb[1:56,]  # First 4 subjects
#' dat$id <- asID(dat$Subject)
#' 
#' attach(dat)
#' 
#' # A single row per subject
#' data.frame(id=levels(id),
#'            Wt=blk_singleValue(Wt, id),
#'            Apgar=blk_singleValue(Apgar, id),
#'            final.dose=blk_singleValue(dose, id, ind=(!is.na(dose)), select="last"))
#' 
#' # Repeat a single value on each row for each subject
#' cbind(dat, data.frame(
#'   first.dose=blk_repeatValue(dose, id, ind=(!is.na(dose))),
#'   final.dose=blk_repeatValue(dose, id, ind=(!is.na(dose)), select="last")
#' ))
#' 
#' detach(dat)
#' 
#' ### Merging a time-fixed covariate (simple left outer join)
#' ### -------------------------------------------------------
#' 
#' # Suppose subjects 1 and 2 are Male, and Subject 4 is Female, but the
#' # gender of subject 3 is not specified.
#' gender <- data.frame(
#'     id=factor(c(1, 2, 4), levels=levels(dat$id)),   # Note: keeping the same factor levels helps
#'     gender=c("Male", "Male", "Female"))
#' 
#' gender
#' 
#' # Now, 'merge' the gender with the rest of the data.
#' # Since subject 3 is absent, it gets the value of 'fill', i.e. NA.
#' dat$gender <- blk_repeatValue(gender$gender, gender$id, dat$id)
#' dat
#' 
#' # Still returns 4 values:
#' blk_singleValue(gender$gender, gender$id)
#' 
#' ### The other way of using blk_repeatValue (without specifying id)
#' ### --------------------------------------------------------------
#' 
#' letter <- LETTERS[1:nlevels(dat$id)]  # Exactly one value per id
#' cbind(dat,    letter=blk_repeatValue(c("A", "B", "C", "D"), id2=dat$id))
#' cbind(gender, letter=blk_repeatValue(c("A", "B", "C", "D"), id2=gender$id))
#' 
#' # EXAMPLE 2
#' id <- gl(4, 4)
#' x <- LETTERS[1:16]
#' y <- Sys.time() + 1:16
#' 
#' data.frame(
#'     id       = levels(id),
#'     first.x  = blk_singleValue(x, id),
#'     first.y  = blk_singleValue(y, id))
#' 
#' data.frame(
#'     id       = id,
#'     x        = x,
#'     first.x  = blk_repeatValue(x, id),
#'     y        = y,
#'     first.y  = blk_repeatValue(y, id))
#' 
#' target.id <- gl(4, 6)
#' data.frame(
#'     id       = target.id,
#'     first.x  = blk_repeatValue(x, id, target.id),
#'     first.y  = blk_repeatValue(y, id, target.id))
#' 
#' @export
blk_singleValue <- function(x, id=noID(x), ind=NULL, select=c("first", "last"), fill=NA) {
    .checkID(id)
    if (length(x) == 1) {
        x <- rep(x, length(id))
    }
    if (length(x) != length(id)) {
        stop("Lengths of x and id must match.")
    }
    if (is.null(ind)) {
        ind <-rep(TRUE, length(id))
    }
    if (!is.logical(ind)) {
        stop("Expecting a logical ind.")
    }
    if (length(ind) != length(id)) {
        stop("Lengths of ind and id must match.")
    }

    select <- match.arg(select)

    if (is.character(select) && select == "first") {
        apply.fun <- function(z, ...) z[1]
    } else if (is.character(select) && select == "last") {
        apply.fun <- function(z, ...) z[length(z)]
    } else {
        stop(sprintf("Unrecognized value for 'select' argument: %s", select))
    }


    yl <- tapply(.myseqalong(id)[ind], id[ind], apply.fun, simplify=FALSE)

    use.fill <- sapply(yl, function(z) (is.null(z) || length(z) == 0 ))

    yl <- lapply(yl, function(z) {if (is.null(z) || length(z) == 0 ) NA else z})

    if (any(sapply(yl, length) != 1)) {
        # Unreachable code
        warning("Something is wrong.  Only one value per ID expected.  The first value will be used.")
        yl <- lapply(yl, function(z) z[1])
    }

    y <- x[unlist(yl)]
    y[use.fill] <- fill
    y
}
#' @export
singleValue <- blk_singleValue


#' @rdname blk_singleValue
#' @export
blk_repeatValue <- function(x, id=noID(x), id2=id, ind=NULL, select=c("first", "last"), fill=NA) {
    .checkID(id2)
    if (missing(id)) {
        if (length(x) == nlevels(id2)) {
            id <- factor(levels(id2), levels=levels(id2))
        } else {
            stop("Lengths of x and id must match. When id is missing, its default value is asID(levels(id2)).")
        }
    }

    x1 <- blk_singleValue(x, id, ind, select, fill=fill)

    x2 <- x1[factor(id2, levels=levels(id))]
    x2[is.na(factor(id2, levels=levels(id)))] <- fill
    x2
}
#' @export
repeatValue <- blk_repeatValue


#' @export
blk_applyNto1 <- function(x, id=noID(x), apply.fun, ind=NULL, fill=NA, ...) {
    .checkID(id)
    if (length(x) != length(id)) {
        stop("Lengths of x and id must match.")
    }
    if (is.null(ind)) {
        ind <-rep(TRUE, length(id))
    }
    if (!is.logical(ind)) {
        stop("Expecting a logical ind.")
    }
    if (length(ind) != length(id)) {
        stop("Lengths of ind and id must match.")
    }

    if (is.character(apply.fun) && apply.fun == "first") {
        apply.fun <- function(z, ...) z[1]
    } else if (is.character(apply.fun) && apply.fun == "last") {
        apply.fun <- function(z, ...) z[length(z)]
    } else {
        apply.fun <- .mymatchfun(apply.fun)
    }

    yl <- tapply(x[ind], id[ind],
        function(z, apply.fun, ...) {
            res <- apply.fun(z, ...)
            if (length(res) != 1) {
                warning(sprintf("Expected 1 value, got %s.", length(res)))
                res <- rep(res, length.out=1)
            }
            res
        },
        simplify=FALSE, apply.fun=apply.fun, ...)

    yl <- lapply(yl, function(z) {if (is.null(z) || length(z) == 0 ) fill else z})

    if (any(sapply(yl, length) != 1)) {
        # Unreachable code
        warning("Something is wrong.  Only one value per ID expected.  The first value will be used.")
        yl <- lapply(yl, function(z) z[1])
    }

    y <- do.call(c, yl)
    y
}
#' @export
applyNto1 <- blk_applyNto1


#' @export
blk_applyNtoN <- function(x, id=noID(x), apply.fun, ind=NULL, fill=NA, ...) {
    .checkID(id)
    if (length(x) != length(id)) {
        stop("Lengths of x and id must match.")
    }
    if (is.null(ind)) {
        ind <-rep(TRUE, length(id))
    }
    if (!is.logical(ind)) {
        stop("Expecting a logical ind.")
    }
    if (length(ind) != length(id)) {
        stop("Lengths of ind and id must match.")
    }

    if (is.character(apply.fun) && apply.fun == "first") {
        apply.fun <- function(z, ...) z[1]
    } else if (is.character(apply.fun) && apply.fun == "last") {
        apply.fun <- function(z, ...) z[length(z)]
    } else {
        apply.fun <- .mymatchfun(apply.fun)
    }

    yl <- tapply(x[ind], id[ind],
        function(z, apply.fun, ...) {
            res <- apply.fun(z, ...)
            if (length(z) %% length(res) != 0) {
                warning(sprintf("Expected length of vector to be a factor of %s, but it was %s.", length(z), length(res)))
            }
            rep(res, length.out=length(z))
        },
        simplify=FALSE, apply.fun=apply.fun, ...)

    y <- do.call(c, yl)
    z <- rep(y[1], length(x))
    z[ind] <- y
    z[!ind] <- fill
    z
}
#' @export
applyNtoN <- blk_applyNtoN


#' @export
blk_firstOnly <- function(id=noID(ind), ind=NULL) {
    .checkID(id)
    if (is.null(ind)) {
        ind <-rep(TRUE, length(id))
    }
    if (!is.logical(ind)) {
        warning("Expecting a logical ind")
        ind <- as.logical(ind)
    }
    if (length(ind) != length(id)) {
        stop("Lengths of ind and id must match")
    }
    (ind * unlist(tapply(ind, id, cumsum, simplify=FALSE))) == 1
}
#' @export
firstOnly <-  blk_firstOnly
#' @export
blk_onlyFirst <- blk_firstOnly
#' @export
onlyFirst <- blk_onlyFirst


#' @export
blk_lastOnly <- function(id=noID(ind), ind=NULL) {
    .checkID(id)
    if (is.null(ind)) {
        ind <-rep(TRUE, length(id))
    }
    if (!is.logical(ind)) {
        warning("Expecting a logical ind")
        ind <- as.logical(ind)
    }
    if (length(ind) != length(id)) {
        stop("Lengths of ind and id must match")
    }
    unlist(tapply(ind, id, function(x) {
        if (!any(x)) x
        else {
            y <- cumsum(x)
            (x * y) == max(y)
        }
    }, simplify=FALSE))
}
#' @export
lastOnly <-  blk_lastOnly
#' @export
blk_onlyLast <- blk_lastOnly
#' @export
onlyLast <- blk_onlyLast


#' @export
blk_firstOnwards <- function(id=noID(ind), ind, include.first=TRUE, fill=FALSE) {
    .checkID(id)
    if (!is.logical(ind)) {
        warning("Expecting a logical ind")
        ind <- as.logical(ind)
    }
    if (length(ind) != length(id)) {
        stop("Lengths of ind and id must match")
    }
    unlist(tapply(ind, id, function(x, include.first, fill) {
        if (!any(x)) rep(fill, length(x))
        else {
            y <- cumsum(x)
            if (include.first) {
                y > 0
            } else {
                (y > 0) & (x * y != 1)
            }
        }
    }, include.first=include.first, fill=fill, simplify=FALSE))
}
#' @export
firstOnwards <- blk_firstOnwards


#' @export
blk_lastOnwards <- function(id=noID(ind), ind, include.last=TRUE, fill=FALSE) {
    .checkID(id)
    if (!is.logical(ind)) {
        warning("Expecting a logical ind")
        ind <- as.logical(ind)
    }
    if (length(ind) != length(id)) {
        stop("Lengths of ind and id must match")
    }
    unlist(tapply(ind, id, function(x, include.last, fill) {
        if (!any(x)) rep(fill, length(x))
        else {
            y <- cumsum(x)
            if (include.last) {
                y == max(y)
            } else {
                (y == max(y)) & (!x)
            }
        }
    }, include.last=include.last, fill=fill, simplify=FALSE))
}
#' @export
lastOnwards <- blk_lastOnwards


#' @export
blk_untilFirst <- function(id=noID(ind), ind, include.first=TRUE, fill=FALSE) {
    !blk_firstOnwards(id, ind, !include.first, !fill)
}
untilFirst <- blk_untilFirst


#' @export
blk_untilLast <- function(id=noID(ind), ind, include.last=TRUE, fill=FALSE) {
    !blk_lastOnwards(id, ind, !include.last, !fill)
}
#' @export
untilLast <- blk_untilLast


#' @export
blk_locf <- function(x, id=noID(x), na.action=c("fill", "carry.back"), fill=NA) {
    .checkID(id)
    if (length(x) != length(id)) {
        stop("Lengths of x and id must match")
    }
    
    na.action <- match.arg(na.action)

    ii <- ifelse(is.na(x), NA, .myseqalong(x))
    ii <- unlist(tapply(ii, id, function(y, na.action) {
        if (all(is.na(y))) return (y)
        z <- cumsum(!is.na(y))
        if (na.action == "carry.back") {
            z[z==0] <- 1
        }
        ifelse(z==0, NA, y[!is.na(y)][ifelse(z>0, z, 1)])
    }, na.action=na.action, simplify=FALSE))
    #ifelse(is.na(ii), fill, x[ii])  # Note: doesn't work with class POSIXt
    y <- x
    y[!is.na(ii)] <- x[ii[!is.na(ii)]]
    y[is.na(ii)] <- fill
    y
}
#' @export
locf <- blk_locf


#' @export
blk_nocb <- function(x, id=noID(x), na.action=c("fill", "carry.forward"), fill=NA) {
    .checkID(id)
    if (length(x) != length(id)) {
        stop("Lengths of x and id must match")
    }
    
    na.action <- match.arg(na.action)
    if (na.action == "carry.forward") na.action <- "carry.back"

    rid <- factor(rev(id), levels=rev(levels(id)))
    rev(blk_locf(rev(x), rid, na.action, fill))
}
#' @export
nocb <- blk_nocb


#' @export
blk_sequentialID <- function(id) {
    cumsum(blk_firstOnly(id=id))
}
#' @export
sequentialID <- blk_sequentialID


#' @export
blk_isConstant <- function(x, id=noID(x), ind=NULL) {
    all(blk_applyNto1(x, id, nunique, ind=ind, fill=0) <= 1)
}
#' @export
isConstant <- blk_isConstant

#' @export
blk_isVarying <- function(x, id=noID(x), ind=NULL) {
    any(blk_applyNto1(x, id, nunique, ind=ind, fill=0) > 1)
}
#' @export
isVarying <- blk_isVarying


#' @export
blk_count <- function(id=noID(ind), id2=id, ind=NULL) {
    res <- blk_applyNto1(x=rep(1, length(id)), id=id, apply.fun=sum, ind=ind, fill=0)
    if (!is.null(id2)) {
        res <- blk_repeatValue(res, id=asID(levels(id)), id2=id2, fill=0)
    }
    res
}
#' @export
stacked_count <- blk_count


#' @export
blk_changeFromBaseline <- function(x, id=noID(x), baseline=blk_singleValue(x, id, ind=ind, select=select), ind=NULL, select="first", diff.op="-", ...) {
    baseline <- blk_repeatValue(baseline, id2=id)
    .mymatchfun(diff.op)(x, baseline, ...)
}
#' @export
changeFromBaseline <- blk_changeFromBaseline


#' @export
blk_relativeTime <- function(time, id=noID(time), baseline=blk_singleValue(time, id, ind=ind, select=select), ind=NULL, select="first", diff.op=difftime_default, ...) {
    blk_changeFromBaseline(time, id, baseline, diff.op=diff.op, ...)
}
#' @export
relativeTime <- blk_relativeTime


#' @export
blk_shift <- function(x, id=noID(x), shift.by=1, ind=NULL, fill=NA) {
    .checkID(id)
    if (length(x) != length(id)) {
        stop("Lengths of x and id must match.")
    }
    if (is.null(ind)) {
        ind <-rep(TRUE, length(id))
    }

    if (shift.by == 0) return (x)

    ii <- rep(NA, length(x))
    ii[ind] <- unlist(tapply(.myseqalong(x)[ind], id[ind], function(z, shift.by) {
        if (length(z) <= abs(shift.by)) {
            rep(NA, length(z))
        } else if (shift.by > 0) {
            c(rep(NA, length=shift.by), z[1:(length(z) - shift.by)])
        } else {
            shift.by <- -shift.by
            c(z[(shift.by + 1):length(z)], rep(NA, length=shift.by))
        }
    }, shift.by=shift.by, simplify=FALSE))
    y <- x[NA]
    y[!is.na(ii)] <- x[ii[!is.na(ii)]]
    y[is.na(ii)] <- fill
    #y <- ifelse(is.na(ii), fill, x[ii])
    y

    #y <- rep(fill, length(x))
    #y[ind] <- unlist(tapply(x[ind], id[ind], function(z, shift.by, fill) {
    #    if (length(z) <= shift.by) {
    #        rep(fill, length(z))
    #    } else if (shift.by > 0) {
    #        c(rep(fill, length=shift.by), z[1:(length(z) - shift.by)])
    #    } else {
    #        shift.by <- -shift.by
    #        c(z[(shift.by + 1):length(z)], rep(fill, length=shift.by))
    #    }
    #}, shift.by=shift.by, fill=fill, simplify=FALSE))
    #y
}
#' @export
stacked_shift <- blk_shift


# Re-implement the standard diff function to allow for a diff.op function
.mydiff <- function(x, lag=1, diff.op="-", ...) {
    if (length(lag) > 1L || lag < 1L)  {
        stop("'lag' must be an integer >= 1")
    }

    diff.op <- .mymatchfun(diff.op)
 
    if (lag >= length(x)) {
        return (diff.op(x[1], x[1])[0])  # Zero length
    }

    r <- x
    i1 <- -seq_len(lag)
    r <- diff.op(r[i1], r[-length(r):-(length(r) - lag + 1L)], ...)
    r
}


#' @export
blk_diff <- function(x, id=noID(x), ind=NULL, lag=1, fill=NA, diff.op="-", ...) {
    .checkID(id)

    if (is.null(ind)) {
        ind <-rep(TRUE, length(id))
    }

    y <- rep(fill, length(x))
    y[ind] <- unlist(tapply(x[ind], id[ind], function(y, lag, fill, ...) {
        if (length(y) <= lag) {
            rep(fill, length(y))
        } else {
            c(rep(fill, lag), .mydiff(y, lag=lag, diff.op=diff.op, ...))
        }
    }, lag=lag, fill=fill, simplify=FALSE, ...))
    y
}
#' @export
stacked_diff <- blk_diff

# Alias, with different default for diff.op
#' @export
blk_intereventTime <- function(time, id=noID(time), ind=NULL, lag=1, fill=NA, diff.op=difftime_default, ...) {
    blk_diff(x=time, id=id, ind=ind, lag=lag, fill=fill, diff.op=diff.op, ...)
}

#' @export
intereventTime <- blk_intereventTime


.consecLengths <- function(ind) {
    y <- table(cumsum(!ind)) - 1
    if (ind[1]) {
        y[1] <- y[1] + 1
    }
    y <- y[y>0]
    y
}


# May be useful as well for non-block data?
.findConsecutive <- function(ind, min.consec=2) {
    y <- .consecLengths(ind)
    z <- rep(0, length(ind))
    z[ind] <- rep(ifelse(y >= min.consec, 1, 0), times=y)
    y <- y[y >= min.consec]
    z[z==1] <- rep(.myseqalong(y), times=y)
    z
}


#' @export
blk_findConsecutive <- function(id=noID(ind), ind, min.consec=2, fill=NA) {
    .checkID(id)

    y <- unlist(tapply(ind, id, .consecLengths, simplify=FALSE))
    z <- rep(0, length(ind))
    z[ind] <- rep(ifelse(y >= min.consec, 1, 0), times=y)
    y <- y[y >= min.consec]
    z[z==1] <- rep(.myseqalong(y), times=y)
    #z <- z - blk_applyNtoN(z-1, id, min, ind=(z > 0), fill=0)
    z[z == 0] <- fill
    z
}
#' @export
findConsecutive <- blk_findConsecutive


#' @export
blk_concomitant <- function(time, id=noID(time), start.time, end.time, id2=id, ends=c("[]", "[)", "(]", "()")) {
    .checkID(id)
    .checkID(id2)

    ends <- match.arg(end)

    o <- switch(ends,
        "[]" = c(1, 0, 2),
        "[)" = c(2, 0, 1),
        "(]" = c(0, 1, 2),
        "()" = c(1, 2, 0))

    temp <- rbind(
        data.frame(id=id,  time=time,       code=2, x=0),
        data.frame(id=id2, time=start.time, code=0, x=1),
        data.frame(id=id2, time=end.time,   code=1, x=-1))

    temp <- temp[order(temp$id, temp$time, temp$code),]
    y    <- unlist(tapply(temp$x, temp$id, cumsum, simplify=FALSE))
    y[temp$code==2]
}
#' @export
concomitant <- blk_concomitant


#' @export
blk_countPastEvents <- function(time, id=noID(time), event.t=time, id2=id) {
    .checkID(id)
    .checkID(id2)

    temp <- rbind(
        data.frame(id=id,  time=time,    code=1),
        data.frame(id=id2, time=event.t, code=0))

    temp <- temp[order(temp$id, temp$time, temp$code),]
    unlist(tapply(temp$code==0, temp$id, cumsum, simplify=FALSE))[temp$code==1]
}
#' @export
countPastEvents <- blk_countPastEvents


#' @export
blk_locf2 <- function(time, id=noID(time), observ.x, observ.t, id2=id, na.action=c("fill", "carry.back"), fill=NA) {
    .checkID(id)
    .checkID(id2)

    if (any(is.na(time))) {
        stop("time contains missing values")
    }
    if (any(is.na(observ.t))) {
        stop("observ.t contains missing values")
    }
    temp <- rbind(
        data.frame(id=id,  time=time,     code=1, x=NA),
        data.frame(id=id2, time=observ.t, code=0, x=observ.x))

    temp <- temp[order(temp$id, temp$time, temp$code),]
    blk_locf(temp$x, temp$id, na.action, fill)[temp$code==1]
}
#' @export
locf2 <- blk_locf2


#' @export
blk_nearestMatch <- function(time, id=noID(time), observ.x, observ.t, id2=id, direction=c("forward", "backward", "both"), tol=Inf, fill=NA, na.rm=FALSE, diff.op=difftime_default) {
    .checkID(id)
    .checkID(id2)

    if (any(is.na(time))) {
        stop("time contains missing values")
    }
    if (any(is.na(observ.t))) {
        stop("observ.t contains missing values")
    }
    if (length(time) != length(id)) {
        stop("Lengths of time and id must match")
    }
    if (length(observ.t) != length(id2)) {
        stop("Lengths of observ.t and id2 must match")
    }
    if (length(observ.x) != length(id2)) {
        stop("Lengths of observ.x and id2 must match")
    }
    if (!is_sorted(id, time)) {
        warning("time is not sorted within blocks of id.")
    }
    if (!is_sorted(id2, observ.t)) {
        warning("observ.t is not sorted within blocks of id2.")
    }

    diff.op <- .mymatchfun(diff.op)

    direction <- match.arg(direction)

    if (!is.numeric(tol)) {
        stop("tol must be numeric.")
    }
    if (length(tol) > 2 || length(tol) < 1) {
        stop("tol should have length 1 or 2.")
    }
    if (length(tol) == 1) {
        tol <- rep(tol, 2)
    }
    tol.f <- tol[1]
    tol.b <- tol[2]

    if (na.rm) {
        keep     <- !is.na(observ.x)
        observ.x <- observ.x[keep]
        observ.t <- observ.t[keep]
        id2      <- id2[keep]
    }

    temp <- rbind(
        data.frame(id=id,  time=time,     code=1, ix.f=NA,               time.f=time[NA], ix.b=NA,               time.b=time[NA]),
        data.frame(id=id2, time=observ.t, code=0, ix.f=.myseqalong(id2), time.f=observ.t, ix.b=.myseqalong(id2), time.b=observ.t))

    temp <- temp[order(temp$id, temp$time, temp$code),]

    if (direction=="forward" || direction=="both") {
        temp$ix.f   <- blk_locf(temp$ix.f,   temp$id, fill=NA)
        temp$time.f <- blk_locf(temp$time.f, temp$id, fill=NA)
    }
    if (direction=="backward" || direction=="both") {
        temp$ix.b   <- blk_nocb(temp$ix.b,   temp$id, fill=NA)
        temp$time.b <- blk_nocb(temp$time.b, temp$id, fill=NA)
    }

    if (direction=="forward") {
        dtfor <- ifelse(is.na(temp$time.f), NA, diff.op(temp$time, temp$time.f))
        ix <- ifelse(!is.na(dtfor) & dtfor <= tol.f, temp$ix.f, NA)
    } else if (direction=="backward") {
        dtbac <- ifelse(is.na(temp$time.b), NA, diff.op(temp$time.b, temp$time))
        ix <- ifelse(!is.na(dtbac) & dtbac <= tol.b, temp$ix.b, NA)
    } else {
        dtfor <- ifelse(is.na(temp$time.f), NA, diff.op(temp$time  , temp$time.f))
        dtbac <- ifelse(is.na(temp$time.b), NA, diff.op(temp$time.b, temp$time  ))
        ix <- ifelse(!is.na(dtfor) & (is.na(dtbac) | dtfor <= dtbac),
              ifelse(dtfor <= tol.f, temp$ix.f, NA),
              ifelse(!is.na(dtbac) & dtbac <= tol.b, temp$ix.b, NA))
    }
    ix1 <- ix[temp$code==1]
    res <- observ.x[ix1]
    res[is.na(ix1)] <- fill

    if (!is_sorted(id, time)) {
        res[order(order(id, time))]   # This will restore the order
    } else {
        res
    }
}
#' @export
nearestMatch <- blk_nearestMatch

#blk_totalDuration <- function(start.time, end.time, id=noID(start.time), id2=id, diff.op="-", fill=0) {
#    .checkID(id)
#
#    if (length(start.time) != length(end.time)) {
#        stop("start.time and end.time must have the same length")
#    }
#
#    if (!all(start.time <= end.time)) {
#        stop("All start.time must precede the corresponding end.time")
#    }
#
#    temp <- rbind(
#        data.frame(id=id, time=start.time, code=0, x=1),
#        data.frame(id=id, time=end.time,   code=1, x=-1))
#
#    temp   <- temp[order(temp$id, temp$time, temp$code),]
#
#    temp$y <- blk_applyNtoN(temp$x, temp$id, cumsum)
#    temp$y <- blk_shift(temp$y, id=temp$id, fill=0)
#    temp$z <- blk_diff(temp$time, id=temp$id, diff.op=diff.op, fill=0)
#
#    res <- blk_applyNto1((temp$y > 0) * temp$z, temp$id, sum, fill=fill)
#
#    if (!is.null(id2)) {
#        res <- blk_repeatValue(res, id=asID(levels(id)), id2=id2, fill=fill)
#    }
#    res
#}

#' @export
blk_timeAfter <- function(time, id, ind, na.action=c("neg.time", "fill"), fill=NA, diff.op=difftime_default) {
    .checkID(id)
    if (!is.logical(ind)) {
        warning("Expecting a logical ind")
        ind <- as.logical(ind)
    }
    if (length(ind) != length(id)) {
        stop("Lengths of ind and id must match")
    }
    if (length(time) != length(id)) {
        stop("Lengths of time and id must match")
    }
    na.action <- match.arg(na.action)
    diff.op <- .mymatchfun(diff.op)

    temp <- time
    temp[!ind] <- NA
    temp <- blk_locf(temp, id, "carry.back")
    tdiff <- diff.op(time, temp)

    if (na.action == "fill") {
        tdiff[tdiff < 0] <- fill
    }
    return(tdiff)
}

#' @export
timeAfter <- blk_timeAfter

#' @export
blk_timeAfter2 <- function(time, id, time2, id2, na.action=c("neg.time", "fill"), fill=NA, diff.op=difftime_default) {
    .checkID(id)
    .checkID(id2)

    temp <- rbind(
        data.frame(id=id,  time=time,   code=0),
        data.frame(id=id2, time=time2, code=1))

    temp <- temp[order(temp$id, temp$time, temp$code),]
    blk_timeAfter(temp$time, id=temp$id, ind=temp$code==1, na.action=na.action, fill=fill, diff.op=diff.op)[temp$code==0]
}

#' @export
timeAfter2 <- blk_timeAfter2
