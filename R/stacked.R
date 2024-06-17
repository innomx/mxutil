.myseqalong <- match.fun("seq_along")
.mymatchfun <- match.fun


# Subtracts time2 from time1.
# If both are of class POSIXt, gives the result in hours.
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


difftime_hours <- function(time1, time2) {
    if (length(time1) != length(time2)) {
        stop("Length of time1 and time2 should match.")
    }
    if (!(inherits(time1, "POSIXt") && inherits(time2, "POSIXt"))) {
        stop("Both time1 and time2 should be POSIXt objects")
    }
    as.double(difftime(time1, time2), units="hours")
}


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

restoreOrder <- function(o) {
    invord <- numeric(length(o))
    invord[o] <- 1:length(o)
    invord
}

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
singleValue <- blk_singleValue


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
repeatValue <- blk_repeatValue


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
applyNto1 <- blk_applyNto1


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
applyNtoN <- blk_applyNtoN


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
firstOnly <-  blk_firstOnly
blk_onlyFirst <- blk_firstOnly
onlyFirst <- blk_onlyFirst


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
lastOnly <-  blk_lastOnly
blk_onlyLast <- blk_lastOnly
onlyLast <- blk_onlyLast


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
firstOnwards <- blk_firstOnwards


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
lastOnwards <- blk_lastOnwards


blk_untilFirst <- function(id=noID(ind), ind, include.first=TRUE, fill=FALSE) {
    !blk_firstOnwards(id, ind, !include.first, !fill)
}
untilFirst <- blk_untilFirst


blk_untilLast <- function(id=noID(ind), ind, include.last=TRUE, fill=FALSE) {
    !blk_lastOnwards(id, ind, !include.last, !fill)
}
untilLast <- blk_untilLast


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
locf <- blk_locf


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
nocb <- blk_nocb


blk_sequentialID <- function(id) {
    cumsum(blk_firstOnly(id=id))
}
sequentialID <- blk_sequentialID


blk_isConstant <- function(x, id=noID(x), ind=NULL) {
    all(blk_applyNto1(x, id, nunique, ind=ind, fill=0) <= 1)
}
isConstant <- blk_isConstant

blk_isVarying <- function(x, id=noID(x), ind=NULL) {
    any(blk_applyNto1(x, id, nunique, ind=ind, fill=0) > 1)
}
isVarying <- blk_isVarying


blk_count <- function(id=noID(ind), id2=id, ind=NULL) {
    res <- blk_applyNto1(x=rep(1, length(id)), id=id, apply.fun=sum, ind=ind, fill=0)
    if (!is.null(id2)) {
        res <- blk_repeatValue(res, id=asID(levels(id)), id2=id2, fill=0)
    }
    res
}
count <- blk_count


blk_changeFromBaseline <- function(x, id=noID(x), baseline=blk_singleValue(x, id, ind=bl.ind, select=bl.select), bl.ind=NULL, bl.select="first", diff.op="-", ...) {
    baseline <- blk_repeatValue(baseline, id2=id)
    .mymatchfun(diff.op)(x, baseline, ...)
}
changeFromBaseline <- blk_changeFromBaseline


blk_relativeTime <- function(time, id=noID(time), baseline=blk_singleValue(time, id, ind=bl.ind, select=bl.select), bl.ind=NULL, bl.select="first", diff.op=difftime_default, ...) {
    blk_changeFromBaseline(time, id, baseline, diff.op=diff.op, ...)
}
relativeTime <- blk_relativeTime


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
shift <- blk_shift


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
diff2 <- blk_diff


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


blk_findConsecutive <- function(id=noID(ind), ind, min.consec=2, fill=NA) {
    .checkID(id)

    y <- unlist(tapply(ind, id, .consecLengths, simplify=FALSE))
    z <- rep(0, length(ind))
    z[ind] <- rep(ifelse(y >= min.consec, 1, 0), times=y)
    y <- y[y >= min.consec]
    z[z==1] <- rep(.myseqalong(y), times=y)
    z <- z - blk_applyNtoN(z-1, id, min, ind=(z > 0), fill=0)
    z[z == 0] <- fill
    z
}
findConsecutive <- blk_findConsecutive


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
concomitant <- blk_concomitant


blk_countPastEvents <- function(time, id=noID(time), event.t=time, id2=id) {
    .checkID(id)
    .checkID(id2)

    temp <- rbind(
        data.frame(id=id,  time=time,    code=1),
        data.frame(id=id2, time=event.t, code=0))

    temp <- temp[order(temp$id, temp$time, temp$code),]
    unlist(tapply(temp$code==0, temp$id, cumsum, simplify=FALSE))[temp$code==1]
}
blk_countPastEvents <- blk_countPastEvents


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
locf2 <- blk_locf2


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
    if (!is.sorted(id, time)) {
        warning("time is not sorted within blocks of id.")
    }
    if (!is.sorted(id2, observ.t)) {
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

    if (!is.sorted(id, time)) {
        res[order(order(id, time))]   # This will restore the order
    } else {
        res
    }
}
nearestMatch <- blk_nearestMatch

#blk_totalDuration <- function(start.time, end.time, id=noID(start.time), id2=id, diff.op="-", fill=0) {
#    PCSmisc:::.checkID(id)
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
