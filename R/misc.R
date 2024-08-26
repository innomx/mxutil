#' @export
names_tolower <- function(x) setNames(x, tolower(names(x)))

#' @export
names_toupper <- function(x) setNames(x, toupper(names(x)))

#' @export
read_analysis_dataset <- function(x, path=getOption("path.analysis.dataset", getwd()), ...) {
    if (is.null(path)) {
        path <- getwd()
    }
    file <- file.path(path, x)
    names_tolower(data.table::fread(file, ..., na.strings=c(".", "", "NA")))
}

#' @export
write_analysis_dataset <- function(x, filebase, path=getOption("path.analysis.dataset", getwd()), datestamp=FALSE, ...) {
    quote <- FALSE
    for (i in seq_len(ncol(x))) {
        if (is.character(x[[i]]) || is.factor(x[[i]])) {
            if (any(grepl("\"", as.character(x[[i]])))) {
                warning(paste0("Column ", i, " contains quotes."))
                quote <- TRUE
            }
            if (any(grepl(",", as.character(x[[i]])))) {
                warning(paste0("Column ", i, " contains commas"))
                quote <- TRUE
            }
        }
    }
    if (datestamp) {
        file <- sprintf("%s_%s.csv", filebase, format(Sys.time(), "%Y%m%d"))
    } else {
        file <- sprintf("%s.csv", filebase)
    }
    if (is.null(path)) {
        path <- getwd()
    }
    file <- file.path(path, file)
    if (!dir.exists(dirname(file))) {
        message(paste0("Creating directory ", dirname(file)))
        dir.create(dirname(file), showWarnings=TRUE, recursive=TRUE)
    }
    message(paste0("Writing file ", file))
    data.table::fwrite(names_toupper(x), file, quote=quote, na=".")
}

#' Replace Values
#'
#' @examples
#' x <- c("a", "", NA, "b")
#' Replace(x, "Missing", z=NA)
#' Replace(x, NA, z="")
#' Replace(x, "Missing", z="")
#' Replace(x, "Missing", z=c("", NA))
#' 
#' x %|% "Missing"
#' 
#' @export
Replace <- function(x, y, z=NA) {
    y <- rep(y, len=length(x))
    z <- unique(z)
    i <- if (length(z) > 1) {
        ifelse(is.na(x), any(is.na(z)), x == z[!is.na(z)])
    } else if (is.na(z)) {
        is.na(x)
    } else {
        ifelse(is.na(x), FALSE, x == z)
    }
    base::replace(x, i, y[i])
}

#' @rdname Replace
`%|%` <- function(x, y) {
    Replace(x, y, z=NA)
}

#' Extract Match of Regular Exptression From String
#' @export
extract_match <- function(re, x, capture=T, simplify=T, perl=T, ...) {
    s <- regexpr(re, x, perl=perl, ...)
    myregmatches <- function(x, m) {
        so <- as.numeric(m)
        so[so < 0] <- NA
        eo <- so + attr(m, "match.length") - 1L
        substring(x, so, eo)
    }
    if (!is.null(attr(s, "capture.start")) && isTRUE(capture)) {
        m <- as.matrix(mapply(function(cs, cl) { myregmatches(x, structure(cs, match.length=cl)) },
            as.data.frame(attr(s, "capture.start")), as.data.frame(attr(s, "capture.length"))))
        if (ncol(m) == 1 && isTRUE(simplify)) {
            drop(m)
        } else {
            dimnames(m) <- dimnames(attr(s, "capture.start"))
            m
        }
    } else {
        myregmatches(x, s)
    }
}

#' Extract Number From String
#' @export
extract_number_from_string <- function(x, convert=T) {
    re <- "-?([[:digit:]]+(\\.[[:digit:]]*)?|\\.[[:digit:]]+)"
    re <- paste0(re, "([eE]-?[[:digit:]]+)?")
    re <- paste0("(^|[[:space:]])(", re, ")($|[[:space:]])")
    x <- extract_match(re, x, capture=FALSE)
    if (isTRUE(convert)) {
        as.numeric(x)
    } else {
        x
    }
}

#' Extract Number From String
#' @examples
#' extract_clocktime_from_dtc(c("2024-01-01T12:16", "2024-01-01T12:1x"))
#' @export
extract_clocktime_from_dtc <- function(x) {
    re <- "T([[:digit:]]{2}:[[:digit:]]{2})"
    extract_match(re, x)
}

#' How Many Missing Values?
#'
#' Returns the number of `NA`'s in a vector.
#'
#' @param x A vector.
#' @return An integer equal to `sum(is.na(x))`.
#' @export
nmissing <- function(x) {
    sum(is.na(x))
}

# How many unique values?
#' How Many Unique Values?
#'
#' Returns the number of unique values in a vector.
#'
#' @param x A vector.
#' @return An integer equal to `length(unique(x))`.
#' @export
nunique <- function(x, na.rm=FALSE) {
    if (na.rm) {
        length(unique(x[!is.na(x)]))
    } else {
        length(unique(x))
    }
}
