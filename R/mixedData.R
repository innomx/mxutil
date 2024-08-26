
#' Mixed Numeric/Non-Numeric Data
#'
#' A simple S3 class for a vector that can contain both numeric and non-numeric
#' values.
#' @export
mixedData <- function(x) {
    chr <- as.character(x)
    num <- suppressWarnings(as.numeric(chr))
    obj <- chr
    attr(obj, "numeric") <- num
    class(obj) <- "mixedData"
    obj
}

#' @export
as_character <- function(obs) UseMethod("as_character")

#' @export
as_numeric <- function(obs) UseMethod("as_numeric")

#' @export
as_double <- function(obs) UseMethod("as_double")

#' @export
as_character.default <- function(obj, ...) {
    as.character(obj, ...)
}

#' @export
as_numeric.default <- function(obj, ...) {
    as_double(obj, ...)
}

#' @export
as_double.default <- function(obj, ...) {
    as.double(obj, ...)
}

#' @export
as_character.mixedData <- function(obj, ...) {
    as.character(unclass(obj), ...)
}

#' @export
as_double.mixedData <- function(obj, ..., fill=NA) {
    num <- attr(obj, "numeric", exact=TRUE)
    num[is.na(num) & !is.na(obj)] <- fill
    num
}

#' @export
print.mixedData <- function(x, ...) {
    print(as_character(x), quote=FALSE)
}


#' @export
as.data.frame.mixedData <- function(x, ...) {
    nm <- deparse(substitute(x), width.cutoff = 500L)
    as.data.frame.vector(x, ..., nm = nm)
}

#' @export
unique.mixedData <- function(x, ...) {
    num <- attr(x, "numeric", exact=TRUE)
    chr <- as.character(unclass(x))
    unique(ifelse(is.na(num), chr, "<Number>"), ...)
}

#' @export
table.mixedData <- function(x, ...) {
    num <- attr(x, "numeric", exact=TRUE)
    chr <- as.character(unclass(x))
    table(ifelse(is.na(num), chr, "<Number>"), ...)
}
