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
    if (!file.exists(file)) {
        file <- paste0(file, ".csv")
    }
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



#' Creatinine Clearance Using the Cockroft-Gault Formula
#'
#' Computes estimated creatinine clearance using the Cockroft-Gault
#' formula.
#'
#' @param ... Ignored. This forces all parameters to be named, which is safer.
#' @param scr.mg.dL Numeric. Serum creatinine in mg/dL.
#' @param weight.kg Numeric. Body weight in kilograms.
#' @param age.yr Numeric. Age in years.
#' @param is.female Logical. TRUE for females, FALSE for males.
#' @param scr.umol.L Numeric. Alternative specification of serum creatinine in umol/L.
#' @param scr.conversion.factor Numeric. A factor for converting serum creatinine from umol/L to mg/dL.
#'
#' @note
#'   All arguments can be vectors.
#'
#' @references
#'   - https://en.wikipedia.org/wiki/Glomerular_filtration_rate
#'   - Cockcroft DW, Gault MH (1976). "Prediction of creatinine clearance from
#'   serum creatinine". Nephron. 16 (1): 31–41. doi:10.1159/000180580. PMID
#'   1244564.
#'
#' @return
#'   A numeric vector containing the estimated creatinine clearance in mL/min.
#'
#' @author Benjamin Rich <mail@benjaminrich.net>
#' 
#' @seealso
#'   [egfr.mdrd()]
#'
#' @examples
#' 
#' # With serum creatinine in mg/dL
#' crcl.cg(scr.mg.dL=0.9, weight.kg=112.9, age.yr=53, is.female=FALSE)
#' 
#' # With serum creatinine in umol/L
#' crcl.cg(scr.umol.L=0.9*88.4, weight.kg=112.9, age.yr=53, is.female=FALSE)
#'
#' @export
crcl.cg <- function(..., scr.mg.dL, weight.kg, age.yr, is.female, scr.umol.L, scr.conversion.factor=88.4017) {
    if (missing(scr.mg.dL) || is.null(scr.mg.dL)) {
        if (missing(scr.umol.L) || is.null(scr.umol.L)) {
            stop("One of scr.mg.dL and scr.umol.L must be present")
        } else {
            # Either way, this gives similar (but not identical) results.
            scr.mg.dL <- (scr.umol.L / scr.conversion.factor)
            #(1.23 * (140 - age.yr) * weight.kg * (1.04/1.23)^is.female / (scr.umol.L))
        }
    } else {
        if (!(missing(scr.umol.L) || is.null(scr.umol.L))) {
            stop("Only one of scr.mg.dL and scr.umol.L must be present")
            #scr.mg.dL[is.na(scr.mg.dL)] <- (scr.umol.L[is.na(scr.mg.dL)] / scr.conversion.factor)
        }
    }
    (140 - age.yr) * weight.kg * (0.85)^is.female / (72 * scr.mg.dL)
}


#' Estimated GFR Using the MDRD Formula
#'
#' Computes estimated glomerular filtration rate (GFR) using the
#' Modification of Diet in Renal Disease (MDRD) formula.
#'
#' @param ... Ignored. This forces all parameters to be named, which is safer.
#' @param scr.mg.dL Numeric. Serum creatinine in mg/dL.
#' @param age.yr Numeric. Age in years.
#' @param is.female Logical. TRUE for females, FALSE for males.
#' @param is.black Logical. TRUE if race is 'Black', FALSE otherwise.
#' @param bun.mg.dL Numeric. Blood urea nitrogen (BUN) concentrations in mg/dL (optional).
#' @param albumin.g.dL Numeric. Albumin concentration in g/dL (optional).
#' @param scr.umol.L Numeric. Alternative specification of serum creatinine in umol/L.
#' @param scr.conversion.factor Numeric. A factor for converting serum creatinine from umol/L to mg/dL.
#' @param idms.calibrated Logical. Indicates whether serum creatinine measurements are IDMS calibrated
#'        (see https://en.wikipedia.org/wiki/Glomerular_filtration_rate)
#'
#' @note
#'   All arguments can be vectors.
#'   Albumin and BUN are optional.  If they are specified, a more elaborate version
#'   of the MDRD equation is used.
#'
#' @return
#'   A numeric vector containing the estimated GFR mL/min/1.73m².
#'
#' @references
#'   https://en.wikipedia.org/wiki/Glomerular_filtration_rate
#'
#' @author Benjamin Rich <mail@benjaminrich.net>
#' 
#' @seealso
#'   [crcl.cg()]
#'
#' @examples
#' 
#' # With serum creatinine in mg/dL
#' egfr.mdrd(scr.mg.dL=0.9, age.yr=53, is.female=FALSE, is.black=FALSE)
#' 
#' # With serum creatinine in umol/L
#' egfr.mdrd(scr.umol.L=0.9*88.4, age.yr=53, is.female=FALSE, is.black=FALSE)
#' 
#' @export
egfr.mdrd <- function(..., scr.mg.dL, age.yr, is.female, is.black, bun.mg.dL=NULL, albumin.g.dL=NULL, scr.umol.L, scr.conversion.factor=88.4017, idms.calibrated=TRUE) {
    if (missing(scr.mg.dL) || is.null(scr.mg.dL)) {
        if (missing(scr.umol.L) || is.null(scr.umol.L)) {
            stop("One of scr.mg.dL and scr.umol.L must be present")
        } else {
            scr.mg.dL <- (scr.umol.L / scr.conversion.factor)
        }
    } else {
        if (!(missing(scr.umol.L) || is.null(scr.umol.L))) {
            stop("Only one of scr.mg.dL and scr.umol.L must be present")
            #scr.mg.dL[is.na(scr.mg.dL)] <- (scr.umol.L[is.na(scr.mg.dL)] / scr.conversion.factor)
        }
    }
    if (!is.null(bun.mg.dL) && !is.null(albumin.g.dL)) {
        egfr <- 170 * scr.mg.dL^(-0.999) * age.yr^(-0.176) * (1.18)^is.black * (0.762)^is.female * bun.mg.dL^(-0.17) * albumin.g.dL^(0.318)
    } else {
        egfr <- 186 * scr.mg.dL^(-1.154) * age.yr^(-0.203) * (1.212)^is.black * (0.742)^is.female
    }
    if (idms.calibrated) {
        egfr <- egfr * (175/186)
    }
    egfr
}

#' Ideal Body Weight Using the Devine Formula
#'
#' Computes ideal body weight (IBW) using the Devine (1974) formula.
#'
#' @param ... Ignored. This forces all parameters to be named, which is safer.
#' @param ht.cm Numeric. Height in centimeters.
#' @param is.male Logical. TRUE for males, FALSE for females.
#'
#' @note
#'   All arguments can be vectors.
#'
#' @references
#'   Devine D. Case study number 25 gentamicin therapy. Drug Intell Clin Pharm 1974; 8: 650–5.
#'
#' @return
#'   A numeric vector containing the IBW kg.
#'
#' @author Benjamin Rich <mail@benjaminrich.net>
#' 
#' @examples
#' 
#' ibw.devine(ht.cm=c(183, 165), is.male=c(TRUE, FALSE))
#' 
#' @export
ibw.devine <- function(..., ht.cm, is.male) {
    45.4 + 0.89*(ht.cm - 152.4) + 4.5*is.male
}

#' Body Surface Area Using the DuBois and DuBois Formula
#'
#' Computes estimated body surface area (BSA) using the DuBois and DuBois
#' (1916) formula.
#'
#' @param ... Ignored. This forces all parameters to be named, which is safer.
#' @param wt.kg Numeric. Body weight in kilograms.
#' @param ht.cm Numeric. Height in centimeters.
#'
#' @note
#'   All arguments can be vectors.
#'
#' @references
#'   DuBois D. and DuBois E.F. A formula to estimate the approximate surface
#'   area if height and weight be known. Arch Intern Med, 17:863–871, 1916.
#'
#' @return
#'   A numeric vector containing the estimated BSA in m².
#'
#' @author Benjamin Rich <mail@benjaminrich.net>
#' 
#' @examples
#' 
#' bsa.dubois(wt.kg=c(79, 82), ht.cm=c(183, 165))
#' 
#' @export
bsa.dubois <- function(..., wt.kg, ht.cm) {
    0.007184*(wt.kg^0.425)*(ht.cm^0.725)
}

#' Lean Body Mass Using the James Equation
#'
#' Computes Lean body mass (LBM) Using the James (1977) equation.
#'
#' @param ... Ignored. This forces all parameters to be named, which is safer.
#' @param wt.kg Numeric. Body weight in kilograms.
#' @param ht.cm Numeric. Height in centimeters.
#' @param is.male Logical. TRUE for males, FALSE for females.
#'
#' @note
#'   All arguments can be vectors.
#'
#' @references
#'   - James WP. Research on obesity. Nutrition Bulletin. 1977 Sep;4(3):187-90.
#'   - <https://www.calculator.net/lean-body-mass-calculator.html>
#'
#' @return
#'   A numeric vector containing the estimated LBM kg.
#'
#' @author Benjamin Rich <mail@benjaminrich.net>
#' 
#' @examples
#' 
#' lbm.james(wt.kg=c(79, 82), ht.cm=c(183, 165), is.male=c(TRUE, FALSE))
#' 
#' @export
lbm.james <- function(..., wt.kg, ht.cm, is.male) {
    ifelse(is.male,
        1.10*wt.kg - 128*(wt.kg/ht.cm)^2,
        1.07*wt.kg - 148*(wt.kg/ht.cm)^2
    )
}

#' Skeletal Muscle Mass Using the Lee Formula
#'
#' Computes skeletal muscle mass (SM) using the Lee et al. (2000) formula.
#'
#' @param ... Ignored. This forces all parameters to be named, which is safer.
#' @param wt.kg Numeric. Body weight in kilograms.
#' @param ht.cm Numeric. Height in centimeters.
#' @param age.yr Numeric. Age in years.
#' @param is.male Logical. TRUE for males, FALSE for females.
#' @param is.asian Logical. TRUE if race is 'Asian', FALSE otherwise.
#' @param is.black Logical. TRUE if race is 'Black', FALSE otherwise.
#'
#' @note
#'   All arguments can be vectors.
#'
#' @references
#'   Lee RC, Wang Z, Heo M, Ross R, Janssen I, Heymsfield SB. Total-body
#'   skeletal muscle mass: development and cross-validation of anthropometric
#'   prediction models. The American journal of clinical nutrition. 2000 Sep
#'   1;72(3):796-803.
#'
#' @return
#'   A numeric vector containing the estimated SM in.
#'
#' @author Benjamin Rich <mail@benjaminrich.net>
#' 
#' @examples
#' 
#' sm.lee(wt.kg=79, ht.cm=183, age.yr=28, is.male=TRUE, is.asian=FALSE, is.black=TRUE)
#' 
#' @export
sm.lee <- function(..., wt.kg, ht.cm, age.yr, is.male, is.asian, is.black) {

    if (any(is.asian & is.black)) {
        stop("An individual cannot be both Asian and Black")
    }

    0.226*wt.kg + 13.0*(ht.cm/100) - 0.089*age.yr + 6.3*is.male - 1.6*is.asian + 1.9*is.black - 11.0
}


#' Compute Age From Date-of-Birth
#'
#' Computes the age of a person (in years) on a specified date from the
#' person's date of birth.
#'
#' @param on.this.date A vector of [`base::Date`] objects.  The
#'   date on which the age is to be computed.
#' @param birth.date A vector of [`base::Date`] objects.  Th
#'   date-of-birth.
#' @return
#'   An integer vector containing the age in years.
#'
#' @author Benjamin Rich <mail@benjaminrich.net>
#' 
#' @seealso
#'   [base::difftime()]
#'
#' @examples
#' age.in.years(
#'   as.Date(c("2012-03-01", "2012-02-29", "2012-02-28")),
#'   as.Date(c("1980-02-03", "1980-02-29", "1980-05-28")))
#'
#' @export
age.in.years <- function(on.this.date, birth.date) {
    this.yr <- as.numeric(format(on.this.date, "%Y"))
    birth.yr <- as.numeric(format(birth.date, "%Y"))
    elap.yr <- this.yr - birth.yr

    bday.this.year <- as.Date(paste(this.yr, format(birth.date, "%m-%d"), sep="-"))

    as.integer(ifelse(on.this.date < bday.this.year, elap.yr - 1, elap.yr))
}



#' @export
blk_isSteadyState <- function(time, id=noID(time), dose.ind, dose.interval, tol, min.time.ss, dose=NULL, diff.op=difftime_default) {
    .checkID(id)

    if (is.null(dose)) {
        delta.dose <- 0
    } else {
        if (!is.numeric(dose)) {
            stop("'dose' must be numeric times.")
        }
        delta.dose <- blk_diff(dose, id=id, ind=dose.ind, fill=0)
    }
    itime <- blk_intereventTime(time, id=id, ind=dose.ind, fill=Inf, diff.op=diff.op)

    cumtime <- cumsum(ifelse(dose.ind, dose.interval, 0))
    flag <- !dose.ind | ((delta.dose == 0) & (abs(itime - dose.interval) <= tol))
    cumtime2 <- blk_timeAfter(cumtime, id, !flag, fill=0)
    ss <- dose.ind & (cumtime2 >= min.time.ss)
    ss

    #min.consec <- ceiling(min.time.ss/dose.interval)
    #cons <- blk_findConsecutive(id, ind=((delta.dose == 0) & (abs(itime - dose.interval) <= tol)), min.consec=min.consec, fill=0)
    #ss <- blk_shift(cons > 0, id=id, shift.by=min.consec-1, ind=cons > 0, fill=FALSE)
    #ss
}

#' @export
isSteadyState <- blk_isSteadyState

#' @export
blk_ii <- function(time, id=noID(time), dose.ind, diff.op=difftime_default) {

    ii <- blk_intereventTime(time, id=id, ind=dose.ind, diff.op=diff.op)
    ii <- blk_shift(ii, id=id, shift.by= -1, ind=dose.ind, fill=0)
    ii[!dose.ind] <- NA
    ii
}

#' @export
interdoseInterval <- blk_ii

#' @export
blk_addl <- function(ii, id=noID(dose.ind), dose.ind, dose, tol.ii=1e-5, tol.dose=1e-5, min.consec=2, include.last=FALSE, diff.op=difftime_default) {

    .checkID(id)

    if (!is.numeric(dose)) {
        stop("'dose' must be numeric times.")
    }
    if (!(length(ii)==length(id) && length(dose.ind)==length(id) && length(dose)==length(id))) {
        stop("'ii', 'id', 'dose.ind' and 'dose' must all have the same length.")
    }

    y <- blk_diff(ii, id=id, ind=dose.ind)
    y <- blk_shift(y, id=id, shift.by= -1, ind=dose.ind, fill=0)

    z <- blk_diff(dose, id=id, ind=dose.ind)
    z <- blk_shift(z, id=id, shift.by= -1, ind=dose.ind, fill=0)

    flag <- dose.ind & (abs(y) < tol.ii) & (abs(z) < tol.dose)

    # Do this twice, because the last ii can be anything
    for (i in 1:2) {
        flag2 <- blk_shift(flag, id=id, fill=FALSE)
        flag[dose.ind & !flag & flag2] <- flag2[dose.ind & !flag & flag2]
    }

    if (!include.last) {
        flag <- flag & blk_shift(dose.ind, id=id, shift.by= -1, fill=FALSE)
    }

    cons <- blk_findConsecutive(id, ind=(dose.ind & flag), min.consec=min.consec, fill=0)
    i <- (cons > 0)
    id2 <- asID(cons[i])
    addl <- ifelse(dose.ind, 0, NA)
    addl[i] <- ifelse(blk_firstOnly(id2), blk_count(id=id2) - 1, -1)
    addl
}

#' @export
addlDoses <- blk_addl

#' @export
blk_noninformativeDose <- function(id=noID(dose.ind), dose.ind, ss.ind=NULL, obs.ind=!dose.ind) {

    # Description: Flags doses preceding a steady state event when no
    # observation event has occurred between the two.

    .checkID(id)
    x <- rep.int(FALSE, length(id))
    if (!is.null(ss.ind)) {
        cons <- blk_findConsecutive(id, ind=!obs.ind, fill=0)
        i <- (cons > 0)
        id2 <- asID(cons[i])
        x[i] <- blk_untilLast(id=id2, ind=ss.ind[i], include.last=FALSE)
    }
    dose.ind & (x | blk_lastOnwards(id=id, ind=obs.ind, include.last=FALSE))
}

#' @export
noninformativeDose <- blk_noninformativeDose
