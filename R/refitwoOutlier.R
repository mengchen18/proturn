#' @title Using subset of rows to refit NLS model
#' @description Used by shinyapp.
#' @param x an element in "list" returned by \code{\link{fitNLSModels}}
#' @param include which rows should be included
#' @param k the synthesis/degradation constant
#' @param t The time point (hours)
#' @param tcc The doubling time of cells. By default this value is Inf, which means the cells are in steady state
#' @param A optinal argument for fixed A, if this argument is given, "A" won't be optimized
#' @param B optional argument for fixed B, if this argument is given, "B" won't be optimized
#' @param par.init The initial values of parameters to be optimized, it should be list of three elements
#'   names as "A", "B" and "k".
#' @param par.lower The lower boundary of parameters to be optimized, it should be a numeric values with
#'   length 3 and named as "A", "B" and "k".
#' @param par.upper The upper boundary of parameters to be optimized, it should be a numeric values with
#'   length 3 and named as "A", "B" and "k".
#' @param message A logical value to indicated if any message should be printed
#' @param fitIndividual A logical value. When multiple lines are combined to fit a single model (see f),
#'   whether each individual lines should also be fitted. Only used when x is an object of class \code{matrix}
#' @return
#'  a vector of optimized parameters, including A, B, kd, confidence intervals (2.5% and 97.5%),
#'  mean square error and r-square values.
#'  In addition, if individual rows are fitted, the object also contains an attribute stores
#'  parameters fitted on each individual row.
#' @export
#' @keywords internal
#'
#' @examples
#' tp <- c(0, 1, 2, 4, 8, 16, 32, 64)
#' A <- runif(5, min = 0.75, max = 0.95)
#' B <- runif(5, min = 0.1, max = 0.15)
#' kd <- runif(5, min = 0.05, max = 0.1)
#' ds <- mapply(function(A, B, kd) {
#'   degCurve(A=A, B=B, kd=kd, tcc=Inf, t = tp) + rnorm(length(tp), sd = 0.03)
#' }, A, B, kd)
#' ds <- t(ds)
#' ds <- rbind(ds, runif(8, 0, 1))
#' rf <- c("p1", "p2", "p2", "p3", "p3", "p3")
#' ss <- fitNLSModels(x = ds, f = rf, t = tp, tcc = Inf, type = "deg",
#'                    par.init = list(A=0.8, B=0.2, kd=0.04),
#'                    par.lower=c(A=0, B=0, kd=0),
#'                    par.upper=c(A=1, B=1, kd=10))
#'
#' plotCurve.comb(x = ss$list$p3, t = tp, tcc = Inf, curve = "deg")
#' #' only use the first two rows
#' res <- refitwoOutlier(x = ss$list$p3, include = 1:2, t = tp, tcc = Inf)
#' plotCurve.comb(x = res, t = tp, tcc = Inf, curve = "deg")
#'
#'
#' # synthesis curve
#' ds2 <- mapply(function(A, B, kd) {
#'   synCurve(A=A, B=B, kd=kd, tcc=Inf, t = tp) + rnorm(length(tp), sd = 0.03)
#' }, A, B, kd)
#' ds2 <- t(ds2)
#' ds2 <- rbind(ds2, runif(8, 0, 1))
#' ss2 <- fitNLSModels(x = ds2, f = rf, t = tp, tcc = Inf, type = "syn",
#'                     par.init = list(A=0.8, B=0.2, ks=0.04),
#'                     par.lower=c(A=0, B=0, ks=0),
#'                     par.upper=c(A=1, B=1, ks=10))
#'
#' plotCurve.comb(x = ss2$list$p3, t = tp, tcc = Inf, curve = "syn")
#' #'#'
#' res <- refitwoOutlier(x = ss2$list$p3, include = 1:2, t = tp, tcc = Inf)
#' plotCurve.comb(x = res, t = tp, tcc = Inf, curve = "syn")
refitwoOutlier <- function(x, include, t, tcc,
                           A = NULL, B = NULL, k = NULL,
                           par.init = list(A=1, B=0.2, k=0.04),
                           par.lower=c(A=0, B=0, k=0),
                           par.upper=c(A=1, B=1, k=10),
                           fitIndividual = TRUE,
                           message = TRUE) {

  ind <- attr(x, "individual")
  mat <- attr(x, "inputmatrix")

  if (attr(x, "type") == "Degradation") {
    func <- fitDegNLS
    names(par.init)[names(par.init) == "k"] <- "kd"
    names(par.lower)[names(par.lower) == "k"] <- "kd"
    names(par.upper)[names(par.upper) == "k"] <- "kd"
  } else if (attr(x, "type") == "Synthesis") {
    func <- fitSynNLS
    names(par.init)[names(par.init) == "k"] <- "ks"
    names(par.lower)[names(par.lower) == "k"] <- "ks"
    names(par.upper)[names(par.upper) == "k"] <- "ks"
  } else {
    stop("unknown type of fit.")

  }

  func(x = mat[include, , drop = FALSE], t = t, tcc = tcc,
       par.init = par.init, par.lower=par.lower,
       par.upper = par.upper,
       fitIndividual = fitIndividual)
}


#' @title rounds the values in a data.frame
#' @description for a data.frame or list, round all numerical variables to specified number of significant digits
#' @param x a \code{data.frame} or \code{list}
#' @param digits teh number of significant digits
#' @return the same object and size as x, but all numerical values are rounded
#' @export
#' @keywords internal
#' @examples 
#' a <- data.frame(num = rnorm(4),
#'                 char = letters[1:4])
#' sigDF(a, 3)

sigDF <- function(x, digits = 3) {
  if (is.matrix(x) || is.numeric(x))
    return(signif(x, digits = digits))
  if (inherits(x, c("data.frame", "list"))) {
    i <- sapply(x, is.numeric)
    i <- i & !sapply(x, is.integer)
    x[i] <- lapply(x[i], signif, digits = digits)
  } else
    stop("unknown data type")
  x
}
