#' @title Fitting synthesis curve using NLS algorithm
#' @description fit protein synthesis curve using NLS algorithm
#' @param x A numeric vector
#' @param t The time point (hours)
#' @param tcc The doubling time of cells. By default this value is Inf, which means the cells are in steady state
#' @param A optinal argument for fixed A, if this argument is given, "A" won't be optimized
#' @param B optional argument for fixed B, if this argument is given, "B" won't be optimized
#' @param ks optional argument for fixed ks, if this argument is given, "ks" won't be optimized
#' @param par.init The initial values of parameters to be optimized, it should be list of three elements
#'   names as "A", "B" and "ks".
#' @param par.lower The lower boundary of parameters to be optimized, it should be a numeric values with
#'   length 3 and named as "A", "B" and "ks".
#' @param par.upper The upper boundary of parameters to be optimized, it should be a numeric values with
#'   length 3 and named as "A", "B" and "ks".
#' @param message A logical value to indicated if any message should be printed
#' @param fitIndividual A logical value, whether each individual row should also be fitted.
#'   Only used when x is an object of class \code{matrix}
#' @details
#'   More information about the fitted model could be find in \code{\link{fitNLSModels}}.
#' @export
#' @return
#'  a vector of optimized parameters, including A, B, ks, confidence intervals (2.5% and 97.5%),
#'  mean square error and r-square values.
#'  In addition, if individual rows are fitted, the object also contains an attribute stores
#'  parameters fitted on each individual row.
#' @examples
#' # synthesi curve
#' tp <- c(0, 1, 2, 4, 8, 16, 32, 64)
#' ratios <- synCurve(A=0.85, B = 0.1, kd=0.5, tcc=Inf, t = tp) + rnorm(length(tp), sd = 0.05)
#'
#' r <- fitSynNLS(ratios, t = tp, tcc = Inf)
#' plotCurve(ratios, tp, tcc = Inf, A = r[["A"]], B = r[["B"]], k = r[["ks"]],
#'           add = FALSE, curve = "syn", err.x = log(2)/r[c("ci025", "ci975")],
#'           err.y = synCurve(A = r[["A"]], B = r[["B"]], kd = r[["ks"]], t = log(2)/ r[["ks"]], tcc = Inf))
#'
#'
#' ratio2 <- rbind(p1 = ratios + rnorm(length(ratios), sd = 0.2),
#'                 p2 = ratios + rnorm(length(ratios), sd = 0.2))
#' r.mat <- fitSynNLS(ratio2, t = tp, tcc = Inf)
#' plotCurve(ratio2, rep(tp, 2), tcc = Inf, A = r.mat[["A"]], B = r.mat[["B"]], k = r.mat[["ks"]],
#'           add = FALSE, curve = "syn", err.x = log(2)/r.mat[c("ci025", "ci975")],
#'           err.y = degCurve(A = r.mat[["A"]], B = r.mat[["B"]], kd = r.mat[["ks"]],
#'                            t = log(2)/ r.mat[["ks"]], tcc = Inf))
#'
#' r.mat.ind <- fitSynNLS(ratio2, t = tp, tcc = Inf, fitIndividual = TRUE)
#' plotCurve.comb(x = r.mat.ind, t = tp, tcc = Inf, curve = "syn")

fitSynNLS <- function(x, t, tcc=Inf,
                      A = NULL, B = NULL, ks = NULL,
                      fitIndividual = FALSE,
                      par.init = list(A=0.9, B=0.1, ks=0.04),
                      par.lower=c(A=0, B=0, ks=0),
                      par.upper=c(A=1, B=1, ks=10),
                      message = TRUE) {

  if (any(!names(par.init) %in% c("A", "B", "ks")) ||
      any(!names(par.lower) %in% c("A", "B", "ks")) ||
      any(!names(par.upper) %in% c("A", "B", "ks")))
    stop("Names of 'par' should be 'A', 'B' and 'ks'.")

  # swap A and B, rename ks to kd
  swp.A <- B
  swp.B <- A
  swp.par.init <- list(A = par.init$B, B = par.init$A, kd = par.init$ks)
  swp.par.lower <- c(A = par.lower[["B"]], B = par.lower[["A"]], kd = par.lower[["ks"]])
  swp.par.upper <- c(A = par.upper[["B"]], B = par.upper[["A"]], kd = par.upper[["ks"]])

  if (inherits(x, "vector") || is.vector(x)) {
    r <- fitDegNLS(x=x, t=t, tcc=tcc,
                   A = swp.A, B = swp.B, kd = ks,
                   par.init = swp.par.init,
                   par.lower = swp.par.lower,
                   par.upper = swp.par.upper,
                   message = TRUE, kd2ks = TRUE)
  } else if (inherits(x, "matrix")) {
    r <- fitDegNLS(x=x, t=t, tcc=tcc,
                   A = swp.A, B = swp.B, kd = ks,
                   fitIndividual = fitIndividual,
                   par.init = swp.par.init,
                   par.lower = swp.par.lower,
                   par.upper = swp.par.upper,
                   message = TRUE, kd2ks = TRUE)
  } else
    stop ("x should be either a vector or a matrix")
  r
}
