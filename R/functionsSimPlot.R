#' @title Calculating ratio of proteins turnover using theoretical formula
#' @description Generate ratio of proteins turnover using theoretical formula. This used for control/unit test.
#' @param A model parameter A (amplitude), see details
#' @param B model parameter B (offset), see details.
#' @param kd model parameter kd (degradation/synthesis constant), see details.
#' @param tcc model parameter tcc (doubling time), see details.
#' @param t The time points (hours)
#' @details Generate
#'  withMathJax("$$f(t)=(A-B)\\cdot e^{-(k_d+\\frac{ln2}{t_{cc}})\\cdot t}+B$$")
#'  "The model fitted for the synthesis curve is as:"
#'  withMathJax("$$f(t)=(B-A)\\cdot e^{-(k_s+\\frac{ln2}{t_{cc}})\\cdot t}+A$$")
#'  "where t is the time points (given in hours)
#' @examples
#'  times <- 2^(0:6)
#'  drat <- degCurve(A = 1, B = 0, kd = 0.05, tcc = Inf, times)
#'  plot(times, drat)
#' @return A vector of ratios of protein turned over
#' @keywords internal
#' @export
degCurve <- function(A, B, kd, tcc, t) {
  (A - B) * exp(-(kd+log(2)/tcc) * t) + B
}

#' @rdname degCurve
#' @export
#' @keywords internal
#' @examples
#'  times <- 2^(0:6)
#'  srat <- synCurve(A = 1, B = 0, kd = 0.05, tcc = Inf, times)
#'  plot(times, srat)
synCurve <- function(A, B, kd, tcc, t) {
  (B - A) * exp(-(kd+log(2)/tcc) * t) + A
}


#' @title draw horizontal error bars
#' @description draw horizontal error bar, not called by users
#' @param x a vector of length 2
#' @param y a single value or a numerical vector specify the position on y-axis
#' @param lwd line width
#' @param length the length of vertical bars at the end of arrow bars
#' @param col color
#' @param lty the lty
#' @param pch shape for points of mean value
#' @import graphics
#' @return no value to be returned
#' @examples
#' # plot(0, 0)
#' # ebar(c(-0.1, 0.1), c(0), lty = 2, lwd = 2, col = "red")
#' # ebar(c(-0.1, 0.1), c(0.1, 0, -0.4), lty = 2)
#'
ebar <- function(x, y, lwd = 2, col = 1, lty = 1, length = 0.1, pch = 18) {
  if (is.list(x))
    x <- unlist(x)
  x1 <- rep(x, length(y))
  y <- rep(y, each = 2)
  mp <- rep(mean(x), length(y))
  arrows(mp, y, x1, y, code = 2, angle = 90, length = length, col = col, lty = lty, lwd = lwd)
  points(mp, y, pch = pch, col = col)
}

#' @title plot protein degration/synthesis curve
#' @description plot protein degration/synthesis curve with origin data point and fitted curve
#' @param x The data point, ratio, if \code{lineOnly=TRUE}, this argument is ignored
#' @param t time points (hours)
#' @param tcc doubling time of cell lines
#' @param A model parameter A (amplitude), see details
#' @param B model parameter B (offset), see details.
#' @param k model parameter kd (degradation constant) or sd (synthesis constant), see details.
#' @param add if the plot should be added on top of another plot
#' @param lineOnly logical, whether only draw the fitted line, ie. the input data point is ignored
#' @param curve which curve want to draw, should be either "degradation" or "synthesis"
#' @param pch passed to \code{plot}
#' @param lty passed to \code{line}
#' @param ylim passed to \code{plot}
#' @param err.x draw error bar
#' @param err.y draw error bar
#' @param err.lwd the line width of error bars
#' @param main the main title of plot
#' @param col The color of line
#' @return no value to be returned
#' @examples
#'  times <- 2^(0:6)
#'  drat <- degCurve(A = 1, B = 0, kd = 0.05, tcc = Inf, times)
#'  drat <- drat + rnorm(length(times), sd = 0.2)
#'  plot(times, drat)
#'  plotCurve(drat, t = times, tcc = Inf, A = 1, B = 0, k = 0.05, add=FALSE, col="red",
#'        lineOnly = FALSE, curve = "degradation", ylim = NULL,
#'        err.x = c(11, 20), err.y = 0.5, main = "toy curve fitting")
#' @export
#'
plotCurve <- function(x, t, tcc, A, B, k, add=TRUE, col="red", lineOnly = FALSE,
                      curve = c("degradation", "synthesis")[1], pch=20, lty = 1,
                      ylim = NULL, err.x = NULL, err.y = NULL, err.lwd = 2, main = "") {

  curve <- match.arg(curve[1], c("degradation", "synthesis"))
  if (is.null(ylim))
    ylim <- c(0, max(1, max(x)))
  fo <- switch (curve,
    "degradation" = degCurve,
    "synthesis"  = synCurve
  )
  if (!lineOnly) {
    if (add)
      points(t, x, pch = pch, col = col) else
        plot(t, x, xlab="time points (Hours)", ylab="Ratio", pch = pch, col = col, ylim = ylim, main = main)
  }
  tp <- 0:max(t)
  y <- fo(A, B, k, tcc, tp)
  lines(tp, y, col=col, lwd=2, lty = lty)

  if (!is.null(err.x) && !is.null(err.y)) {
    ebar(x = err.x, y = err.y, lwd = err.lwd, col = col, lty = lty)
  }
}


#' @title plot a single fit of degradation or synthesis curve
#' @description plot a single fit of degradation or synthesis curve
#' @param x an object returned by \code{\link{fitSynNLS}} or \code{\link{fitDegNLS}} with vector input
#' @param t time points
#' @param tcc cell doubling time
#' @param col color, passed to plot
#' @param curve what type of curve fit, if is not null, should be either "degradation" or "synthesis"
#' @param pch passed to plot
#' @param lty passed to plot
#' @param err.lwd the line width of error bars
#' @param main main title, passed to plot
#' @note This function won't work if input of \code{\link{fitSynNLS}} or \code{\link{fitDegNLS}} is a one row matrix
#' @return a plot generated, returns nothing
#' @examples 
#'   tp <- c(0, 1, 2, 4, 8, 16, 32, 64)
#'   ratios <- degCurve(A=0.85, B = 0.1, kd=0.5, tcc=Inf, t = tp) + rnorm(length(tp), sd = 0.05)
#'   # vector input
#'   r <- fitDegNLS(ratios, t = tp, tcc = Inf)
#'   plotCurve.one(r, tp)
#' @export
#' 
plotCurve.one <- function(x, t, tcc = Inf, col="red", curve = NULL, 
                          pch=20, lty = 1, err.lwd = 2, main = "") {
  
  if (is.null(curve))
    curve <- tolower(attr(x, "type"))
  curve <- match.arg(curve[1], c("degradation", "synthesis"))
  
  k <- paste0("k", substr(curve, 1, 1))
  plotCurve(attr(x, "x"), t = t, tcc = tcc, A = x["A"], B = x["B"], k = x[k], 
            add=FALSE, col=col,
            lineOnly = FALSE, curve = curve, 
            err.x = log(2)/x[c("ci025", "ci975")], 
            err.y = degCurve(A = x[["A"]], B = x[["B"]], kd = x[[k]],
                             t = log(2)/ x[[k]], tcc = Inf), 
            main = main)
}


#' plot curve from combined fitting
#' @param x - an element of list returned by fitDegNLS or fitSynNLS
#' @param t time points (hours)
#' @param leg.vec a character vector for legend, the name of the vector should be the same as in x
#' @param curve which curve want to draw. If it's not \code{NULL}, 
#'   should be either "degradation" or "synthesis"
#' @param tcc doubling time of cell lines
#' @param add if the plot should be added on top of another plot
#' @param pch passed to \code{plot}
#' @param lty passed to \code{line}
#' @param legend logical value, whether legend should be generated
#' @param err a logical value, whether error bar should be plotted
#' @param main the main title of plot
#' @param leg.cex the cex of points in legend
#' @export
#' @return no value to be returned 
#' @examples
#' # see \code{fitSynNLS}, \code{fitDegNLS}, \code{fitNLSModels}
#' # simulating data
#' tp <- c(0, 1, 2, 4, 8, 16, 32, 64)
#' ratios <- degCurve(A=0.85, B = 0.1, kd=0.5, tcc=Inf, t = tp) + rnorm(length(tp), sd = 0.05)
#'
#' #' vector input
#' r <- fitDegNLS(ratios, t = tp, tcc = Inf)
#' plotCurve(ratios, tp, tcc = Inf, A = r[["A"]], B = r[["B"]], k = r[["kd"]],
#'           add = FALSE, curve = "deg", err.x = log(2)/r[c("ci025", "ci975")],
#'           err.y = degCurve(A = r[["A"]], B = r[["B"]], kd = r[["kd"]], t = log(2)/ r[["kd"]], tcc = Inf))
#'
#' #' matrix input, fit a single model
#' ratio2 <- rbind(p1 = ratios + rnorm(length(ratios), sd = 0.4),
#'                 p2 = ratios + rnorm(length(ratios), sd = 0.4))
#' r.mat <- fitDegNLS(ratio2, t = tp, tcc = Inf)
#' plotCurve(x = ratio2, t = rep(tp, nrow(ratio2)),
#'           tcc = Inf, A = r.mat[["A"]], B = r.mat[["B"]], k = r.mat[["kd"]],
#'           add = FALSE, curve = "deg", err.x = log(2)/r.mat[c("ci025", "ci975")],
#'           err.y = degCurve(A = r.mat[["A"]], B = r.mat[["B"]], kd = r.mat[["kd"]],
#'                            t = log(2)/ r.mat[["kd"]], tcc = Inf))
#'
#' #' matrix input, fit a single model, in addition, each individual row should also be fitted
#' r.mat.ind <- fitDegNLS(ratio2, t = tp, tcc = Inf, fitIndividual = TRUE)
#' plotCurve.comb(x = r.mat.ind, t = tp, tcc = Inf,
#'                leg.vec = c(p1="peptide 1", p2 = "peptide 2"), curve = "deg")

plotCurve.comb <- function(x, t, leg.vec = NULL, curve = NULL,
                           tcc = Inf, add = FALSE, pch = 20, lty = 1, legend = TRUE,
                           err = TRUE, main = "", leg.cex = 1) {
  if (is.null(x) || all(is.na(x))) {
    cat("no available values\n")
    return(invisible())
  }
  if (is.null(curve))
    curve <- tolower(attr(x, "type"))
  curve <- match.arg(curve[1], c("degradation", "synthesis"))
  fo <- switch (curve,
                "degradation" = degCurve,
                "synthesis"  = synCurve)

  curve <- match.arg(curve[1],  c("degradation", "synthesis"))
  k <- ifelse(curve == "degradation", "kd", "ks")
  ind <- attr(x, "individual")
  m <- attr(x, "inputmatrix")
  nr <- nrow(m)
  err.x <- NULL
  err.y <- NULL
  if (err) {
    err.x <- log(2)/x[c("ci025", "ci975")]
    err.y <- fo(A = x["A"], B = x["B"], x[k], tcc = tcc, t = log(2)/x[k])
  }

  for (i in 1:nrow(m))
    plotCurve(m[i, ], t = t, tcc = tcc, A = ind[i, "A"], B = ind[i, "B"], ylim =c(0, max(1, max(m, na.rm = TRUE))),
              k = ind[i, k], add = i!=1 || add,  main = main,
              col = i+1, curve = curve, pch = pch, lty = lty)

  plotCurve(NA, t = t, tcc = tcc, A = x["A"], B = x["B"], k = x[k],
            col = "black", lineOnly = TRUE, curve = curve, lty = lty,
            err.x = err.x, err.y = err.y)
  if (legend) {
    if (is.null(leg.vec))
      leg <- c(rownames(m), "combined") else
        leg <- c(leg.vec[rownames(m)], "combined")

    legend("right", legend = leg, col = c((1:nr)+1, 1), cex = leg.cex,
           pch = rep(c(20, NA), c(nr, 1)), lty = lty, bty = "n")
  }
}
