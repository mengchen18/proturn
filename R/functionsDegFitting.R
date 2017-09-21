#' @title Fitting degration curve using NLS algorithm
#' @description fitting degration curve using NLS algorithm
#' @param x a numeric matrix or vector
#' @param ... other arguments
#' @param t The time point (hours)
#' @param tcc The doubling time of cells. By default this value is Inf, which means the cells are in steady state
#' @param A optinal argument for fixed A, if this argument is given, "A" won't be optimized
#' @param B optional argument for fixed B, if this argument is given, "B" won't be optimized
#' @param kd optional argument for fixed kd, if this argument is given, "ks" won't be optimized
#' @param par.init The initial values of parameters to be optimized, it should be list of three elements
#'   names as "A", "B" and "kd".
#' @param par.lower The lower boundary of parameters to be optimized, it should be a numeric values with
#'   length 3 and named as "A", "B" and "kd".
#' @param par.upper The upper boundary of parameters to be optimized, it should be a numeric values with
#'   length 3 and named as "A", "B" and "kd".
#' @param message A logical value to indicated if any message should be printed
#' @param fitIndividual A logical value, whether each individual row should also be fitted.
#'   Only used when x is an object of class \code{matrix}
#' @param kd2ks Should not be changed by user. A logical value, whether should be transformed to
#'   fit synthesis curve.
#' @return
#'  a vector of optimized parameters, including A, B, kd, confidence intervals (2.5% and 97.5%),
#'  mean square error and r-square values.
#'  In addition, if individual rows are fitted, the object also contains an attribute stores
#'  parameters fitted on each individual row.
#' @export
#' @import methods stats
#' @examples
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
#'
#'
setGeneric("fitDegNLS", function(x, ...) {
  standardGeneric("fitDegNLS")
})

#' @describeIn fitDegNLS Fitting degradation curve given a vector
setMethod(
  "fitDegNLS", signature(x = "vector"),
  function(x, t, tcc=Inf, A = NULL, B = NULL, kd = NULL,
           par.init = list(A=0.9, B=0.1, kd=0.04),
           par.lower=c(A=0, B=0, kd=0),
           par.upper=c(A=1, B=1, kd=10),
           message = TRUE, kd2ks = FALSE) {
    
    call <- match.call()
    if (any(!names(par.init) %in% c("A", "B", "kd")) ||
        any(!names(par.lower) %in% c("A", "B", "kd")) ||
        any(!names(par.upper) %in% c("A", "B", "kd")))
      stop("Names of 'par' should be 'A', 'B' and 'kd'.")

    if (!is.null(A))
      par.init$A <- par.lower["A"] <- par.upper["A"] <- A
    if (!is.null(B))
      par.init$B <- par.lower["B"] <- par.upper["B"] <- B
    if (!is.null(kd))
      par.init$kd <- par.lower["kd"] <- par.upper["kd"] <- kd
    k <- ifelse(kd2ks, "ks", "kd")
    names <- c("A", "B", k, "ci025", "ci975", "mse", "rsq")
    err.return <- structure(rep(NA, 7), names = names)

    # define formula
    form <- x ~ (A - B) * exp(-(kd+log(2)/tcc) * t) + B

    # remove NA ratios, ratio higher than 1, less than 0 and corresponding time points
    to <- t
    xo <- x
    idx <- !is.na(x) # & x >= 0 # & x <= 1
    if (sum(idx) <= 2) return( err.return )
    x <- x[idx]
    t <- t[idx]
    alg <- NA

    # fitting first try the 'port' algorithm
    fit <- try(nls(form, start = par.init, algorithm="port",
                   lower = par.lower, upper = par.upper,
                   control = list(warnOnly = TRUE)),
               silent = TRUE)
    conv <- fit$convInfo$isConv
    if (!inherits(fit, "nls")) {
      res <- err.return
    } else {
      use.unconv.port <- FALSE

      if (!conv) {
        fit.pl <- try(nls(form, start = par.init, algorithm="plinear"), silent = TRUE)
        if (inherits(fit.pl, "nls")) {
          if (message)
            message("'plinear' method used, upper and lower limits of parameters are ignored.")
          res <- coef(fit.pl)[-4]
          alg <- "plinear"
        } else
          use.unconv.port <- TRUE
      }
      if (conv || use.unconv.port) {
        res <- coef(fit)
        alg <- paste("port", fit$convInfo$stopMessage)
      }
      # }
      #
      #
      # if (inherits(fit, "nls")) {
      mse <- try(mean(residuals(fit)^2), silent = TRUE)
      if (!is.numeric(mse)) {
        mse <- NA
        rsq <- NA
      }
      rsq <- 1 - var(residuals(fit))/var(x)
      ci <- try(confint(fit, "kd"), silent = TRUE)
      if (inherits(ci, "try-error"))
        ci <- c(NA, NA)
      res <- c(res, ci, mse, rsq)
    }
    out <- structure(res, names = names, algorithm = alg)
    if (kd2ks) {
      out.before.swp <- out
      out[["A"]] <- out.before.swp[["B"]]
      out[["B"]] <- out.before.swp[["A"]]
    }
    structure(out, x = x, type = ifelse(kd2ks, "Synthesis", "Degradation"),
              call = call)
  })

#' @describeIn fitDegNLS Fitting degradation curve given a matrix
setMethod(
  "fitDegNLS", signature(x = "matrix"),
  function(x, t, tcc=Inf,
           A = NULL, B = NULL, kd = NULL,
           fitIndividual = FALSE,
           par.init = list(A=0.9, B=0.1, kd=0.04),
           par.lower=c(A=0, B=0, kd=0),
           par.upper=c(A=1, B=1, kd=10),
           message = TRUE, kd2ks = FALSE) {

    stopifnot(ncol(x) == length(t))
    fitIndInput <- fitIndividual
    nr <- nrow(x)
    if (nr == 1)
      fitIndividual <- FALSE

    A1 <- length(A) == 1
    B1 <- length(B) == 1
    k1 <- length(kd) == 1
    r <- NA
    if (fitIndividual) {
      r <- sapply(1:nrow(x), function(i) {
        iA <- iB <- ik <- i
        if (A1) iA <- 1
        if (B1) iB <- 1
        if (k1) ik <- 1
        fitDegNLS(x[i, ], t=t, tcc=tcc,
                  A = A[iA], B = B[iB], kd = kd[ik],
                  par.init = par.init,
                  par.lower = par.lower,
                  par.upper = par.upper,
                  kd2ks = kd2ks)
      })
      r <- t(r)
      k <- ifelse(kd2ks, "ks", "kd")
      colnames(r) <- c("A", "B", k, "ci025", "ci975",  "mse", "rsq")
      rownames(r) <- rownames(x)
    }

    xx <- c(t(x))
    tt <- rep(t, nr)
    mA <- mB <- mk <- NULL
    if (!is.null(A)) mA <- mean(A, na.rm = TRUE)
    if (!is.null(B)) mB <- mean(B, na.rm = TRUE)
    if (!is.null(kd)) mk <- mean(kd, na.rm = TRUE)
    rcomb <- fitDegNLS(xx, tt, tcc=tcc,
                       A = mA, B = mB, kd = mk,
                       par.init = par.init,
                       par.lower = par.lower,
                       par.upper = par.upper,
                       message = message,
                       kd2ks = kd2ks)


    if (nr == 1)
      r <- matrix(rcomb, nrow = 1, dimnames = list(rownames(x), names(rcomb)))
      

    structure(rcomb, individual = r, inputmatrix = x,
              type = ifelse(kd2ks, "Synthesis", "Degradation"),
              call = call)
  })
