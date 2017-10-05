#' @title Fitting protein degradation/synthesis curves using nonlinear 
#'   least square methods (NLS)
#' @description Main function to be called by users to fit
#'   protein degradation or synthesis curves using 
#'   nonlinear least square methods (NLS).
#' @param x A numeric matrix where rows are variables (e.g. peptides or proteins)
#'   and columns are different time points.
#' @param f A factor or vector, its length should be the same as the number of 
#'   rows in \code{x}. It is possible to fit a single model using multiple lines, this is 
#'   useful, for example, when multiple lines are from the peptide or protein. This
#'   is implemented by multiple rows in \code{x} corresponds to the same with in \code{f}.
#' @param t The time point (hours)
#' @param tcc The doubling time of cells. By default this value is Inf, 
#'   which means the cells are in steady state. 
#' @param type which curve should be fitted, should be either "deg" (degradation curve) or "syn" (synthesis curve)
#' @param A optinal argument for fixed A, if this argument is given, "A" won't be optimized
#' @param B optional argument for fixed B, if this argument is given, "B" won't be optimized
#' @param par.init The initial values of parameters to be optimized, it should be list of three elements
#'   names as "A", "B" and "k".
#' @param par.lower The lower boundary of parameters to be optimized, it should be a numeric values with
#'   length 3 and named as "A", "B" and "k".
#' @param par.upper The upper boundary of parameters to be optimized, it should be a numeric values with
#'   length 3 and named as "A", "B" and "k".
#' @param ncore the number of cores to be used, passed to \code{mclapply}
#' @details The function fit the following models:
#' for degradation
#' for synthesis
#' Using nls algorithm.
#'
#' @return A list of three elmenets
#'   \itemize{
#'     \item mat
#'     \item list
#'     \item type could be either combine, individual
#'   }
#'
#' @importFrom parallel mclapply
#' @export
#' @examples
#' tp <- c(0, 1, 2, 4, 8, 16, 32, 64)
#' A <- runif(5, m                                         
  n 
#' B <- runif(5, min = 0.1, max = 0.15)
#' kd <- runif(5, min = 0.05, max = 0.1)
#' ds <- mapply(function(A, B, kd) {
#'   degCurve(A=A, B=B, kd=kd, tcc=Inf, t = tp) + rnorm(length(tp), sd = 0.03)
#' }, A, B, kd)
#' ds <- t(ds)
#' ds <- rbind(ds, runif(8, 0, 1))
#'
#' rf <- c("p1", "p2", "p2", "p3", "p3", "p3")
#'
#' ss <- fitNLSModels(x = ds, f = rf, t = tp, tcc = Inf, type = "deg",
#'                    par.init = list(A=0.8, B=0.2, kd=0.04),
#'                    par.lower=c(A=0, B=0, kd=0),
#'                    par.upper=c(A=1, B=1, kd=10))
#'
#' plotCurve.comb(x = ss$list$p1, t = tp, tcc = Inf, curve = "deg")
#' plotCurve.comb(x = ss$list$p2, t = tp, tcc = Inf, curve = "deg")
#' plotCurve.comb(x = ss$list$p3, t = tp, tcc = Inf, curve = "deg")
#'
#'
#' ds2 <- mapply(function(A, B, kd) {
#'   synCurve(A=A, B=B, kd=kd, tcc=Inf, t = tp) + rnorm(length(tp), sd = 0.03)
#' }, A, B, kd)
#' ds2 <- t(ds2)
#' ds2 <- rbind(ds2, runif(8, 0, 1))
#'
#'
#' ss2 <- fitNLSModels(x = ds2, f = rf, t = tp, tcc = Inf, type = "syn",
#'                     par.init = list(A=0.8, B=0.2, ks=0.04),
#'                     par.lower=c(A=0, B=0, ks=0),
#'                     par.upper=c(A=1, B=1, ks=10))
#'
#' plotCurve.comb(x = ss2$list$p1, t = tp, tcc = Inf, curve = "syn")
#' plotCurve.comb(x = ss2$list$p2, t = tp, tcc = Inf, curve = "syn")
#' plotCurve.comb(x = ss2$list$p3, t = tp, tcc = Inf, curve = "syn")


fitNLSModels <- function(x, f, t, tcc=Inf, type = c("deg", "syn"),
                         A = NULL, B = NULL,
                         par.init = list(A=0.9, B=0.1, k=0.04),
                         par.lower=c(A=0, B=0, k=0),
                         par.upper=c(A=1, B=1, k=10),
                         ncore = 1) {
  f <- as.character(f)
  type <- match.arg(type, c("deg", "syn"))
  if (type == "deg") {
    names(par.init)[names(par.init) == "k"] <- "kd"
    names(par.lower)[names(par.lower) == "k"] <- "kd"
    names(par.upper)[names(par.upper) == "k"] <- "kd"
  } else {
    names(par.init)[names(par.init) == "k"] <- "ks"
    names(par.lower)[names(par.lower) == "k"] <- "ks"
    names(par.upper)[names(par.upper) == "k"] <- "ks"
  }
  
  if (is.null(rownames(x)))
    rownames(x) <- paste("X", 1:nrow(x), sep = "")
  type <- match.arg(type[1], choices = c("deg", "syn"))
  fun <- switch (type,
                 "deg" = fitDegNLS,
                 "syn" = fitSynNLS)
  
  diffB <- length(B) == nrow(x) && !is.null(B)
  diffA <- length(A) == nrow(x) && !is.null(A)
  
  ll <- mclapply(unique(f), function(item) {
    i <- f == item
    m <- x[i, , drop = FALSE]
    
    bb <- B[1]
    aa <- A[1]
    if (diffB) {
      bb <- B[i] 
      if (is.na(bb) || !is.numeric(bb)) {
        bb <- NULL
        cat("NA is given to B, use NULL instead.\n")
      } 
    }
    if (diffA) {
      aa <- A[i]
      if (is.na(aa) || !is.numeric(aa)) {
        aa <- NULL
        cat("NA is given to A, use NULL instead.\n")
      }
    }
    
    fun(x = m, t = t, tcc=tcc, A = aa, B = bb,
        par.init = par.init,
        par.lower = par.lower,
        par.upper = par.upper,
        message = TRUE,
        fitIndividual = TRUE)
  }, mc.cores = ncore)
  names(ll) <- unique(f)
  #
  fit.mat <- do.call(rbind, lapply(ll, attr, "individual"))
  fit.mat <- fit.mat[rownames(x), ]
  fit.mat <- cbind(f, fit.mat)
  
  type <- "individual"
  comb.fit <- max(table(f)) > 1
  if (comb.fit) {
    l3 <- lapply(ll, function(x) {
      ind <- attr(x, "individual")
      nr <- nrow(ind)
      if (is.null(nr)) nr <- 1
      matrix(rep(x, nr), nrow = nr, byrow = TRUE,
             dimnames = list(rownames(ind), 
                             paste0("comb.", names(x))))
    })
    comb.fit.mat <- do.call(rbind, l3)
    comb.fit.mat <- comb.fit.mat[rownames(x), ]
    fit.mat <- cbind(fit.mat, comb.fit.mat)
    type <- "combined"
  }
  
  fit.mat <- data.frame(fit.mat, stringsAsFactors = FALSE)
  fit.mat[-1] <- lapply(fit.mat[-1], as.numeric)
  colnames(fit.mat)[1] <- "collapsed.factor"
  
  list(mat = fit.mat, list = ll, type = type)
}



