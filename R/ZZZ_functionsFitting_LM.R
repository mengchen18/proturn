#' # this file is not updated
#'
#' #' @title Fitting protein degradation curve using linear method (LM)
#' #'
#' #' 'lm' used to fit the curve, after log transform. A is not used because it
#' #' will be estimated by the function (the intecept)
#' #' @param x The ratio
#' #' @param t The time points (hours)
#' #' @param B The B
#' #' @param tcc The doubling time
#' #' @param intercept0 A logical to indicate if the intercept should be 0,
#' #' namely A = 1+B
#' #' @param outlier.rm A logical to indicate if outlier in the linear
#' #' regression should be removed. The outlier is identified by cook's
#' #' distance > 1. Up to two iterations could be applied.
#' #' @param message Logical to indicate if the message should be printed.
#' #' @author Chen Meng
#' #' @exportMethod
#'
#' setGeneric("fitDegLM", function(x, ...) {
#'   standardGeneric("fitDegLM")
#' })
#'
#' setMethod("fitDegLM", signature(x = "vector"),
#'           function(x, t, B=0, tcc=Inf, intercept0 = FALSE,
#'                    outlier.rm = FALSE, message = TRUE) {
#'
#'             stopifnot(B >= 0 & B <= 1)
#'             err.return <- structure(rep(NA, 4), names = c("A", "B", "kd", "r.sq"))
#'             xb <- x-B
#'             idx <- xb > 0
#'
#'             if (sum(idx) < 2) {
#'               if (message)
#'                 message("less than 2 values are valid for fitting, NA returned.")
#'               return(err.return)
#'             }
#'
#'             xb <- xb[idx]
#'             t <- t[idx]
#'             R <- log(xb)
#'             xr <- x[idx]
#'
#'             if (intercept0)
#'               form <- -R-t*(log(2)/tcc) ~ -1 + t else
#'                 form <- -R-t*(log(2)/tcc) ~ t
#'
#'             fit <- try(lm(form), silent = TRUE)
#'
#'             if (!inherits(fit, "lm")) {
#'               if (message)
#'                 message("Error in fitting model.")
#'               return(err.return)
#'             }
#'
#'             iter <- 1
#'             while (outlier.rm && iter <= 2) {
#'               cd <- cooks.distance(fit)
#'               icd <- cd < 1
#'               if (all(icd))
#'                 break
#'               iter <- iter+1
#'               if (message)
#'                 message(paste("outlier point",
#'                               paste(which(!icd), collapse = ","),
#'                               "is removed."))
#'               R <- R[icd]
#'               t <- t[icd]
#'               xr <- xr[icd]
#'               fit <- try(lm(form), silent = TRUE)
#'
#'               if (!inherits(fit, "lm")) {
#'                 if (message)
#'                   message("Error in fitting model.")
#'                 return(err.return)
#'               }
#'             }
#'
#'             res <- fit[[1]]
#'             r.sq <- summary(fit)$r.squared
#'             A <- ifelse(intercept0, 1+B, exp(-res["(Intercept)"]) + B)
#'             structure(c(A, B, res["t"], r.sq),
#'                       names = c("A", "B", "kd", "r.sq"),
#'                       x.used = xr, t.used = t)
#'           })
#'
#' setMethod("fitDegLM", signature(x = "matrix"),
#'           function(x, t, B=0, tcc=Inf, fitIndividual = FALSE,
#'                    intercept0 = FALSE, outlier.rm = FALSE, message = TRUE) {
#'
#'             stopifnot(ncol(x) == length(t))
#'
#'             nr <- nrow(x)
#'             xx <- c(t(x))
#'             tt <- rep(t, nr)
#'
#'             if (fitIndividual || outlier.rm) {
#'               r <- lapply(1:nrow(x), function(r) {
#'                 fitDegLM(x[r, ], t, B=B, tcc=tcc, intercept0 = intercept0,
#'                          outlier.rm = outlier.rm, message = message)
#'               })
#'               rr <- do.call("rbind", r)
#'               xx <- unlist(lapply(r, attr, "x.used"))
#'               tt <- unlist(sapply(r, attr, "t.used"))
#'             }
#'
#'             rcomb <- fitDegLM(xx, tt, B=B, tcc=tcc,
#'                               intercept0 = intercept0,
#'                               outlier.rm = outlier.rm,
#'                               message = message)
#'             structure(rcomb, individual = rr)
#'           })
#'
#'
#' # not used
#' # matchScale <- function(deg, syn, turnover.rm = FALSE, match.bottom = NA,
#' #                        torm.first = FALSE) {
#' #
#' #   nc.deg <- ncol(deg)
#' #   nc.syn <- ncol(syn)
#' #   if (nc.deg != nc.syn)
#' #     stop ("The number of columns of degn and syn should be the same.")
#' #   if (is.data.frame(deg))
#' #     deg <- as.matrix(deg, rownames.force = FALSE)
#' #   if (is.data.frame(syn))
#' #     syn <- as.matrix(syn, rownames.force = FALSE)
#' #   it <- intersect(rownames(deg), rownames(syn))
#' #   if (length(it) == 0)
#' #     stop("no intersect rownames between deg and syn. Check row names.")
#' #
#' #   if (is.na(match.bottom[1])) {
#' #     a <- log10(deg[it, 1])
#' #     b <- log10(syn[it, nc.syn])
#' #     fab <- !(is.infinite(a) | is.infinite(b) | is.na(a) | is.na(b))
#' #     a <- a[fab]
#' #     b <- b[fab]
#' #     seqa <- seq(min(a), max(a), length.out = 1000)
#' #     corpart.a <- sapply(seqa, function(i) {
#' #       idx <- a < i
#' #       cor(a[idx], b[idx])
#' #     })
#' #     seqb <- seq(min(b), max(b), length.out = 1000)
#' #     corpart.b <- sapply(seqb, function(i) {
#' #       idx <- b < i
#' #       cc <- cor(a[idx], b[idx])
#' #     })
#' #     cuta <- seqa[ max(which(corpart.a < 0.05)) ]
#' #     cutb <- seqb[ max(which(corpart.a < 0.05)) ]
#' #     match.bottom <- c(cuta, cutb)
#' #   } else if (length(match.bottom) == 1) {
#' #     match.bottom <- rep(match.bottom, 2)
#' #   }
#' #
#' #   on.exit(par)
#' #   ### plot 1
#' #   op <- par(no.readonly = TRUE)
#' #   on.exit(par(op), add = TRUE)
#' #   layout(matrix(1:3, 1, 3))
#' #   plot(a, b, xlab = "Intensity of deg time point 0",
#' #        ylab = paste("Intensity of Syn time point", nc.syn))
#' #   abline(h = cutb, v = cuta)
#' #   mtext(3, at = cuta, text = signif(cuta, 3))
#' #   mtext(4, at = cutb, text = signif(cutb, 3))
#' #
#' #   # NA, 0 should be excluded
#' #   iit <- deg[it, 1] > match.bottom[1]
#' #   iit <- iit & syn[it, nc.syn] > match.bottom[2]
#' #   iit <- iit & !is.na(deg[it, 1])
#' #   iit <- iit & !is.na(syn[it, nc.syn])
#' #   it <- it[iit]
#' #
#' #
#' #   # correct turnover
#' #   if (turnover.rm && torm.first) {
#' #     dsum <- colSums(deg[it, ], na.rm = TRUE)
#' #     dsum <- dsum/dsum[1]
#' #     ssum <- colSums(syn[it, ], na.rm = TRUE)
#' #     ssum <- ssum/ssum[nc.syn]
#' #
#' #     deg.off <- sapply(ssum, function(x) x*deg[, nc.deg])
#' #     syn.off <- sapply(dsum, function(x) x*syn[, 1])
#' #
#' #     deg <- deg - deg.off
#' #     syn <- syn - syn.off
#' #     deg[deg < 0] <- 0
#' #     syn[syn < 0] <- 0
#' #   }
#' #
#' #   # correct each row
#' #   mc <- (deg[it, 1] + syn[it, nc.syn])/2
#' #   deg[it, ] <- deg[it, ] * mc / deg[it, 1]
#' #   syn[it, ] <- syn[it, ] * mc / syn[it, nc.syn]
#' #
#' #   deg[is.infinite(deg)] <- NaN
#' #   syn[is.infinite(syn)] <- NaN
#' #
#' #   # correct summed intensity
#' #   sfsep <- rbind(colSums(deg[it, ], na.rm = TRUE),
#' #                  colSums(syn[it, ], na.rm = TRUE))
#' #   sf1 <- colSums(sfsep)
#' #   sf1 <- sfsep/mean(sf1)
#' #   sf <- colSums(sf1)
#' #   ### plot 2
#' #   vv <- barplot(sf1, ylab = "scaling factor")
#' #   abline(h = 1)
#' #   mtext(side = 3, at = vv, text = signif(sf, digits = 3), las = 2)
#' #   deg <- sweep(deg, 2, sf, "/")
#' #   syn <- sweep(syn, 2, sf, "/")
#' #
#' #   # correct turnover
#' #   if (turnover.rm && !torm.first) {
#' #     dsum <- colSums(deg[it, ], na.rm = TRUE)
#' #     dsum <- dsum/dsum[1]
#' #     ssum <- colSums(syn[it, ], na.rm = TRUE)
#' #     ssum <- ssum/ssum[nc.syn]
#' #
#' #     deg.off <- sapply(ssum, function(x) x*deg[, nc.deg])
#' #     syn.off <- sapply(dsum, function(x) x*syn[, 1])
#' #
#' #     deg <- deg - deg.off
#' #     syn <- syn - syn.off
#' #     deg[deg < 0] <- 0
#' #     syn[syn < 0] <- 0
#' #   }
#' #
#' #   ## plot 3
#' #   mm <- rbind(colSums(deg[it, ], na.rm = TRUE), colSums(syn[it, ], na.rm = TRUE))
#' #   barplot(mm, ylab = "Summed intensty (intersect) after removing turnover.")
#' #
#' #   ## calculate ratio
#' #   deg.ratio <- deg/deg[, 1]
#' #   deg.ratio[is.infinite(deg.ratio)] <- NA
#' #   syn.ratio <- syn/syn[, nc.syn]
#' #   syn.ratio[is.infinite(syn.ratio)] <- NA
#' #
#' #   list(deg=deg, syn=syn, deg.ratio = deg.ratio, syn.ratio = syn.ratio,
#' #        mappedid = it, mc = structure(mc, names = it), scale.factor = sf)
#' # }
#'
