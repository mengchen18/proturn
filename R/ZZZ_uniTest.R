# library(RUnit)
#
# library(parallel)
#
# source("/home/chen/CloudChen/Projects/ProteinTurnOverRate/MSPD/R/functionsFitting.R")
# source("/home/chen/CloudChen/Projects/ProteinTurnOverRate/MSPD/R/functionsSimPlot.R")
#
# # tp = c(0, 1, 2, 4, 8, 16, 32, 64)
# # rd <- degCurve(A=0.6, B = 0.1, kd=0.5, tcc=Inf, t = tp)
# # rs <- synCurve(A=0.6, B = 0.1, kd=0.5, tcc=Inf, t = tp)
# #
# # plot(NA, xlim = c(0, 64), ylim = c(0, 1))
# # points(tp, rd, lty = 1)
# # lines(tp, rd, lty = 1)
# # points(tp, rs)
# # lines(tp, rs, lty = 1)
# #
# # test.fitDegNLS()
# # test.fitDegNLS()
# #
# #
# # ###
# # tp = c(0, 1, 2, 4, 8, 16, 32, 64)
# # ratios <- degCurve(A=0.85, B = 0.1, kd=0.5, tcc=12, t = tp)
# # r <- fitDegNLS(ratios, t = tp, tcc = Inf)
# # plotCurve(ratios,  t = tp, tcc = Inf, A = r[["A"]], B = r[["B"]], k = r[["kd"]], curve = "deg", add = FALSE)
# #
# #
# # tp = c(0, 1, 2, 4, 8, 16, 32, 64)
# # ratios <- synCurve(A=0.85, B = 0.1, kd=0.5, tcc=Inf, t = tp)
# # r <- fitSynNLS(ratios, t = tp, tcc = Inf)
# # r
# # plotCurve(ratios,  t = tp, tcc = Inf, A = r[["A"]], B = r[["B"]], k = r[["ks"]], curve = "syn", add = TRUE, pch = 1)
#
#
# ## ========================= test functions ==============================
# test.fitDegNLS <- function() {
#   tp = c(0, 1, 2, 4, 8, 16, 32, 64)
#   ratios <- degCurve(A=0.85, B = 0.1, kd=0.5, tcc=Inf, t = tp) + rnorm(length(tp), sd = 0.05)
#   r <- fitDegNLS(ratios, t = tp, tcc = Inf)
#   r.mat <- fitDegNLS(matrix(ratios, nrow = 1), t = tp, tcc = Inf)
#   r.mat.ind <- fitDegNLS(matrix(ratios, nrow = 1, dimnames = list("A", NULL)),
#                          t = tp, tcc = Inf, fitIndividual = TRUE)
#   # checkEquals(r[1:3] - c(A = 0.85, B = 0.1, kd = 0.5))
#   # checkEquals(r.mat[1:3], c(A = 0.85, B = 0.1, kd = 0.5))
#   checkEquals(attr(r.mat, "individual"), NA)
#   # checkEquals(r.mat.ind[1:3], c(A = 0.85, B = 0.1, kd = 0.5))
#   checkEquals(class(attr(r.mat.ind, "individual")), "matrix")
#   checkEquals(rownames(attr(r.mat.ind, "individual")), "A")
#   checkEqualsNumeric(attr(r.mat.ind, "individual")[1, ], r.mat.ind)
#   at <- attributes(r.mat.ind)
#   checkEquals(at$names, c("A", "B", "kd", "ci025", "ci975", "mse", "rsq"))
#   checkEquals(dim(at$individual), c(1, 7))
#   checkEquals(dim(at$inputmatrix), c(1, 8))
#
#   plotCurve(ratios, tp, tcc = Inf, A = r[["A"]], B = r[["B"]], k = r[["kd"]],
#             add = FALSE, curve = "deg", err.x = log(2)/r[c("ci025", "ci975")],
#             err.y = degCurve(A = r[["A"]], B = r[["B"]], kd = r[["kd"]], t = log(2)/ r[["kd"]], tcc = Inf))
#   plotCurve.comb(x = r.mat.ind, t = tp, tcc = Inf, leg.vec = c(A="A"), curve = "deg")
# }
#
# test.fitSynNLS <- function() {
#   tp = c(0, 1, 2, 4, 8, 16, 32, 64)
#   ratios <- synCurve(A=0.85, B = 0.1, kd=0.5, tcc=Inf, t = tp) + rnorm(length(tp), sd = 0.1)
#   r <- fitSynNLS(ratios, t = tp, tcc = Inf)
#   r.mat <- fitSynNLS(matrix(ratios, nrow = 1), t = tp, tcc = Inf)
#   r.mat.ind <- fitSynNLS(matrix(ratios, nrow = 1, dimnames = list("A", NULL)),
#                          t = tp, tcc = Inf, fitIndividual = TRUE)
#   # checkEquals(r[1:3], c(A = 0.85, B = 0.1, ks = 0.5))
#   # checkEquals(r.mat[1:3], c(A = 0.85, B = 0.1, ks = 0.5))
#   checkEquals(attr(r.mat, "individual"), NA)
#   # checkEquals(r.mat.ind[1:3], c(A = 0.85, B = 0.1, ks = 0.5))
#   checkEquals(class(attr(r.mat.ind, "individual")), "matrix")
#   checkEquals(rownames(attr(r.mat.ind, "individual")), "A")
#   checkEqualsNumeric(attr(r.mat.ind, "individual")[1, ], r.mat.ind)
#   at <- attributes(r.mat.ind)
#   checkEquals(at$names, c("A", "B", "ks", "ci025", "ci975", "mse", "rsq"))
#   checkEquals(dim(at$individual), c(1, 7))
#   checkEquals(dim(at$inputmatrix), c(1, 8))
#
#   plotCurve(ratios, tp, tcc = Inf, A = r[["A"]], B = r[["B"]], k = r[["ks"]],
#             add = FALSE, curve = "syn", err.x = log(2)/r[c("ci025", "ci975")],
#             err.y = degCurve(A = r[["A"]], B = r[["B"]], kd = r[["ks"]], t = log(2)/ r[["ks"]], tcc = Inf))
#
#   plotCurve.comb(x = r.mat.ind, t = tp, tcc = Inf, leg.vec = c(A="A"), curve = "syn")
# }
#
#
#
# source("/media/msdata5/users_files/Chen/Projects/ProteinTurnOverRate/MSPD/R/functionsFitting.R")
# source("/media/msdata5/users_files/Chen/Projects/ProteinTurnOverRate/MSPD/R/functionsSimPlot.R")
#
# library(RUnit)
# test.fitNLSModels()
#
# test.fitNLSModels <- function() {
#
#   tp <- c(0, 1, 2, 4, 8, 16, 32, 64)
#   A <- runif(5, min = 0.75, max = 0.95)
#   B <- runif(5, min = 0.1, max = 0.15)
#   kd <- runif(5, min = 0.05, max = 0.1)
#   ds <- mapply(function(A, B, kd) {
#     degCurve(A=A, B=B, kd=kd, tcc=Inf, t = tp) + rnorm(length(tp), sd = 0.03)
#   }, A, B, kd)
#   ds <- t(ds)
#   ds <- rbind(ds, runif(8, 0, 1))
#
#   rf <- c("p1", "p2", "p2", "p3", "p3", "p3")
#
#   ss <- fitNLSModels(x = ds, f = rf, t = tp, tcc = Inf, type = "deg",
#                      par.init = list(A=0.8, B=0.2, kd=0.04),
#                      par.lower=c(A=0, B=0, kd=0),
#                      par.upper=c(A=1, B=1, kd=10))
#
#   checkEquals(class(ss), "list")
#   checkEquals(names(ss), c("mat", "list", "type"))
#   checkEquals(dim(ss$mat), c(6, 15))
#   checkEquals(rownames(ss$mat), paste("X", 1:6, sep = ""))
#
#   plotCurve.comb(x = ss$list$p2, t = tp, tcc = Inf, curve = "deg")
#   plotCurve.comb(x = ss$list$p3, t = tp, tcc = Inf, curve = "deg")
#
#
#   ds2 <- mapply(function(A, B, kd) {
#     synCurve(A=A, B=B, kd=kd, tcc=Inf, t = tp) + rnorm(length(tp), sd = 0.01)
#   }, A, B, kd)
#   ds2 <- t(ds2)
#   ds2 <- rbind(ds2, runif(8, 0, 1))
#
#   ss2 <- fitNLSModels(x = ds, f = rf, t = tp, tcc = Inf, type = "syn",
#                      par.init = list(A=0.8, B=0.2, ks=0.04),
#                      par.lower=c(A=0, B=0, ks=0),
#                      par.upper=c(A=1, B=1, ks=10))
#
#   ss2$list
# }
#
#
#
#
#
# r <- refitwoOutlier(x = ss$list$p3, include = attr(x, "individual")[, "mse"] < 0.01, t, tcc = Inf)
#
