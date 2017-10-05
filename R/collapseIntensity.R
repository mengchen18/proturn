#' @title Collapse rows of an intensity matrix according to a vector/factor
#' @description Collapse rows of an intensity matrix according to a vector/factor. 
#' @param x a \code{data.frame} containing at least the intensity columns and 
#'   a column according to which intensities are collapsed.
#' @param intensity.col the columns of intensities, could be either column names
#'   or column indices
#' @param collapse.col the column accroding to which intensity are collapsed. Most of 
#'   the cases, this column should be peptide sequences or protein IDs, see \code{details}. 
#' @param annot.col the annotation columns, which will be retained in the output. 
#' @param n.cores the number of cores to be used, it is passed to \code{mclapply}
#' @param collapse.method the collapse method, see \code{details} section.
#' @details This function prepares input for \code{\link{fitDegNLS}}. Two typical 
#'   usages are: 1) collapse intensitys on peptide level, i.e. if multiple rows 
#'   corresponds to identical peptide, combine them into one rows in some way (see below).
#'   In this case, the peptide sequence should be used as \code{collapse.col}. 
#'   2) collapse from peptide level to protein levels, so the protein name or 
#'   ID could be used as \code{collapse.col}. 
#'   
#'   The following collapse methods are avaliable:
#'   \itemize{
#'     \item sum sum the intensities;
#'     \item mean mean of intensities;
#'     \item median median of intensities;
#'     \item max.range retain the one with max range;
#'     \item max.range01 retain the one with max difference between first and last intensity columns.
#'   }
#' @author Chen Meng
#' @examples
#'   df <- data.frame(
#'     sequence = c("CCAD", "CSDS", "IIURE", "IIURE", "XMNT"),
#'     t0 = rnorm(5),
#'     t1 = rnorm(5),
#'     t2 = rnorm(5),
#'     name = LETTERS[1:5]
#'   )
#'   collapseIntesity(x = df, collapse.col = "sequence", intensity.col = c("t0", "t1", "t2"), annot.col = "name")
#'   collapseIntesity(x = df, collapse.col = 1, intensity.col = 2:4, annot.col = 5, collapse.method = "max.range01")
#' @return a \code{data.frame} of collapsed intensities, containing the collapsed.col, intensity.col and annot.cols.
#' @export
#' @importFrom parallel mclapply
#' @importFrom matrixStats colMedians rowMaxs rowMins
#'

collapseIntesity <- function(x, intensity.col, collapse.col = "Sequence", annot.col, n.cores = 1,
                    collapse.method = c("sum", "mean", "median", "max.range", "max.range01")[1]) {

  collapse.method <- match.arg(collapse.method[1], c("sum", "mean", "median", "max.range", "max.range01"))

  collapseFun.core <- switch (collapse.method,
    "sum" = colSums,
    "mean" = colMeans,
    "median" = colMedians,
    "max.range" = function(x, na.rm = FALSE) {
      i <- which.max(rowMaxs(x, na.rm = na.rm) - rowMins(x, na.rm = na.rm))
      x[i, ]
    },
    "max.range01" = function(x, na.rm = FALSE) {
      i <- which.max(abs(x[, 1] - x[, ncol(x)]))
      x[i, ]
    }
  )
  collapseFun <- function(x, na.rm = FALSE) {
    if (is.vector(x))
      r <- x else if (is.matrix(x)) {
        if (nrow(x) == 1)
          r <- x[1, ] else
            r <- collapseFun.core(x, na.rm = na.rm)
      } else
        stop ("Un-recognized data type (should be vector or matrix)")
    r
  }

  if (inherits(intensity.col, c("numeric", "integer")))
    intensity.col <- colnames(x)[intensity.col]
  if (inherits(collapse.col, c("numeric", "integer")))
    collapse.col <- colnames(x)[collapse.col]
  if (inherits(annot.col, c("numeric", "integer")))
    annot.col <- colnames(x)[annot.col]

  intensity <- apply(x[, intensity.col, drop = FALSE], 2, function(x) as.numeric(as.character(x)))
  retain <- x[, annot.col, drop = FALSE]
  collapse <- as.character(x[, collapse.col])

  vecs <- mclapply(unique(collapse), function(key) {
    i <- key == collapse
    c(key=key,
      sapply(retain[i, , drop = FALSE], function(x) paste(unique(x), collapse = "||")),
      collapseFun(intensity[i, ], na.rm = TRUE)
      )
  }, mc.cores = n.cores)

  mat <- do.call(rbind, vecs)
  df <- data.frame(mat, stringsAsFactors = FALSE)

  # colnames(df)[-(1:(length(annot.col)+1))] <- intensity.col
  colnames(df) <- c(collapse.col, annot.col, intensity.col)
  df[intensity.col] <- lapply(df[, intensity.col, drop = FALSE], as.numeric)
  df
}
