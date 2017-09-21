#' @title Lanch shiny application for protein degradation/synthesis curve fitting
#' @description Lanch shiny application for protein degradation/synthesis curve fitting
#' @param maxFileSize max size of files can be uploaded
#' @param figureFolder a character string to specify a directory where the figures to be saved
#' @return A shiny application will be launched, not values to be returned
#' @export
#' @importFrom parallel mclapply
#' @examples 
#' # not run
#' 1
#' # runShiny()

runShiny <- function(maxFileSize = 100*1024^2, figureFolder = getwd()) {
  # find and launch the app
  
  s <- Sys.info()[["sysname"]]
  absPath <- FALSE
  if (s == "Windows") {
    substr(figureFolder, 1, 2) %in% paste0(LETTERS, ":")
    absPath <- TRUE
  } else {
    substr(figureFolder, 1, 1) %in% c("/", "~")
    absPath <- TRUE
  }
  
  if (!absPath) {
    stop("figureFolder need an ABSOLUTE path!")
  }
  
  options(shiny.maxRequestSize=maxFileSize)
  assign(".__proturn.shiny.var__", new.env(), envir = globalenv())
  .__proturn.shiny.var__$path <- figureFolder
  appDir <- system.file("shinyApp", package = "proturn")
  shiny::runApp(appDir, display.mode = "normal")
}


