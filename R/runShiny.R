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

runShiny <- function(maxFileSize = 100*1024^2, figureFolder = "~") {
  # find and launch the app
  options(shiny.maxRequestSize=maxFileSize)
  .__shiny.var__ <- environment()
  .__shiny.var__$path <- figureFolder
  appDir <- system.file("shinyApp", package = "proturn")
  shiny::runApp(appDir, display.mode = "normal")
}



