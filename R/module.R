
# ========================================================================
#
#                          fitting Server
#
# ========================================================================

#' @title shiny module - fitting in server
#' @description Should not be called by users. shiny module for server, curve fitting
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#' @param x a reactive object, matrix, passed to fitNLSModels
#' @param tcc a reactive object, cell doubling time
#' @param f a reactive object, f passed to fitNLSModels
#' @param time a reactive object, time points passed to fitNLSModels
#' @param type a reactive object, the type of fitting (syn or deg) passed to fitNLSModels
#' @param A a reactive object, passed to fitNLSModels
#' @param B a reactive object, passed to fitNLSModels
#' @param par.init a reactive object, passed to fitNLSModels
#' @param par.lower a reactive object, passed to fitNLSModels
#' @param par.upper a reactive object, passed to fitNLSModels
#' @param pre.col a reactive object containing a data.frame, for the annotation columns
#' @param ncore a reactive object, passed to fitNLSModels
#' @param resultPath where the figures to be saved
#' @return a reactive value
#'   reactive( list(pre.col = pre.col(), mat = r()$mat, list = r()$list, type = r()$type) )
#' @import shiny shinyBS 
#' @importFrom grDevices dev.off tiff
#' @importFrom utils write.table
#' @include runShiny.R
#' @keywords internal


fmod <- function(input, output, session, x, tcc = reactive(Inf), f, time, type, A, B,
                 par.init, par.lower, par.upper, pre.col, ncore, resultPath) {

  r <- reactive({
    fitNLSModels(x = x(), f = f(), t = time(), type = type(),
                 A = A(), B = B(), tcc = tcc(),
                 par.init = par.init(), par.lower = par.lower(),
                 par.upper = par.upper(), ncore = ncore())
  })

  k <- reactive( ifelse(type() == "syn", "ks", "kd") )
  ot <- reactive({ cbind(pre.col(), sigDF(r()$mat)) })

  output$tab <- DT::renderDataTable({
    DT::datatable(ot(), filter = "bottom", class = list("nowrap"),
                  selection = "single", options = list(scrollX = TRUE, autoWidth = TRUE),
                  rownames= FALSE)
  }, server = TRUE)

  output$downloadData <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".txt", sep="")
    },
    content = function(file) {
      write.table(ot(), file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    }
  )

  figInd <- reactive({
    req(input$tab_rows_selected)
    i <- input$tab_rows_selected

    err.x <- log(2)/ r()$mat[i, c("ci025", "ci975")]
    err.y <- ifelse(type() == "syn",
                    synCurve(A = r()$mat[i, "A"], B =  r()$mat[i, "B"], r()$mat[i, k()],
                             tcc = tcc(), t = log(2)/ r()$mat[i, k()]),
                    degCurve(A = r()$mat[i, "A"], B =  r()$mat[i, "B"], r()$mat[i, k()],
                             tcc = tcc(), t = log(2)/ r()$mat[i, k()]))

    list(x = x()[i, ], t = time(), tcc = tcc(), A = r()$mat[i, "A"],
         B = r()$mat[i, "B"], k = r()$mat[i, k()], curve = type(),
         err.x = err.x, err.y = err.y)
  })

  figComb <- reactive({
    req(input$tab_rows_selected)
    updateOutlierBox <- TRUE
    if (r()$type == "individual")
      return(NULL)
    i <- input$tab_rows_selected
    ii <- r()$mat[i, 1]

    combList <- r()$list[[ii]]
    if (!is.null(input$outlier))  {
      if (input$outlier[1] != "") {
        iir <- !rownames(attr(combList, "individual")) %in% input$outlier
        if (any(iir)) {
          isolate({
            updateOutlierBox <- FALSE
            combList  <- refitwoOutlier(
              x = combList, include = iir,
              t = time(), A = A(), B = B(), tcc = tcc(), par.init = par.init(),
              par.lower = par.lower(), par.upper = par.upper() )
          })
        }
      }
    }

    list(x = combList, t = time(), tcc = tcc(),
         leg.vec = structure(rownames(x()), names = rownames(x())),
         curve = type(), updateOutlierBox = updateOutlierBox)
  })

  saveCombFigure <- function(ind, comb, file) {
    tiff(file, width = 24, height = 14, res = 150, units = "cm")
    layout(matrix(1:2, 1, 2))
    plotCurve(ind$x, t = ind$t, tcc = ind$tcc, A = ind$A, B = ind$B, k = ind$k,
              add = FALSE, col="red", lineOnly = FALSE, curve = ind$curve, pch=20,
              err.x = ind$err.x, err.y = ind$err.y)

    plotCurve.comb(x = comb$x, t = comb$t, tcc = comb$tcc,
                   leg.vec = comb$leg.vec, curve = comb$curve)
    dev.off()
  }

  output$downloadFigures <- downloadHandler(
    filename = function() {
      paste("figure-", Sys.Date(), ".tif", sep="")
    },
    content = function(file) {
      saveCombFigure(ind = figInd(), comb = figComb(), file = file)
    }
  )

  output$plotInd <- renderPlot({
    req(figInd()$A)
    req(figInd()$B)
    req(figInd()$k)
    plotCurve(figInd()$x, t = figInd()$t, tcc = figInd()$tcc, A = figInd()$A, B = figInd()$B, k = figInd()$k,
              add = FALSE, col="red", lineOnly = FALSE, curve = figInd()$curve, pch=20,
              err.x = figInd()$err.x, err.y = figInd()$err.y)
  })

  observeEvent(input$tab_rows_selected, {
      updateSelectInput(session, "outlier", choices = rownames(attr(figComb()$x, "inputmatrix")), selected = NULL)
  })

  output$plotComb <- renderPlot({
    plotCurve.comb(x = figComb()$x, t = figComb()$t, tcc = figComb()$tcc,
                   leg.vec = figComb()$leg.vec, curve = figComb()$curve)
  })

  observeEvent(input$save2folder, {
    withProgress(message = "making plot", value = 0, {
      unif <- unique(r()$mat[, "f"])
      n <- length(unif)
      dir <- file.path(resultPath, gsub("-| |:|", "", Sys.time()))
      if (!dir.exists(dir))
        dir.create(dir)

      for (isave in unif) {
        incProgress(1/n)
        combList <- r()$list[[isave]]
        comb <- list(x = combList, t = time(), tcc = tcc(),
                     leg.vec = structure(rownames(x()), names = rownames(x())),
                     curve = type())
        file <- paste(dir, "/", isave, ".tiff", sep = '')
        tiff(file, width = 500, height = 400)
        plotCurve.comb(x = comb$x, t = comb$t, tcc = comb$tcc, main = isave,
                       leg.vec = comb$leg.vec, curve = comb$curve, leg.cex = 0.7)
        dev.off()
      }
    })
  })
  reactive( list(pre.col = pre.col(), mat = r()$mat, list = r()$list, type = r()$type) )
}


# ========================================================================
#
#                          Fitting UI
#
# ========================================================================

#' @title shiny module, curve fitting UI
#' @description Should not be called by users. shiny module, curve fitting UI
#' @param id id for namespace
#' @return no value to be returned
#' @keywords internal

fmodUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12, DT::dataTableOutput(ns("tab")))
    ),
    fixedRow(
      column(6,
             wellPanel(
               tags$b("Single fitting"),
               tags$br(),
               tags$br(),
               tags$br(),
               plotOutput(ns("plotInd"))
               )
      ),
      column(6,
             wellPanel(
               fixedRow(
                 column(6, tags$b("Combined fitting")),
                 tags$head(
                   tags$style("label {display:inline;}")
                 ),
                 column(6, selectizeInput(ns("outlier"),"", choices=NULL, multiple=TRUE,
                                          options = list(placeholder = 'Select outlier to exclude ...')))
               ),
               plotOutput(ns("plotComb"))
             )
      )
    ),
    fixedRow(
      column(2, downloadButton(ns("downloadData"), "Download Table")),
      column(3, downloadButton(ns("downloadFigures"), "Download Current Figures")),
      column(3, actionButton(ns("save2folder"), "Save All Comb Figures", icon = icon("download")))
      # column(3, shinyDirButton(ns("save2folder"), "Save all fig to folder", "Select a folder"))

    )
  )
}

# ========================================================================
#
#                          View server
#
# ========================================================================

#' @title shiny module server - view curves
#' @description Should not be called by users. shiny module server - view curves
#' @param input shiny input
#' @param output shiny output
#' @param session shiny session
#' @param deg a reactive object, containing an object returned by fitNLSModels
#' @param syn a reactive object, containing an object returned by fitNLSModels
#' @param deg.ratio a reactive object, a matrix of degradation ratios at differnet time points
#' @param syn.ratio a reactive object, a matrix of synthesis ratios at differnet time points
#' @param deg.t a reactive object, times points for degradation ratios
#' @param syn.t a reactive object, time points for synthesis ratios
#' @param tcc a reactive object, cell doubling time
#' @return write table to disk
#' @keywords internal
#' @importFrom grDevices dev.off tiff
#' @importFrom utils write.table

#'
combView <- function(input, output, session, deg, syn,
                     deg.ratio, syn.ratio, deg.t, syn.t, tcc) {
  tab <- reactive({
    cbind(deg()$pre.col,
          sigDF(deg()$mat),
          sigDF(syn()$mat))
  })

  output$tab <- DT::renderDataTable({
    DT::datatable(tab(), filter = "bottom", class = list("nowrap"),
                  selection = "single", options = list(scrollX = TRUE, autoWidth = TRUE),
                  rownames= FALSE)
  }, server = TRUE)

  output$figInd <- renderPlot({
    req(input$tab_rows_selected)
    i <- input$tab_rows_selected
    req(syn()$mat[i, "B"])
    req(syn()$mat[i, "A"])
    req(deg()$mat[i, "ks"])
    req(deg()$mat[i, "B"])
    req(deg()$mat[i, "A"])
    req(deg()$mat[i, "kd"])
    
    ylim <- c(0, max(c(deg.ratio()[i, ], syn.ratio()[i, ]), na.rm = TRUE))
    plotCurve(deg.ratio()[i, ], t = deg.t(), tcc = tcc(), ylim = ylim,
              A = deg()$mat[i, "A"], B = deg()$mat[i, "B"], k = deg()$mat[i, "kd"],
              add = FALSE, col="red", lineOnly = FALSE, curve = "deg", pch=20)
    plotCurve(syn.ratio()[i, ], t = syn.t(), tcc = tcc(), lty = 2,
              A = syn()$mat[i, "A"], B = syn()$mat[i, "B"], k = syn()$mat[i, "ks"],
              add = TRUE, col="red", lineOnly = FALSE, curve = "syn", pch=1)
  })

  output$figComb <- renderPlot({
    req(input$tab_rows_selected)
    if (deg()$type == "individual")
      return(NULL)
    i <- input$tab_rows_selected
    ii <- deg()$mat[i, 1]
    lv <- structure(rownames(deg.ratio()), names = rownames(deg.ratio()))
    plotCurve.comb(x = deg()$list[[ii]], t = deg.t(), tcc = tcc(),
                   leg.vec = lv, curve = "deg")
    plotCurve.comb(x = syn()$list[[ii]], t = syn.t(), tcc = tcc(),
                   leg.vec = lv, curve = "syn", add = TRUE, pch = 1, lty = 2, legend = FALSE)
  })

  output$downloadFigures <- downloadHandler(
    filename = function() {
      paste("figure-", Sys.Date(), ".tif", sep="")
    },
    content = function(file) {
      tiff(file, width = 24, height = 14, res = 150, units = "cm")
      layout(matrix(1:2, 1, 2))
      i <- input$tab_rows_selected
      ylim <- c(0, max(c(deg.ratio()[i, ], syn.ratio()[i, ]), na.rm = TRUE))
      plotCurve(deg.ratio()[i, ], t = deg.t(), tcc = tcc(), ylim = ylim,
                A = deg()$mat[i, "A"], B = deg()$mat[i, "B"], k = deg()$mat[i, "kd"],
                add = FALSE, col="red", lineOnly = FALSE, curve = "deg", pch=20)
      plotCurve(syn.ratio()[i, ], t = syn.t(), tcc = tcc(), lty = 2,
                A = syn()$mat[i, "A"], B = syn()$mat[i, "B"], k = syn()$mat[i, "ks"],
                add = TRUE, col="red", lineOnly = FALSE, curve = "syn", pch=1)
      if (deg()$type != "individual") {
        ii <- deg()$mat[i, 1]
        lv <- structure(rownames(deg.ratio()), names = rownames(deg.ratio()))
        plotCurve.comb(x = deg()$list[[ii]], t = deg.t(), tcc = tcc(),
                       leg.vec = lv, curve = "deg")
        plotCurve.comb(x = syn()$list[[ii]], t = syn.t(), tcc = tcc(),
                       leg.vec = lv, curve = "syn", add = TRUE, pch = 1, lty = 2, legend = FALSE)
      }
      dev.off()
    }
  )

  output$downloadData <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".txt", sep="")
    },
    content = function(file) {
      write.table(tab(), file, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    }
  )

}

# ========================================================================
#
#                          View UI
#
# ========================================================================

#' @title shiny module, view UI
#' @description Should not be called by users. shiny module, view UI
#' @param id id for namespace
#' @return no value to be returned
#' @keywords internal

combViewUI <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(12, DT::dataTableOutput(ns("tab")))
    ),
    fixedRow(
      column(6,
             wellPanel(
               tags$b("Single fitting"),
               plotOutput(ns("figInd"))
             )
      ),
      column(6,
             wellPanel(
               tags$b("Combined fitting"),
               # column(6, selectizeInput(ns("outlier"), "outlier", choices=NULL, multiple=TRUE,
               #                          options = list(placeholder = 'Select outlier to exclude ...'))),
               plotOutput(ns("figComb"))
             )
      )
    ),
    fixedRow(
      column(2, downloadButton(ns("downloadData"), "Download Table")),
      column(2, downloadButton(ns("downloadFigures"), "Download Figures"))
    )
  )
}
