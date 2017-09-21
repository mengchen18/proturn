server <- function(input, output, session) {
  figureFolder = get("path", .__shiny.var__)
  # ============== Update inputs ===================
  dataInput <- reactive({
    req(input$file)
    data <- read.delim(input$file$datapath, as.is = TRUE, header = TRUE)
    header <- colnames(data)
    headerNum <- header[sapply(data, is.numeric)]
    list(data=data, header=header, headerNum=headerNum)
  })

  # update
  observe( updateSelectizeInput(session, "id.col", choices = dataInput()$header, selected = dataInput()$header[1],
                                options = list(placeholder = 'Select a column ...')))
  observe( updateSelectizeInput(session, "export.col", choices = dataInput()$header,
                                selected = c(input$id.col, input$collapse.col, isolate(input$export.col)),
                                options = list(placeholder = 'Select columns ...')))
  observe( updateSelectizeInput(session, "collapse.col", choices = dataInput()$header, selected = "",
                                options = list(placeholder = 'Select a column ...')))
  observe( updateSelectizeInput(session, "deg.col", choices = dataInput()$headerNum, selected = NULL,
                                options = list(placeholder = 'Select columns ...')))
  observe( updateSelectizeInput(session, "B.col", choices = dataInput()$headerNum, selected = "",
                                options = list(placeholder = 'Select a column ...')))
  observe( updateSelectizeInput(session, "syn.col", choices = dataInput()$headerNum, selected = NULL,
                                options = list(placeholder = 'Select columns ...')))
  observe( updateSelectizeInput(session, "B.col.syn", choices = dataInput()$headerNum, selected = "",
                                options = list(placeholder = 'Select a column ...')))

  observe({
    if (length(timePointsDeg()) == length(input$deg.col))
      updateCollapse(session, "fit.deg.param",
                     style = list("Basic settings" = "success",
                                  "Advanced settings" = "success")) else
                                    updateCollapse(session, "fit.deg.param",
                                                   style = list("Basic settings" = "default",
                                                                "Advanced settings" = "default"))
  })
  observe({
    if (length(timePointsSyn()) == length(input$syn.col))
      updateCollapse(session, "fit.syn.param",
                     style = list("Basic settings" = "success",
                                  "Advanced settings" = "success")) else
                                    updateCollapse(session, "fit.syn.param",
                                                   style = list("Basic settings" = "default",
                                                                "Advanced settings" = "default"))
  })


  # ==================== intermediate calculation =====================

  timePointsSyn <- reactive({
    req(input$timePointsSyn)
    as.numeric(strsplit(input$timePointsSyn, ',')[[1]])
  })
  timePointsDeg <- reactive({
    req(input$timePointsDeg)
    as.numeric(strsplit(input$timePointsDeg, ',')[[1]])
  })
  degRatio <- reactive({
    req(input$deg.col)
    x <- dataInput()$data[, input$deg.col, drop=FALSE]
    rownames(x) <- dataInput()$data[, input$id.col]
    x
  })
  synRatio <- reactive({
    req(input$syn.col)
    x <- dataInput()$data[, input$syn.col, drop=FALSE]
    rownames(x) <- dataInput()$data[, input$id.col]
    x
  })

  fac <- reactive({
    if (is.null(input$collapse.col) || input$collapse.col == "") {
      if (is.null(input$id.col) || input$id.col == "")
        r <- rownames(dataInput()$data) else
          r <- dataInput()$data[, input$id.col]
    } else {
      r <- dataInput()$data[, input$collapse.col]
    }
    r
  })

  pre.col <- reactive({
    if (is.null(input$export.col) || input$export.col == "") {
      if (!is.null(rownames(dataInput()$data)))
        r <- data.frame(row.names = rownames(dataInput()$data), stringsAsFactors = FALSE)
    }
    r <- dataInput()$data[, input$export.col, drop = FALSE]
    r
  })

  ### ===================== output ======================

  res <- reactiveValues(deg = NULL, syn = NULL)
  observeEvent(input$go, {
    B <- NULL
    if (input$B.col != "")
      B <- dataInput()$data[, input$B.col]

    init <- list(A = input$deg.Ainit,
                 B = input$deg.Binit,
                 kd = input$deg.kdinit)
    low <- c(A = input$deg.Arange[1],
             B = input$deg.Brange[1],
             kd = input$deg.kdrange[1])
    up <- c(A = input$deg.Arange[2],
            B = input$deg.Brange[2],
            kd = input$deg.kdrange[2])
    
    res$deg <- callModule(proturn:::fmod, "DEG", x = reactive(as.matrix(degRatio(), rownames.force = TRUE)),
                          f = fac, time = timePointsDeg, type = reactive("deg"),
                          tcc = reactive(as.numeric(input$tcc)),
                          A = reactive(NULL), B = reactive(B),
                          par.init = reactive(init),
                          par.lower = reactive(low),
                          par.upper = reactive(up),
                          pre.col = pre.col,
                          ncore = reactive(input$ncore), 
                          resultPath = figureFolder)
  })

  observeEvent(input$go, {
    B <- NULL
    if (input$B.col.syn != "")
      B <- dataInput()$data[, input$B.col.syn]

    init <- list(A = input$syn.Ainit,
                 B = input$syn.Binit,
                 ks = input$syn.kdinit)
    low <- c(A = input$syn.Arange[1],
             B = input$syn.Brange[1],
             ks = input$syn.kdrange[1])
    up <- c(A = input$syn.Arange[2],
            B = input$syn.Brange[2],
            ks = input$syn.kdrange[2])
    res$syn <- callModule(proturn:::fmod, "SYN", x = reactive(as.matrix(synRatio(), rownames.force = TRUE)),
                          f = fac, time = timePointsSyn, type = reactive("syn"),
                          tcc = reactive(as.numeric(input$tcc)),
                          A = reactive(NULL), B = reactive(B),
                          par.init = reactive(init),
                          par.lower = reactive(low),
                          par.upper = reactive(up),
                          pre.col = pre.col,
                          ncore = reactive(input$ncore),
                          resultPath = figureFolder)
  })

  observeEvent(input$go, {
    req(res$deg)
    req(res$syn)
    callModule(proturn:::combView, "CMB", deg = res$deg, syn = res$syn,
               deg.ratio = reactive(as.matrix(degRatio(), rownames.force = TRUE)),
               syn.ratio = reactive(as.matrix(synRatio(), rownames.force = TRUE)),
               syn.t = timePointsSyn, deg.t = timePointsDeg,
               tcc = reactive(as.numeric(input$tcc)))
  })

  output$dt.degRatio <- DT::renderDataTable({
    DT::datatable(sigDF(degRatio()), filter = "bottom", selection = "single", 
                  options = list(scrollX = TRUE), caption = "Degradation ratios")
  })
  
  output$dt.synRatio <- DT::renderDataTable({
    DT::datatable(sigDF(synRatio()), filter = "bottom", selection = "single", 
                  options = list(scrollX = TRUE), caption = "Synthesis ratios")
  })
}


