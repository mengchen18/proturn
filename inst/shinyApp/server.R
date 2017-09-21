server <- function(input, output, session) {
  figureFolder <- get("path", .__proturn.shiny.var__)
  
  # ============== Update inputs ===================
  
  fpath <- reactiveValues(path = NULL)
  observeEvent(input$exData, {
    fpath$path <- system.file("example.data.txt", package = "proturn")
  })
  observeEvent(input$file, {
    req(input$file)
    fpath$path <- input$file$datapath
  })
  
  dataInput <- reactive({
    req(fpath$path)
    data <- read.delim(fpath$path, as.is = TRUE, header = TRUE)
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
  output$dt.degRatio <- DT::renderDataTable({
    DT::datatable(sigDF(degRatio()), filter = "bottom", selection = "single", 
                  options = list(scrollX = TRUE), caption = "Degradation ratios")
  })
  
  output$dt.synRatio <- DT::renderDataTable({
    DT::datatable(sigDF(synRatio()), filter = "bottom", selection = "single", 
                  options = list(scrollX = TRUE), caption = "Synthesis ratios")
  })
  
  
  ###
  modinputf <- eventReactive(input$go, fac())
  modinputTcc <- eventReactive(input$go, as.numeric(input$tcc))
  modinputNcore <- eventReactive(input$go, input$ncore)
  
  ###################### degradation #####################
  modinputDeg <- eventReactive(input$go, as.matrix(degRatio(), rownames.force = TRUE))
  modinputDegTime <- eventReactive(input$go, timePointsDeg())
  modinputDegB <- eventReactive(input$go, {
    B <- NULL
    if (input$B.col != "")
      B <- dataInput()$data[, input$B.col]
    B
  })
  modinputDegInit <- eventReactive(input$go, {
    list(A = input$deg.Ainit, B = input$deg.Binit, kd = input$deg.kdinit)
  })
  modinputDegLow <- eventReactive(input$go, {
    c(A = input$deg.Arange[1], B = input$deg.Brange[1], kd = input$deg.kdrange[1])
  })
  modinputDegUp <- eventReactive(input$go, {
    c(A = input$deg.Arange[2], B = input$deg.Brange[2], kd = input$deg.kdrange[2])
  })
  deg <- callModule(proturn:::fmod, "DEG", 
                        x = modinputDeg, time = modinputDegTime,
                        f = modinputf, type = reactive("deg"),
                        tcc = modinputTcc,
                        A = reactive(NULL), B = modinputDegB,
                        par.init = modinputDegInit,
                        par.lower = modinputDegLow,
                        par.upper = modinputDegUp,
                        pre.col = pre.col,
                        ncore = modinputNcore, 
                        resultPath = figureFolder)
  
  ###################### synthesis #####################
  modinputSyn <- eventReactive(input$go, as.matrix(synRatio(), rownames.force = TRUE))
  modinputSynTime <- eventReactive(input$go, timePointsSyn())
  modinputSynB <- eventReactive(input$go, {
    B <- NULL
    if (input$B.col.syn != "")
      B <- dataInput()$data[, input$B.col.syn]
    B
  })
  modinputSynInit <- eventReactive(input$go, {
    list(A = input$syn.Ainit, B = input$syn.Binit, ks = input$syn.kdinit)
  })
  modinputSynLow <- eventReactive(input$go, {
    c(A = input$syn.Arange[1], B = input$syn.Brange[1], ks = input$syn.kdrange[1])
  })
  modinputSynUp <- eventReactive(input$go, {
    c(A = input$syn.Arange[2], B = input$syn.Brange[2], ks = input$syn.kdrange[2])
  })
  
  syn <- callModule(proturn:::fmod, "SYN", f = modinputf, 
                    x = modinputSyn, time = modinputSynTime, 
                    type = reactive("syn"),
                    tcc = modinputTcc,
                    A = reactive(NULL), B = modinputSynB,
                    par.init = modinputSynInit,
                    par.lower = modinputSynLow,
                    par.upper = modinputSynUp,
                    pre.col = pre.col,
                    ncore = modinputNcore,
                    resultPath = figureFolder)
  
  ###################### degradation #####################
  observeEvent(input$go, {
    req(deg())
    req(syn())
    callModule(proturn:::combView, "CMB", deg = deg, syn = syn,
               deg.ratio = reactive(as.matrix(degRatio(), rownames.force = TRUE)),
               syn.ratio = reactive(as.matrix(synRatio(), rownames.force = TRUE)),
               syn.t = timePointsSyn, deg.t = timePointsDeg,
               tcc = reactive(as.numeric(input$tcc)))
  })
  
  session$onSessionEnded(function() {
    rm(.__proturn.shiny.var__)
  })
}


