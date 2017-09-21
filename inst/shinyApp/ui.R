ui <- shinyUI(
  pageWithSidebar(
    titlePanel("PROTURN: Protein Turnover Curve Fitting Tool"),
    sidebarPanel(
      tabsetPanel(
        tabPanel("General settings",
                 fileInput("file", "Upload file",
                           multiple=FALSE,
                           accept = c('text/csv',
                                      'text/comma-separated-values',
                                      'text/tab-separated-values',
                                      'text/plain',
                                      '.csv',
                                      '.tsv')),
                 selectizeInput("id.col", "ID column",
                                choices=NULL, multiple=FALSE,
                                options = list(placeholder = 'Waiting data ...')),
                 selectizeInput("collapse.col", "Collapse on column",
                                choices=NULL, multiple=TRUE,
                                options=list(placeholder = "waiting data ...")),
                 selectizeInput("export.col", "Export columns",
                                choices=NULL, multiple=TRUE,
                                options=list(placeholder = "waiting data ...")),
                 sliderInput("ncore", "Number of cores to use", 1, min(20, detectCores()), value = 4)
        ),
        tabPanel("Model parameters",
                 p("This application using nonlinear least square method
               to fit protein turnover curves. For the degradation curve, the following
               model is fitted:"),
                 p(withMathJax("$$f(t)=(A-B)\\cdot e^{-(k_d+\\frac{ln2}{t_{cc}})\\cdot t}+B$$")),
                 p("The model fitted for the synthesis curve is as:"),
                 p(withMathJax("$$f(t)=(B-A)\\cdot e^{-(k_s+\\frac{ln2}{t_{cc}})\\cdot t}+A$$")),
                 p("where t is the time points (given in hours);
               kd is the degradation consistant to be estimated.
               In addtion, two other parameters are included:
               A is the (normalized) amplitude, B is accounting
               for the offset seen in data, which could be attributed
               by the recycling of amino acids."),
                 hr(),
                 textInput("tcc", "Doubling Time (Hour)", value="Inf"),
                 tags$b("Fitting degradation curves"),
                 bsCollapse(id = "fit.deg.param", open = "Basic settings", multiple = TRUE,
                            bsCollapsePanel("Basic settings",
                                            selectizeInput("deg.col", "degradation columns",
                                                           choices=NULL, multiple=TRUE,
                                                           options = list(placeholder = 'Waiting data ...')),
                                            fixedRow(
                                              column(7, textInput("timePointsDeg", "Time points (Hour)", value="0,1,3,6,10,16,24,34,48")),
                                              column(5, selectizeInput("B.col", "B column",
                                                                       choices=NULL, multiple=FALSE,
                                                                       options = list(placeholder = 'Waiting data ...')))
                                            ),
                                            style = "default"),
                            bsCollapsePanel("Advanced settings",
                                            fixedRow(
                                              column(4, h5("Initial A"), numericInput("deg.Ainit", NULL, value=0.9)),
                                              conditionalPanel("input[['B.col']] == ''",
                                                               column(4, h5("Initial B"), numericInput("deg.Binit", NULL, value=0.1))
                                              ),
                                              column(4, h5("Initial Kd"), numericInput("deg.kdinit", NULL, value=0.04))
                                            ),
                                            fixedRow(
                                              column(4, h5("Range A"),
                                                     sliderInput("deg.Arange", NULL, min = 0, max = 2, step=0.01, value=c(0, 1))),
                                              conditionalPanel("input[['B.col']] == ''",
                                                               column(4, h5("Range B"),
                                                                      sliderInput("deg.Brange", NULL, min = -1, max = 1, step = 0.1, value=c(0, 1)))
                                              ),
                                              column(4, h5("Range Kd"),
                                                     sliderInput("deg.kdrange", NULL, min = 0, max = 10, step=0.05, value=c(0, 10)))
                                            ),
                                            style = "default")
                 ),
                 tags$b("Fitting synthesis curves"),
                 bsCollapse(id = "fit.syn.param", open = NULL, multiple = TRUE,
                            bsCollapsePanel("Basic settings",
                                            selectizeInput("syn.col", "Synthesis columns",
                                                           choices=NULL, multiple=TRUE,
                                                           options = list(placeholder = 'Waiting data ...')),
                                            fixedRow(
                                              column(7, textInput("timePointsSyn", "Time points (Hour)", value="1,3,6,10,16,24,34,48,96")),
                                              column(5, selectizeInput("B.col.syn", "B column",
                                                                       choices=NULL, multiple=FALSE,
                                                                       options = list(placeholder = 'Waiting data ...'))
                                              )
                                            ),
                                            style = "default"),
                            bsCollapsePanel("Advanced settings",
                                            fixedRow(
                                              conditionalPanel("input[['B.col.syn']] == ''",
                                                               column(4, h5("Initial A"), numericInput("syn.Ainit", NULL, value=0.9))
                                              ),
                                              column(4, h5("Initial B"), numericInput("syn.Binit", NULL, value=0.1)),
                                              column(4, h5("Initial Kd"), numericInput("syn.kdinit", NULL, value=0.04))
                                            ),
                                            fixedRow(
                                              conditionalPanel("input[['B.col.syn']] == ''",
                                                               column(4, h5("Range A"),
                                                                      sliderInput("syn.Arange", NULL, min = 0, max = 2, step=0.01, value=c(0, 1)))
                                              ),
                                              column(4, h5("Range B"),
                                                     sliderInput("syn.Brange", NULL, min = -1, max = 1, step = 0.1, value=c(0, 1))),
                                              column(4, h5("Range Kd"),
                                                     sliderInput("syn.kdrange", NULL, min = 0, max = 10, step=0.05, value=c(0, 10)))
                                            ),
                                            style = "default")
                 ),
                 actionButton("go", "Run!", icon = icon("calculator"))
        )
      )
    ),
    mainPanel(
      navbarPage("Results",
                 # tabPanel("Introduction",
                 #          tags$video(src = "/home/chen/CloudChen/Projects/ProteinTurnOverRate/MSPD/vid/shinyKD.webm",
                 #                     type = "video/webm", autoplay = NA, height = "100%", width = "600")
                 #          ),
                 tabPanel("Ratio Tables",
                          DT::dataTableOutput("dt.degRatio"),
                          DT::dataTableOutput("dt.synRatio")),
                 tabPanel("Degradation Curve Fitting",
                          proturn:::fmodUI("DEG")),
                 tabPanel("Synthesis Curve Fitting",
                          proturn:::fmodUI("SYN")),
                 tabPanel("Synthesis & Degradation",
                          proturn:::combViewUI("CMB"))
      )
    )
  )
)

