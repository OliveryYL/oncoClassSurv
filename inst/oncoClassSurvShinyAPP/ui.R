library(shiny)
library(shinyWidgets)
library(shinyjs)#(重置更新)
library(shinycssloaders)#github开发版#https://github.com/daattali/shinycssloaders

# Define the js method that resets the page
#定义函数
jsResetCode <- "shinyjs.resetPage = function() {history.go(0)}"

shinyUI(fluidPage(
  theme = shinythemes::shinytheme("cerulean"),#superhero,united,cerulean,...
  #shinythemes::themeSelector(),
  #使用themeSelector自定义主题风格
  tags$div(id = "selector", shinythemes::themeSelector()),
  #自定义的CSS样式，用来设置下拉框的初始状态为隐藏，并且设置一个类名，例如"show"，用来表示下拉框的显示状态
  tags$style(HTML("
  #selector {
    display: none;
  }
  .show {
    display: block;
  }
")),
#监听按钮的点击事件，并且切换下拉框的显示或隐藏状态
tags$script(HTML("
  $(document).ready(function() {
    $('#toggle').click(function() {
      $('#selector').toggleClass('show');
    });
  });
")),


  useShinyjs(), # Include shinyjs in the UI (重置更新)
  extendShinyjs(text = jsResetCode, functions = "resetPage"), # Add the js code to the page

  # 创建一个自定义的CSS样式
  tags$style("
    .btn {
      color: black;
    }
    .btn:hover {
      color: red;
    }
  "),

  titlePanel("oncoClassSurv prediction"),
  sidebarLayout(
    sidebarPanel(
                 tabsetPanel(
                   tabPanel("Data",

                            fluidRow(fileInput(inputId="file.train.exp", label = h4("Expression profile for training"))),
                            fluidRow(fileInput(inputId="file.input.exp", label = h4("Expression profile for input"))),

                            fluidRow(fileInput(inputId="file.train.clin", label = h4("Clinical data for training"))),
                            fluidRow(fileInput(inputId="file.input.clin", label = h4("Clinical data for input"))),
                            fluidRow(fileInput(inputId="file.train.marker", label = h4("Marker genes"))),#train_cluster.feature.path
                            fluidRow(fileInput(inputId="file.train.Prog.features", label = h4("Prognostic features"))),#train_survival.feature.path


                            # Create an action button with the label "Run"
                            actionButton("run", "Run",icon("paper-plane",lib = "font-awesome"),
                                         style="color: red; background-color: #337ab7;
                                         border-color: #2e6da4")),

                   tabPanel("Parameter",
                            fluidRow(column(4, offset = 0,
                                            selectInput("s1.task", label = h5("Task"),
                                                        choices = list("1-Classification" = 1, "2-Prognosis" = 2,
                                                                       "3-Both" = 3), selected = 1),
                                            selectInput("s2.rm.batch.effect", label = h5("BatchEffect"),
                                                        choices = list("Remove" = 1, "Ignore" = 2), selected = 1),
                                            selectInput("s3.plot.combatch", label = h5("PlotBatchPCA"),
                                                        choices = list("TRUE" = 1, "FALSE" = 2), selected = 1)#,
                                            # selectInput("sx.miss_go.on", label = h5("miss_go.on"),
                                            #             choices = list("TRUE" = 1, "FALSE" = 2), selected = 1)
                            ),
                            column(4, offset = 1,
                                   selectInput("s4.cluster.method", label = h5("ClassMethod"),
                                               choices = list("RF" = 1, "SVM" = 2), selected = 2),
                                   selectInput("s5.kernel", label = h5("Kernel"),
                                               choices = list("radial" = 1, "linear" = 2,
                                                              "polynomial" = 3,"sigmoid" = 4),
                                               selected = 1),
                                   selectInput("s6.plot.surv.curve", label = h5("PlotSurvCurve"),
                                               choices = list("TRUE" = 1, "FALSE" = 2), selected = 1)#,
                                   # selectInput("sy.importance", label = h5("importance"),
                                   #             choices = list("TRUE" = 1, "FALSE" = 2), selected = 1)
                            )),

                            fluidRow(column(3, offset = 0,radioButtons("exp.type", label = h5("ExprFormat"),
                                                                       choices = list("FPKM" = 1, "TPM" = 2),
                                                                       selected = 1)),
                                     column(3, offset = 1,radioButtons("limmaNorm.train", label = h5("TrainNorm"),
                                                                       choices = list("Yes" = 1, "No" = 2),
                                                                       selected = 2)),
                                     column(3, offset = 1,radioButtons("limmaNorm.input", label = h5("InputNorm"),
                                                                       choices = list("YES" = 1, "No" = 2),
                                                                       selected = 2))),
                            fluidRow(column(3, offset = 0,radioButtons("RandomSamTrain", label = h5("RandomSamTrain"),
                                                                       choices = list("Yes" = 1, "No" = 2),
                                                                       selected = 2)),
                                     column(3, offset = 1,numericInput("TrainSeed", label = h5("TrainSeed"), value = 1234,
                                                                       min = 1, max = 9999, step = 1)),
                                     column(3, offset = 1,numericInput("SamplingProb", label = h5("SamplingProb"), value = 0.7,
                                                                       min = 0.01, max = 1, step = 0.01))),
                            fluidRow(fileInput(inputId="file.PartSamples.Pred", label = h5("Predict parts of samples"))),#PartSamples.Pred.Path

                            fluidRow(column(12, offset = 0,sliderInput("ntree", label = h5("ntree"),
                                                                       min = 0, max = 5000, value = c(3000)))),
                            fluidRow(column(12, offset = 0,sliderInput("nodesize", label = h5("nodesize"),
                                                                       min = 0, max = 100, value = c(3)))),
                            fluidRow(column(12, offset = 0,sliderInput("mtry", label = h5("mtry"),
                                                                       min = 0, max = 100, value = c(5)))),
                            fluidRow(column(12, offset = 0,sliderInput("cost", label = h5("cost"),
                                                                       min = 0, max = 100, value = c(4)))),

                            fluidRow(column(12, offset = 0,pickerInput("plot.samples", label = h6("SurvSample"),
                                                                      choices = 1:10, multiple = TRUE,
                                                                      options = list(`actions-box` = TRUE)))),

                            fluidRow(column(3, offset = 0,numericInput("survcurve.break.x.by", label = h6("SurvBreak"),
                                                                       value = 12,  min = 0, max = 365, step = 1)),
                                     column(3, offset = 0,textInput(inputId="survtime.unit", label = h5("SurvUnit"),
                                                                    value = "months",placeholder ="months...")),
                                     column(3, offset = 0,numericInput("leg.size", label = h6("LegendSize"),
                                                                       value = 8,  min = 1, max = 20, step = 1))),


                            #自定义感兴趣的基因（自定义勾选(clustergenes/immunegene)/输入基因名称/上传基因集文件），绘制热图
                            fluidRow(
                              column(width=5,offect=0,fileInput(inputId="selgene.uplod",
                                                                label = h5("Upload Geneset"))),
                              column(width=5,offect=1,
                                     selectInput("sel.geneset", label = h5("Default Geneset"),
                                                 choices = list("Marker" = 1,
                                                                "Immune" = 2),
                                                 selected = 2))
                              ),

                            fluidRow(textInput(inputId="text.input.gene.explore", label = h5("Explore interesting genes"),
                                               placeholder ="TP53...")),#input text,e.g., gene names
                            fluidRow(column(width=12, verbatimTextOutput("value.gene"))),
                            fluidRow(column(width=5.5,offect=0,selectInput("select.hp.gene", label = h5("Heat Map Geneset"),
                                                                  choices = list("Default" = 1, "Upload" = 2, "TextInput"=3),
                                                                  selected = 1))),
                            checkboxGroupInput("hp.sample", label = h5("Heat Map Options"),
                                               choices = list("Selected Samples" = 1, "Show Sample Name" = 2),
                                               selected = NULL),
                            actionButton("reset_button", "Reset",icon("redo",lib = "font-awesome"),
                                                  style="color: red; background-color: #337ab7;
                                         border-color: #2e6da4")
                            ),

                   tabPanel("Usage",
                            p(span("The flowchart and examples are as follows.",style = "color:blue")),
                            img(src="flowchart_example.png", height = "40%", width = "100%"),
                            p(span("The usage of parameters is as follows.",style = "color:blue")),
                            img(src="parameter_usage.png", height = "100%", width = "100%"),
                            p(span(em("Claim: For academic research only. Contact us for cooperation."),style = "color:red; font-size:10px")),
                            br(),
                            p("The ",span("oncoClassSurv", style = "color:blue"),"application is based on the R package called",
                              span("oncoClassSurv", style = "color:blue")," on GitHub",
                              ' so you can visit the ',
                              a("GitHub repository",
                                href = "https://github.com/OliveryYL/oncoClassSurv")," for more information.")
                    ),
                   tabPanel("Theme",
                            p(span(em("Click the button below to switch your favorite theme style."),style = "color:black; font-size:16px")),
                            br(),
                            br(),
                            actionButton("toggle", "Show/Hide Theme Selector"))
                   )
                 ),
    mainPanel(
      tabsetPanel(
        tabPanel("Classification",
                 textOutput("class.method"), # 添加一个文本输出元素
                 br(),
                 fluidRow(
                   column(6, offset = 0, plotOutput("Plot_class.summary")),
                   column(6, offset = 0, plotOutput("Plot_hp"))),
                 br(),
                 fluidRow(
                   column(6, offset = 0, downloadButton(outputId="download.class.count",
                                                        label ="Download classifications",
                                                        style = "margin-bottom: 10px; align: center;")),
                   column(6, offset = 0, downloadButton(outputId="download.class.hp",
                                                        label ="Download heat map",
                                                        style = "margin-bottom: 10px; align: center;")))),

        tabPanel("Prognosis",
                 fluidRow(
                   column(12, plotOutput("ggsurv.curve"))),
                 fluidRow(
                   column(12, offset = 0, downloadButton(outputId="download.survplot",
                                                        label ="Download survival curves",
                                                        style = "margin-bottom: 10px; align: center;")))),

        tabPanel("Table",
                 fluidRow(
                   column(6, offset = 0, tableOutput("class.table")),
                   column(6, offset = 0, tableOutput("surv.table"))),
                 br(),
                 fluidRow(
                   column(6, offset = 0, downloadButton(outputId="download.class",
                                               label ="Download classification results",
                                               style = "margin-bottom: 10px; align: center;")),
                   column(6, offset = 0, downloadButton(outputId="download.surv",
                                                      label ="Download prognostic results",
                                                      style = "margin-bottom: 10px; align: center;")))),

        tabPanel("More",
                 h4(em("Check batch effect between the train and input:")),
                 br(),
                 fluidRow(column(12, offset = 0,plotOutput("Combat"))),
                 br(),
                 fluidRow(column(12, offset = 0, downloadButton(outputId="download.combat",
                                label ="Download Combat PCA",
                                style = "margin-bottom: 10px; align: center;"))))

        )
      )
    )
  )
)
