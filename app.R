# load libraries --------------------
library(shiny)
library(DT)
library(data.table)
library(plotly)
library(shinyjs)
library(stringr)
library(dplyr)
library(paletteer) # for colors see: https://r-charts.com/color-palettes/

# load Freyja pipeline preprocessing and graph functions
source("Freyja_functions.R")

# Define UI
ui <- fluidPage(
  
  # Resetable with shinyjs
  useShinyjs(),
  
  # -------------------------------------------------------------------------
  
  id = "mainFluidPage",
  includeCSS("www/custom.css"),
  
  # -------------------------------------------------------------------------
  
  sidebarLayout(
    sidebarPanel = sidebarPanel(
      width = 2,
      h3('SARS-CoV-2 dashboard'),
      hr(),
      # Upload user data (compatible with Galaxy output)
      fileInput(inputId = "upload_data", label = "Upload data"),
      hr(),
      # Upload user metadata (compatible with Galaxy output)
      fileInput("upload_metadata", label = "Upload metadata"),
      hr(),
      downloadButton("downloadSummarized", "Download summarized", class = "customButton"),
      downloadButton("downloadLineages", "Download lineages"),
      hr(),
      # Plot by Date or Samples (x-axis)
      radioButtons(
        inputId = "barplotToolSelection",
        label = "Software tool",
        # h4("Χ-axis values"),
        choices = list("Freyja" = 1,
                       "lineagespot" = 2),
        selected = 1
      ),
      hr(),
      # Plot by Date or Samples (x-axis)
      radioButtons(
        inputId = "barplotXaxisSelection",
        label = "Χ-axis values",
        # h4("Χ-axis values"),
        choices = list("Date" = 1,
                       "Sample" = 2),
        selected = 2
      ),
      hr(),
      # Date range
      dateRangeInput(
        inputId = 'barplotDateRange',
        label = "Date range",
        # h4("Date range"),
        start = NULL,
        end = NULL,
        min = NULL,
        max = NULL
      ),
      hr(),
      # Cluster variable
      selectInput(
        inputId = "barplotSelectCluster",
        label = "Cluster variable",
        # h4("Cluster variable"),
        choices = list(
          "None",
          "Sample",
          "Variant",
          "Percentage",
          "Date",
          "pseudo_clusters"
        ),
        # fix with columns from outputs as choices (not from the initial table)
        selected = "None"
      ),
      hr(),
      # Percentage
      sliderInput(
        inputId = "barplotPercentage",
        label = "Percentage (%)",
        min = 0,
        max = 100,
        value = 0
      ),
      hr(),
      tags$span(
        # Reset
        actionButton(
          inputId = "barplotResetButton",
          label = "Reset" ,
          icon = icon("broom")
        ),
        # Update
        actionButton(
          inputId = "barplotUpdateButton",
          label = "Update graph" ,
          icon = icon("wand-magic-sparkles")
        )
      ),
      hr(),
      # Session information
      tags$strong("Session information"),
      tags$div(
        class = "sessionInfo",
        textOutput("pipeline"),
        textOutput("fileName"),
        textOutput("fileDate"),
        textOutput("NoSamples"),
        textOutput("NoLineages")
      ),
      hr()
    ),
    mainPanel = mainPanel(
      width = 10,
      # Reactive graphs
      tabsetPanel(
        id = "graphPanel",
        tabPanel(
          title = "Freyja",
          plotlyOutput(outputId = 'freyjaGraph',
                       height = "700px")
        ),
        tabPanel(
          title = "lineagespot",
          plotlyOutput(outputId = 'lineagespotGraph',
                       height = "700px")
        )
      ),
      # DataTables
      tabsetPanel(
        id = "DTpanel",
        tabPanel(title = "Summarized",
                 DT::dataTableOutput(outputId = 'barplotDTsummarized')),
        tabPanel(title = "Lineages",
                 DT::dataTableOutput(outputId = 'barplotDTlineages'))
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {

  # Process Galaxy SARS-CoV-2 pipeline results ------------------------------

  # generate list object (server-side, once per user)
  # [[1]] summarized results
  # [[2]] lineage results

  # reactive values
  freyjaSummarized_initial = reactiveVal(NA)
  freyjaSummarized = reactiveVal(NA)
  freyjaLineages = reactiveVal(NA)
  
  # initial data
  observeEvent(
    c(input$upload_data, input$upload_metadata),
    {
      # require both data & metadata
      req(input$upload_data)
      req(input$upload_metadata)

      # access file name and path
      inputData = fread(input$upload_data$datapath)
      inputMetadata = fread(input$upload_metadata$datapath)

      # preprocess Freyja pipeline results
      result = freyja_results_preprocessing(
        data = inputData,
        metadata = inputMetadata
      )

      # return merged data tables
      freyjaSummarized_initial(result[[1]])
      freyjaSummarized(result[[1]])
      freyjaLineages(result[[2]])
    }
  )
  
  # -------------------------------------------------------------------------

  # barplot Tab - graph
  output$freyjaGraph <- renderPlotly({
    
    # require input data & metadata upload
    req(input$upload_data)
    req(input$upload_metadata)
    
    # generate barplot
    res_barplot = freyja_barplot(
      x = freyjaSummarized(), 
      from = input$barplotDateRange[1], # starting from
      to = input$barplotDateRange[2], # ending on
      xAxisSelection = input$barplotXaxisSelection, # select to plot by Date or by Sample (default)
      percentageThreshold = input$barplotPercentage, # % MIN threshold
      samplesCluster = input$barplotSelectCluster # cluster variable
    )
    
    # # by default:
    # # hide display mode bar
    # # compare data on hover
    # res_barplot <- res_barplot  |>
    #   layout(
    #     hovermode = 'x'
    #   )
    # # config(displayModeBar = FALSE) |>
    
    # return plot
    return(res_barplot)
  })

  # -------------------------------------------------------------------------

  # Tab - DT - summarized
  output$barplotDTsummarized <- DT::renderDataTable({
    # require input data & metadata upload
    req(input$upload_data)
    req(input$upload_metadata)
    # generate DT
    res_dt = freyjaSummarized()
    # User input selection
    # date range
    res_dt = res_dt[Date >= input$barplotDateRange[1] & Date <= input$barplotDateRange[2]]
    # MIN percentage threshold for report
    res_dt = res_dt[Percentage >= input$barplotPercentage/100]
    # return DT
    return(res_dt)
  })
  
  # -------------------------------------------------------------------------

  # download summarized DT as CSV
  output$downloadSummarized <- downloadHandler(
    filename = "summarized.csv",
    content = function(file) {
      # show searched data from barplot
      # req(input$barplotDTsummarized_rows_all)
      searched_rows <- input$barplotDTsummarized_rows_all
      # check for selection event
      if(!is.null(searched_rows)){
        # subset for search results
        write.csv(freyjaSummarized_initial()[searched_rows,], file)  
      } else {
        # write the entire processed table
        write.csv(freyjaSummarized_initial(), file)  
      }
    }
  )
  
  # -------------------------------------------------------------------------

  # Tab - DT - lineages
  output$barplotDTlineages <- DT::renderDataTable({
    # require input data & metadata upload
    req(input$upload_data)
    req(input$upload_metadata)
    # generate DT
    res_dt = freyjaLineages()
    # User input selection
    # date range
    res_dt = res_dt[Date >= input$barplotDateRange[1] & Date <= input$barplotDateRange[2]]
    # MIN percentage threshold for report
    res_dt = res_dt[abundance >= input$barplotPercentage/100]
    # return DT
    return(res_dt)
  })
  
  # -------------------------------------------------------------------------
  
  # download summarized DT as CSV
  output$downloadLineages <- downloadHandler(
    filename = "lineages.csv",
    content = function(file) {
      # show searched data from barplot
      # req(input$barplotDTlineages_rows_all)
      searched_rows <- input$barplotDTlineages_rows_all
      # check for selection event
      if(!is.null(searched_rows)){
        # subset for search results
        write.csv(freyjaLineages()[searched_rows,], file)  
      } else {
        # write the entire processed table
        write.csv(freyjaLineages(), file)  
      }
    }
  )
  
  # -------------------------------------------------------------------------

  # observe Update button and change dataset for graph
  observeEvent(
    input$barplotUpdateButton,
    {
      # show searched data from barplot
      req(input$barplotDTsummarized_rows_all)
      searched_rows <- input$barplotDTsummarized_rows_all

      # check for selection event
      if(!is.null(searched_rows)){
        freyjaSummarized(
          freyjaSummarized_initial()[searched_rows,]
        )
      }
    }
  )

  # -------------------------------------------------------------------------

  # Session information
  output$pipeline <- renderText({
    # require both data & metadata
    req(input$upload_data)
    req(input$upload_metadata)
    paste0("Galaxy pipeline: ", "Freyja")
  })
  output$fileName <- renderText({
    # require both data & metadata
    req(input$upload_data)
    req(input$upload_metadata)
    paste0("File: ", input$upload_data$name)
  })
  output$fileDate <- renderText({
    # require both data & metadata
    req(input$upload_data)
    req(input$upload_metadata)
    paste0("Date: ", file.info(input$upload_data$datapath)$mtime)
  })
  output$NoSamples <- renderText({
    # require both data & metadata
    req(input$upload_data)
    req(input$upload_metadata)
    paste0("Samples: ", freyjaSummarized_initial()$Sample |> unique() |> length() |> as.character() )
  })
  output$NoLineages <- renderText({
    # require both data & metadata
    req(input$upload_data)
    req(input$upload_metadata)
    paste0("Lineages: ", freyjaLineages()$lineage |> unique() |> length() |> as.character())
  })
  
  # reactive values
  freyjaSummarized_initial = reactiveVal(NA)
  freyjaSummarized = reactiveVal(NA)
  freyjaLineages = reactiveVal(NA)
  
  # -------------------------------------------------------------------------
  
  observeEvent(
    c(input$upload_data, input$upload_metadata),
    {
      # require both data & metadata
      req(input$upload_data)
      req(input$upload_metadata)
      
      # calculate values
      startVal = freyjaSummarized()$Date |> min()
      endVal = freyjaSummarized()$Date |> max()
      
      # update dateRangeInput widget
      updateDateRangeInput(
        inputId = 'barplotDateRange',
        start = startVal,
        min = startVal, 
        end = endVal, 
        max = endVal
      )
    }
  )
  
  # Reset graph -------------------------------------------------------------
  observeEvent(
    input$barplotResetButton,
    {
      # calculate values
      startVal = freyjaSummarized_initial()$Date |> min()
      endVal = freyjaSummarized_initial()$Date |> max()
      # update dateRangeInput widget
      updateDateRangeInput(
        inputId = 'barplotDateRange',
        start = startVal,
        min = startVal, 
        end = endVal, 
        max = endVal
      )
      reset("barplotSelectCluster")
      reset("barplotPercentage")
      freyjaSummarized( freyjaSummarized_initial() )
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)
