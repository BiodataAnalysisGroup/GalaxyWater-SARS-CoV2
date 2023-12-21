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
source("lineagespot_functions.R")

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
      # Create graph
      actionButton(
        inputId = "barplotCreateButton",
        label = "Create graph" ,
        icon = icon("chart-column")
      ),
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
      # Software pipeline used
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
        # Update -> Subset
        actionButton(
          inputId = "barplotUpdateButton",
          label = "Subset graph" ,
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
      plotlyOutput(outputId = 'dataGraph',
                   height = "700px"),
      
      # tabsetPanel(
      #   id = "graphPanel",
      #   tabPanel(
      #     title = "Freyja",
      #     plotlyOutput(outputId = 'freyjaGraph',
      #                  height = "700px")
      #   ),
      #   tabPanel(
      #     title = "lineagespot",
      #     plotlyOutput(outputId = 'lineagespotGraph',
      #                  height = "700px")
      #   )
      # ),
      
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
  # Freyja:
  # [[1]] summarized results
  # [[2]] lineage results
  # lineagespot:
  # [[1]] lineage results

  # reactive values
  dataInitial = reactiveVal(NA)
  dataSummarized = reactiveVal(NA)
  dataLineages = reactiveVal(NA)

  # -------------------------------------------------------------------------

  # initial data
  observeEvent(
    input$barplotCreateButton, # c(input$upload_data, input$upload_metadata),
    {
      # require both data & metadata
      req(input$upload_data)
      req(input$upload_metadata)

      # access file name and path
      inputData = fread(input$upload_data$datapath)
      inputMetadata = fread(input$upload_metadata$datapath)
      
      if(input$barplotToolSelection == 1){ # Freyja
        
        # preprocess pipeline results
        result = freyja_results_preprocessing(
          data = inputData,
          metadata = inputMetadata
        )
        
        # return merged data tables
        dataInitial(result[[1]])
        dataSummarized(result[[1]])
        dataLineages(result[[2]])
        
      } else if(input$barplotToolSelection == 2){ # lineagespot
        
        # preprocess pipeline results
        result = lineagespot_results_preprocessing(
          data = inputData,
          metadata = inputMetadata
        )
        
        # return merged data tables
        dataInitial(result)
        dataSummarized(data.frame())
        dataLineages(result)
         
      }
      
      
    }
  )
  
  # -------------------------------------------------------------------------

  # barplot Tab - graph
  output$dataGraph <- renderPlotly({
    
    # require create graph button
    req(input$barplotCreateButton)
    
    # req(
    #   dataInitial,
    #   dataSummarized,
    #   dataLineages
    # )
    # 
    # # require create graph button
    # req(input$barplotCreateButton)
    # # require input data & metadata upload
    # req(input$upload_data)
    # req(input$upload_metadata)
    
    if(input$barplotToolSelection == 1){ # Freyja
      
      # generate barplot
      res_barplot = freyja_barplot(
        x = dataSummarized(), 
        from = input$barplotDateRange[1], # starting from
        to = input$barplotDateRange[2], # ending on
        xAxisSelection = input$barplotXaxisSelection, # select to plot by Date or by Sample (default)
        percentageThreshold = input$barplotPercentage, # % MIN threshold
        samplesCluster = input$barplotSelectCluster # cluster variable
      )
      
    } else if(input$barplotToolSelection == 2){ # lineagespot
      
      # generate barplot
      res_barplot = lineagespot_barplot(
        x = dataLineages(), 
        from = input$barplotDateRange[1], # starting from
        to = input$barplotDateRange[2], # ending on
        xAxisSelection = input$barplotXaxisSelection, # select to plot by Date or by Sample (default)
        percentageThreshold = input$barplotPercentage, # % MIN threshold
        samplesCluster = input$barplotSelectCluster # cluster variable
      )
      
    }
    
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
    
    # require create graph button
    req(input$barplotCreateButton)
    
    # req(
    #   dataInitial,
    #   dataSummarized,
    #   dataLineages
    # )
    # 
    # # require create graph button
    # req(input$barplotCreateButton)
    # # require input data & metadata upload
    # req(input$upload_data)
    # req(input$upload_metadata)
    
    if(input$barplotToolSelection == 1){
      # generate DT
      res_dt = dataSummarized()
      # User input selection
      # date range
      res_dt = res_dt[Date >= input$barplotDateRange[1] & Date <= input$barplotDateRange[2]]
      # MIN percentage threshold for report
      res_dt = res_dt[Percentage >= input$barplotPercentage/100]
    } else if(input$barplotToolSelection == 2){
      # generate DT
      res_dt = dataSummarized()
    }
    
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
        write.csv(dataInitial()[searched_rows,], file)  
      } else {
        # write the entire processed table
        write.csv(dataInitial(), file)  
      }
    }
  )
  
  # -------------------------------------------------------------------------

  # Tab - DT - lineages
  output$barplotDTlineages <- DT::renderDataTable({
    
    # require create graph button
    req(input$barplotCreateButton)
    
    # req(
    #   dataInitial,
    #   dataSummarized,
    #   dataLineages
    # )
    # 
    # # require create graph button
    # req(input$barplotCreateButton)
    # # require input data & metadata upload
    # req(input$upload_data)
    # req(input$upload_metadata)
    
    # generate DT
    res_dt = dataLineages()
    # User input selection
    # date range
    res_dt = res_dt[Date >= input$barplotDateRange[1] & Date <= input$barplotDateRange[2]]
    
    
    if(input$barplotToolSelection == 1){ # Freyja
      # MIN percentage threshold for report
      res_dt = res_dt[abundance >= input$barplotPercentage/100]
    } else if(input$barplotToolSelection == 2){ # lineagespot
      # MIN percentage threshold for report
      res_dt = res_dt[Percentage >= input$barplotPercentage/100] 
    }
    
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
        write.csv(dataLineages()[searched_rows,], file)  
      } else {
        # write the entire processed table
        write.csv(dataLineages(), file)  
      }
    }
  )
  
  # -------------------------------------------------------------------------

  # observe Update button and change dataset for graph
  observeEvent(
    input$barplotUpdateButton,
    {
      ###
      if(input$barplotToolSelection == 1){ # Freyja
        
        # show searched data from barplot
        req(input$barplotDTsummarized_rows_all)
        searched_rows <- input$barplotDTsummarized_rows_all
        
        # check for selection event
        if(!is.null(searched_rows)){
          dataSummarized(
            dataInitial()[searched_rows,]
          )
        }
      } else if(input$barplotToolSelection == 2){ # lineagespot
        # show searched data from barplot
        req(input$barplotDTlineages_rows_all)
        searched_rows <- input$barplotDTlineages_rows_all
        
        # check for selection event
        if(!is.null(searched_rows)){
          dataLineages(
            dataInitial()[searched_rows,]
          )
        }
      }
      ###
    }
  )

  # -------------------------------------------------------------------------

  # Session information
  output$pipeline <- renderText({
    
    # 
    req(input$barplotCreateButton)
    
    # # require both data & metadata
    # req(input$upload_data)
    # req(input$upload_metadata)
    
    # check tool selected
    if(input$barplotToolSelection == 1){ # Freyja
      tool = "Freyja"
    } else if(input$barplotToolSelection == 2){ # lineagespot
      tool = "lineagespot"
    }
    
    # return software tool used
    return(
      paste0("Galaxy pipeline: ", tool)
    )
    
  })
  output$fileName <- renderText({
    
    # 
    req(input$barplotCreateButton)
    
    
    # # require both data & metadata
    # req(input$upload_data)
    # req(input$upload_metadata)
    
    paste0("File: ", input$upload_data$name)
  })
  output$fileDate <- renderText({
    
    # 
    req(input$barplotCreateButton)
    
    
    # # require both data & metadata
    # req(input$upload_data)
    # req(input$upload_metadata)
    
    paste0("Date: ", file.info(input$upload_data$datapath)$mtime)
  })
  output$NoSamples <- renderText({
    
    # 
    req(input$barplotCreateButton)
    
    # # require both data & metadata
    # req(input$upload_data)
    # req(input$upload_metadata)
    
    # check tool selected
    if(input$barplotToolSelection == 1){ # Freyja (uppercase)
      noSamples = paste0("Samples: ", dataInitial()$Sample |> unique() |> length() |> as.character() )
    } else if(input$barplotToolSelection == 2){ # lineagespot (lowercase)
      noSamples = paste0("Samples: ", dataInitial()$sample |> unique() |> length() |> as.character() )
    }
    
    # return number of samples
    return(noSamples)
    
  })
  output$NoLineages <- renderText({
    
    # 
    req(input$barplotCreateButton)
    
    if(input$barplotToolSelection == 1){ # Freyja
      paste0("Lineages: ", dataLineages()$lineage |> unique() |> length() |> as.character())
    } else if(input$barplotToolSelection == 2){ # lineagespot
      paste0("Lineages: ", dataInitial()$lineage |> unique() |> length() |> as.character())
    }
    
    # # require both data & metadata
    # req(input$upload_data)
    # req(input$upload_metadata)
    
    
  })
  
  # -------------------------------------------------------------------------
  
  observeEvent(
    
    input$barplotCreateButton,
    # c(input$upload_data, input$upload_metadata),
    {
      # require both data & metadata
      req(input$upload_data)
      req(input$upload_metadata)
      
      # calculate values
      startVal = dataInitial()$Date |> min()
      endVal = dataInitial()$Date |> max()
      
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
      startVal = dataInitial()$Date |> min()
      endVal = dataInitial()$Date |> max()
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
      
      if(input$barplotToolSelection == 1){ # Freyja
        dataSummarized( dataInitial() )
      } else if(input$barplotToolSelection == 2){ # lineagespot
        dataLineages( dataInitial() )
      }

    }
  )
  
}

# Run the application 
shinyApp(ui = ui, server = server)
