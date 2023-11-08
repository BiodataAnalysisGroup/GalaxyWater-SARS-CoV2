# load libraries
library(shiny)
library(DT)
library(data.table)
library(plotly)
library(shinyjs)
# setwd("CERTH/Visualizations_for_Galaxy_workflows/cojac-freyja-output-examples/shiny_app_galaxy_sars_cov_2")

# load Freyja preprocessing function
source("freyja_output_preprocessing.R")

# Load galaxy workflow results --------------------------------------------
galaxy_res = fread("us-freyja-tsv-mix-dates.tsv") # Demo file??? See server side for Galaxy results upload per user???

# Number of lineages
lineages_number = paste(galaxy_res$lineages, collapse = " ") |>
  strsplit(split = " ") |>
  unlist() |>
  unique() |>
  length()

# input-widget starting values
startDate = min(galaxy_res$Date)
endDate = max(galaxy_res$Date)

# Define UI
ui <- navbarPage(
  
  # Resetable with shinyjs
  useShinyjs(),
  
  # navbarPage options ------------------------------------------------------
  title = "Galaxy SARS-CoV-2",
  id = "mainNavBarPage",
  # header = "Test header",
  # footer = "Test footer",
  # theme = ,
  # icon = ,
  includeCSS("www/custom.css"),

  # barplot tab -------------------------------------------------------------
  tabPanel(
    title = "Barplot",
    id = "barplotTab",
    sidebarLayout(
      sidebarPanel = sidebarPanel(
        # id = "testSideBarID", 
        # shinyjs::useShinyjs(),
        width = 2,
        h3('SARS-CoV-2 barplot'),
        hr(),
        # div(
        #   id = "barplotWidgets",
        #   
        # ),
        # Plot by Date or Samples (x-axis)
        radioButtons(
          inputId = "barplotXaxisSelection", 
          label = h4("Χ-axis values"),
          choices = list(
            "Date" = 1, 
            "Sample" = 2),
          selected = 1),
        hr(),
        # Date range
        dateRangeInput(
          inputId = 'barplotDateRange',
          label = h4("Date range"),
          start = startDate,
          end = endDate, 
          min = startDate, 
          max = endDate
        ),
        hr(), 
        # Cluster variable
        selectInput(
          inputId = "barplotSelectCluster",
          label = h4("Cluster variable"),
          choices = list(  
            "None",
            "Sample",
            "Variant",
            "Percentage",
            "Date",
            "test_clusters"
          ), # fix with columns from outputs as choices (not from the initial table)
          selected = "None"
        ),
        hr(),
        # Percentage
        sliderInput(
          inputId = "barplotPercentage",
          label = h4("Percentage (%)"),
          min = 0,
          max = 100,
          value = 0
        ),
        hr(),
        # Reset
        actionButton(
          inputId = "barplotResetButton", 
          label = "Reset" ,
          icon = icon("wand-magic-sparkles")
        ),
        hr(),
        # Session information
        h4("Session information"),
        textOutput("pipeline"),
        textOutput("fileName"),
        textOutput("fileDate"),
        textOutput("NoSamples"),
        textOutput("NoLineages"),
        hr(),
        # Upload user data (compatible with Galaxy output) (???)
        fileInput("file", label = h4("File input"))
      ),
      mainPanel = mainPanel(
        width = 10,
        # Reactive barplot --------------------------------------------------------
        plotlyOutput(
          outputId = 'barplotGraph',
          height = "700px"
        ),
        # DataTable ---------------------------------------------------------------
        tabsetPanel(
          id = "barplotDTpanel",
          tabPanel(
            title = "Summarized",
            DT::dataTableOutput(
              outputId = 'barplotDTsummarized'
            )
          ),
          tabPanel(
            title = "Lineages",
            DT::dataTableOutput(
              outputId = 'barplotDTlineages'
            )
          )
        )
        # -------------------------------------------------------------------------
      )
    )
  ),
  
  # lineplot tab ------------------------------------------------------------
  tabPanel(
    title = "Lineplot",
    id = "lineplotTab",
    sidebarLayout(
      sidebarPanel = sidebarPanel(width = 2,
                                  h2('test')),
      mainPanel = mainPanel(width = 10,
                            h2('test'))
    )
  )
  #
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  # Reset graph -------------------------------------------------------------
  observeEvent(
    input$barplotResetButton,
    {
      reset("barplotDateRange")
      reset("barplotSelectCluster")
      reset("barplotPercentage")
    }
  )

  # Process Galaxy SARS-CoV-2 pipeline results ------------------------------

  # generate list object (server-side, once per user)
  # [[1]] summarized results
  # [[2]] lineage results
  freyjaList = freyja_results_preprocessing(x = "us-freyja-tsv-mix-dates.tsv")

  # -------------------------------------------------------------------------

  # barplot Tab - graph
  output$barplotGraph <- renderPlotly({
    # generate plot
    res_plot = freyja_results_graph(
      x = freyjaList,
      from = input$barplotDateRange[1], # starting from
      to = input$barplotDateRange[2], # ending on
      xAxisSelection = input$barplotXaxisSelection, # select to plot by Date (default) or by Sample
      percentageThreshold = input$barplotPercentage, # % MIN threshold
      samplesCluster = input$barplotSelectCluster # cluster variable
    )
    #  by default:
    # hide display mode bar
    # compare data on hover
    res_plot <- res_plot  |> 
      layout(
        hovermode = 'x'
      )
    # config(displayModeBar = FALSE) |>
    
    # return plot
    return(res_plot)
  })
  
  # -------------------------------------------------------------------------

  # barplot Tab - DT - summarized
  output$barplotDTsummarized <- DT::renderDataTable({
    # generate DT
    res_dt = freyjaList[[1]]
    # User input selection
    # date range
    res_dt = res_dt[Date >= input$barplotDateRange[1] & Date <= input$barplotDateRange[2]]
    # MIN percentage threshold for report
    res_dt = res_dt[Percentage >= input$barplotPercentage/100]
    # return DT
    return(res_dt)
  })
  
  # barplot Tab - DT - lineages
  output$barplotDTlineages <- DT::renderDataTable({
    # generate DT
    res_dt = freyjaList[[2]]
    # User input selection
    # date range
    res_dt = res_dt[Date >= input$barplotDateRange[1] & Date <= input$barplotDateRange[2]]
    # MIN percentage threshold for report
    res_dt = res_dt[abundance >= input$barplotPercentage/100]
    # return DT
    return(res_dt)
  })
  

  # Shiny callback functions to update DT tables ----------------------------
  
  # Shiny callback function to update summarized DT table based on 
  # selected values of the barplot graph
  
  # observeEvent(event_data("plotly_selected"), {
  #   # extract selected data from barplot
  #   selected_data <- event_data(event = "plotly_selected")
  #   # check for NULL
  #   if (!is.null(selected_data)) {
  #     # point index
  #     selected_indices <- selected_data$pointNumber
  #     # selected_data <- data[selected_indices, ]
  #     output$barplotDTsummarized <- DT::renderDataTable({
  #       # output$barplotDTsummarized[selected_indices, ] # CAN NOT ACCESS IT THIS WAY!
  #     })
  #   }
  # })
  
  # Shiny callback function to update lineages DT table based on 
  # selected observations from the summarized DT table
  
  # # Reactive expression to filter the data based on the selected variable
  # filtered_data <- reactive({
  #   selected_var <- input$input_table_rows_selected
  #   if (length(selected_var) == 0) {
  #     return(NULL)
  #   } else {
  #     selected_var <- input_data[selected_var, "Variable"]
  #     # Filter the output data based on the selected variable
  #     output_data <- input_data[input_data$Variable %in% selected_var, ]
  #     return(output_data)
  #   }
  # })
  # 
  # # Create the output DT table based on the selected variable
  # output$output_table <- renderDT({
  #   datatable(filtered_data(), options = list(pageLength = 5))
  # })
  
  # Shiny callback function to update summarized DT table based on 
  # selected observations from the lineages DT table
  
  # output$distPlot <- renderPlot({
  #     # generate bins based on input$bins from ui.R
  #     x    <- faithful[, 2]
  #     bins <- seq(min(x), max(x), length.out = input$bins + 1)
  # 
  #     # draw the histogram with the specified number of bins
  #     hist(x, breaks = bins, col = 'darkgray', border = 'white',
  #          xlab = 'Waiting time to next eruption (in mins)',
  #          main = 'Histogram of waiting times')
  # })
  
  # Session information
  output$pipeline <- renderText({
    paste0("Galaxy pipeline: ", "Freyja")
  })
  output$fileName <- renderText({
    paste0("File: ", "Test filename")
  })
  output$fileDate <- renderText({
    paste0("Date: ", "Test date")
  })
  output$NoSamples <- renderText({
    paste0("Samples: ", as.character(nrow(galaxy_res)))
  })
  output$NoLineages <- renderText({
    paste0("Lineages: ", as.character(lineages_number))
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
