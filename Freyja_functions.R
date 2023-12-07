# Freyja pipeline results -------------------------------------------------

# Freyja pipeline output preprocessing ------------------------------------

freyja_results_preprocessing = function(data, metadata){
  
  # rename input data 1st column (Sample IDs)
  colnames(data)[1] = 'Sample'
  
  # merge input Data and Metadata based on Sample IDs (Sample (first) column)
  freyja_res = merge(
    x = data, 
    y = metadata,
    by.x = 'Sample',
    by.y = 'Sample',
    all.x = TRUE
  )
  
  # rename date column to 'Date'
  if(
    grepl( fixed = T, pattern = 'date', x = tolower(  colnames(freyja_res))) |> any()
  ){
    date_idx = grep( fixed = T, pattern = 'date', x = tolower(  colnames(freyja_res)))
    colnames(freyja_res)[date_idx] = 'Date'
  }
  
  # remove ".fastq" suffix
  freyja_res$Sample = str_split(freyja_res$Sample, "\\.", simplify = TRUE)[, 1]
  
  # extract freyja summarized results ---------------------------------------
  freyja_summarized = freyja_res$summarized |>
    str_remove_all("\\[|\\]") |>
    str_split("\\)\\,\\ \\(") |>
    lapply(str_remove_all, "\\(|\\)") |>
    lapply(function(x) {
      data.table(
        "Variant" = str_split(x, "\\,", simplify = TRUE)[, 1] |> str_remove_all("\\'"),
        "Percentage" = str_split(x, "\\,", simplify = TRUE)[, 2] |> as.numeric()
      )
    })
  # assign sample IDs
  names(freyja_summarized) = freyja_res$Sample
  # convert to data.table
  freyja_summarized = freyja_summarized |> rbindlist(idcol = "Sample")
  
  # extract Freyja lineage results ------------------------------------------
  freyja_lineages = mapply(
    FUN = function(x, y) {
      # create data.tables
      data.table("lineage"   = x,
                 "abundance" = y |> as.numeric())
    },
    str_split(freyja_res$lineages, "\\ "),
    str_split(freyja_res$abundances, "\\ "),
    SIMPLIFY = FALSE
  )
  # assign sample IDs
  names(freyja_lineages) = freyja_res$Sample
  # convert to data.table
  freyja_lineages = freyja_lineages |> rbindlist(idcol = "Sample")
  
  # Index for additional (unknown) columns
  idx = !colnames(freyja_res) %in% c("Sample",
                                     "summarized",
                                     "lineages",
                                     "abundances",
                                     "resid",
                                     "coverage")
  unknown_colnames = colnames(freyja_res)[idx]
  
  # match with Sample column (IDs) and add to summarized and lineages DTs
  for(j in seq(length(unknown_colnames)) ){
    # column name
    colName = unknown_colnames[j] # colnames(freyja_res[,idx, with = FALSE])[i]
    # column values
    colValues = as.data.frame(freyja_res)[, colName] # freyja_res[,colName, with = FALSE]
    # add column to summarized and lineages DTs
    freyja_summarized = freyja_summarized[
      , (colName) :=  colValues[ freyja_summarized$Sample |> rleid() ]
    ]
    freyja_lineages = freyja_lineages[
      , (colName) :=  colValues[ freyja_lineages$Sample |> rleid() ]
    ]
  }
  
  # return pre-processed data
  freyja_preprocessed = list(
    freyja_summarized,
    freyja_lineages
  )
  names(freyja_preprocessed) = c(
    "summarized", 
    "lineages"
  )
  return(
    freyja_preprocessed
  )
}

# Freyja graph generation -------------------------------------------------

# Barplot generation ----------------------------------------------------

freyja_barplot = function(
    x,    
    from, 
    to, 
    percentageThreshold,
    xAxisSelection,
    samplesCluster = NULL
  ){ 
  
  # User input selection
  # date range
  x = x[Date >= from & Date <= to]
  
  # MIN percentage threshold for report
  x = x[Percentage >= percentageThreshold/100]
  
  # barplot: Date vs Perc ---------------------------------------------------

  # tooltip label 
  # https://stackoverflow.com/questions/34605919/formatting-mouse-over-labels-in-plotly-when-using-ggplotly
  
  # x$tooltip_label = paste0(
  #   x$Variant, ": ", round(100 * x$Percentage, digits = 3), "%"
  # )
  # 
  # # summarize and generate reactive label
  # x = x[, by = Sample, tooltip_label := paste(tooltip_label, collapse = "\n")]
  
  x$Percentage = x$Percentage * 100
  
  # Graph drawing objective
  if( as.integer(xAxisSelection) == 1 & (samplesCluster == "None" | is.null(samplesCluster) )){
    # 1. Time series
    
    # update tooltip label
    # x$tooltip_label = paste(
    #   x$tooltip_label,
    #   paste0("Sample: ", x$Sample),
    #   sep = '\n'
    # )
    
    # static ggplot2
    # gr = ggplot(data = x, aes(text = tooltip_label)) +
    #   
    #   geom_bar(
    #     aes(x = Date, y = Percentage, fill = Variant),
    #     position = "fill", # position_dodge2(preserve = "single"), # "dodge", #
    #     stat = "identity"
    #   ) +
    #   
    #   scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
    #   # scale_x_date(
    #   #   expand = c(0, 0), 
    #   #   date_breaks = paste0("1 ", scaleDate)
    #   #   ) +
    #   
    #   scale_fill_manual(values = paletteer_d("ggsci::default_igv")) +
    #   
    #   theme_minimal() +
    #   
    #   theme(
    #     panel.grid = element_blank(),
    #     
    #     axis.ticks = element_line(linewidth = .3),
    #     axis.text.x = element_text(
    #       size = 8,
    #       angle = 45,
    #       hjust = 1
    #     ),
    #     
    #     plot.margin = margin(20, 20, 20, 20)
    #   ) + 
    #   
    #   labs( x = "", y = "" )
    # 
    # # interactive plotly
    # gr_plotly = ggplotly(gr, tooltip = c("text", "x"), dynamicTicks = TRUE) 
      
      
      
      gr_plotly = plot_ly(
          x, x = ~Date, y = ~Percentage, type = 'bar', 
          color = ~Variant, stroke = ~Variant, strokes = "grey50"
      ) |> 
          layout(
              yaxis = list(ticksuffix = "%"),
              barmode = 'stack', hovermode = "x unified"
          )
    
  } else if( as.integer(xAxisSelection) == 2 & (samplesCluster == "None" | is.null(samplesCluster) )) {
    
    # 2.1. Samples/No clusters

    # update tooltip label
    # x$tooltip_label = paste(
    #   x$tooltip_label,
    #   paste0("Date: ", x$Date),
    #   sep = '\n'
    # )
    
    # static ggplot2
    # gr = ggplot(data = x, aes(text = tooltip_label)) +
    #   
    #   geom_bar(
    #     aes(x = reorder(Sample, Date), y = Percentage, fill = Variant),
    #     position = "fill",
    #     stat = "identity"
    #   ) +
    #   
    #   scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
    #   # scale_x_date(expand = c(0, 0), date_breaks = paste0("1 ", scaleDate)) +
    #   
    #   scale_fill_manual(values = paletteer_d("ggsci::default_igv")) +
    #   
    #   theme_minimal() +
    #   
    #   theme(
    #     panel.grid = element_blank(),
    #     axis.ticks.y = element_line(linewidth = .3),
    #     axis.ticks.x = element_blank(),
    #     axis.text.x = element_blank(),
    #     plot.margin = margin(20, 20, 20, 20)
    #   ) + 
    #   
    #   labs( x = "", y = "" ) # +
    #   # facet_grid(as.formula(paste0(".~ ", samplesCluster)), scales = "free_x")
    # 
    # # interactive plotly
    # gr_plotly = ggplotly(
    #   gr, tooltip = c("text", "x"), dynamicTicks = F
    # )
      
      gr_plotly = plot_ly(
          x, x = ~Sample, y = ~Percentage, type = 'bar', 
          color = ~Variant, stroke = ~Variant, strokes = "grey50"
      ) |> 
          layout(
              yaxis = list(ticksuffix = "%"),
              barmode = 'stack', hovermode = "x unified"
          )
    
  } else {
    # 2.2. Samples/Clusters
    
    # update tooltip label
    # x$tooltip_label = paste(
    #   x$tooltip_label,
    #   paste0("Date: ", x$Date),
    #   sep = '\n'
    # )
    
    # static ggplot2
    # gr = ggplot(data = x, aes(text = tooltip_label)) +
    #   
    #   geom_bar(
    #     aes(x = reorder(Sample, Date), y = Percentage, fill = Variant),
    #     position = "fill",
    #     stat = "identity"
    #   ) +
    #   
    #   scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
    #   # scale_x_date(expand = c(0, 0), date_breaks = paste0("1 ", scaleDate)) +
    #   
    #   scale_fill_manual(values = paletteer_d("ggsci::default_igv")) +
    #   
    #   theme_minimal() +
    #   
    #   theme(
    #     panel.grid = element_blank(),
    #     axis.ticks.y = element_line(linewidth = .3),
    #     axis.ticks.x = element_blank(),
    #     axis.text.x = element_blank(),
    #     plot.margin = margin(20, 20, 20, 20)
    #   ) + 
    #   
    #   labs( x = "", y = "" ) + 
    #   facet_grid(as.formula(paste0(".~ ", samplesCluster)), scales = "free_x")
    # 
    # # interactive plotly
    # gr_plotly = ggplotly(
    #   gr, tooltip = c("text", "x"), dynamicTicks = F
    # )
      
      gr_plotly = plot_ly(
          x, x = ~Sample, y = ~Percentage, type = 'bar', 
          color = ~Variant, stroke = ~Variant, strokes = "grey50"
      ) |> 
          layout(
              yaxis = list(ticksuffix = "%"),
              barmode = 'stack', hovermode = "x unified"
          )
  }
  
  # return result
  return(gr_plotly)
}

# Lineplot generation ----------------------------------------------------

freyja_lineplot = function(
    x,    
    from, 
    to, 
    percentageThreshold,
    xAxisSelection,
    samplesCluster = NULL
){ 
  
  # User input selection
  # date range
  x = x[Date >= from & Date <= to]
  
  # MIN percentage threshold for report
  x = x[Percentage >= percentageThreshold/100]
  
  # barplot: Date vs Perc ---------------------------------------------------
  
  # tooltip label 
  # https://stackoverflow.com/questions/34605919/formatting-mouse-over-labels-in-plotly-when-using-ggplotly
  x$tooltip_label = paste0(
    x$Variant, ": ", round(100 * x$Percentage, digits = 3), "%"
  )
  # summarize and generate reactive label
  x = x[, by = Sample, tooltip_label := paste(tooltip_label, collapse = "\n")]
  
  # Graph drawing objective
  if( as.integer(xAxisSelection) == 1 & (samplesCluster == "None" | is.null(samplesCluster) )){
    # 1. Time series
    
    # update tooltip label
    x$tooltip_label = paste(
      x$tooltip_label,
      paste0("Sample: ", x$Sample),
      sep = '\n'
    )
    
    # static ggplot2
    gr = ggplot(
      data = x,
      aes(x = Date,
          y = Percentage#,
          #text = tooltip_label
          )
      ) +
      geom_point(
        aes(fill = Variant), 
        shape = 21, color = "white", size = 2, stroke = .1
      ) +
      
      geom_smooth(aes(color = Variant, fill = Variant)) +
      
      scale_y_continuous(
        # expand = c(0, 0), 
        # limits = c(0, 1), 
        labels = scales::percent) +
      # scale_x_date(
      #   # expand = c(0, 0), 
      #   date_breaks = "1 week") +
      
      scale_color_manual(values = paletteer_d("ggsci::planetexpress_futurama")) +
      scale_fill_manual(values = paletteer_d("ggsci::planetexpress_futurama")) +
      
      theme_minimal() +
      
      theme(
        panel.grid.major = element_line(linewidth = .3, linetype = "dashed", color = "grey75"),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(linewidth = .3),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        plot.margin = margin(20, 20, 20, 20)
      ) + 
      
      labs( x = "", y = "" )

    # interactive plotly
    # gr_plotly = ggplotly(gr, tooltip = c("x", "text"), dynamicTicks = TRUE) 
    
    gr_plotly = ggplotly(
      gr, dynamicTicks = TRUE, tooltip = c("x", "color")
    )
    
    
  } else if( as.integer(xAxisSelection) == 2 & (samplesCluster == "None" | is.null(samplesCluster) )) {
    
    # 2.1. Samples/No clusters
    
    # update tooltip label
    x$tooltip_label = paste(
      x$tooltip_label,
      paste0("Date: ", x$Date),
      sep = '\n'
    )
    
    # static ggplot2
    gr = ggplot(
      data = x,
      aes(
        x = reorder(Sample, Date),
        y = Percentage,
        fill = Variant,
        text = tooltip_label
      )) +
      geom_point(
        shape = 21,
        color = "white",
        size = 2,
        stroke = .1
      ) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 1), labels = scales::percent) +
      # scale_x_date(expand = c(0, 0), date_breaks = "1 day") +
      scale_color_manual(values = paletteer_d("ggsci::planetexpress_futurama")) +
      scale_fill_manual(values = paletteer_d("ggsci::planetexpress_futurama")) +
      theme_minimal() +
      theme(
        # panel.grid.major = element_line(linewidth = .3, linetype = "dashed", color = "grey75"),
        # panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(20, 20, 20, 20)
      ) + 
      labs( x = "", y = "" )
    
    # interactive plotly
    gr_plotly = ggplotly(
      gr, tooltip = c("text", "x"), dynamicTicks = F
    )
    
  } else {
    # 2.2. Samples/Clusters
    
    # update tooltip label
    x$tooltip_label = paste(
      x$tooltip_label,
      paste0("Date: ", x$Date),
      sep = '\n'
    )
    
    # static ggplot2
    gr = ggplot(
      data = x,
      aes(
        x = reorder(Sample, Date),
        y = Percentage,
        fill = Variant,
        text = tooltip_label
      )) +
      geom_point(
        shape = 21,
        color = "white",
        size = 2,
        stroke = .1
      ) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 1), labels = scales::percent) +
      # scale_x_date(expand = c(0, 0), date_breaks = "1 day") +
      scale_color_manual(values = paletteer_d("ggsci::planetexpress_futurama")) +
      scale_fill_manual(values = paletteer_d("ggsci::planetexpress_futurama")) +
      theme_minimal() +
      theme(
        # panel.grid.major = element_line(linewidth = .3, linetype = "dashed", color = "grey75"),
        # panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(20, 20, 20, 20)
      ) + 
      labs( x = "", y = "" )+
      facet_grid(as.formula(paste0(".~ ", samplesCluster)), scales = "free_x")
    
    # interactive plotly
    gr_plotly = ggplotly(
      gr, tooltip = c("text", "x"), dynamicTicks = F
    )
  }
  
  # return result
  return(gr_plotly)
}
