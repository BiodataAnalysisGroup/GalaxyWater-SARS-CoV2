# lineagespot pipeline results -------------------------------------------------

# lineagespot pipeline output preprocessing ------------------------------------

lineagespot_results_preprocessing = function(data, metadata){
  
  # merge input Data and Metadata based on Sample IDs 
  
  # order based on Data Sample IDs
  idx_sample = match(data$sample, metadata$sample)
  metadata = metadata[idx_sample,]
  
  # cbind
  lineagespot_res = cbind(
    data,
    metadata[,-which(colnames(metadata) == 'sample'),with = F]
  )
  
  # rename date column to 'Date'
  if(
    grepl( fixed = T, pattern = 'date', x = tolower(  colnames(lineagespot_res))) |> any()
  ){
    date_idx = grep( fixed = T, pattern = 'date', x = tolower(  colnames(lineagespot_res)))
    colnames(lineagespot_res)[date_idx] = 'Date'
  }
  
  # percenatages list
  perc_list = list()
  
  # adjust lineage prop. to percentage (%)
  for(i in unique(lineagespot_res$sample) |> length() |> seq()){
    # subset
    perc_values = subset(lineagespot_res, sample == unique(lineagespot_res$sample)[i])$`lineage prop.`
    lineages = subset(lineagespot_res, sample == unique(lineagespot_res$sample)[i])$lineage
    # convert to percentages
    perc_values = 100*perc_values/sum(perc_values)
    # save
    perc_list[[i]] = data.frame(
      sample = unique(lineagespot_res$sample)[i],
      lineage = lineages,
      percentages = perc_values
    )
  }
  # to data.frame
  perc_list = do.call(rbind, perc_list)
  
  # match and add percentages 
  idx_match = match(
    paste0(lineagespot_res$sample, "_", lineagespot_res$lineage),
    paste0(perc_list$sample, "_", perc_list$lineage)
  )
  lineagespot_res$Percentage = perc_list[idx_match,]$percentages
  
  # return
  return(
    lineagespot_res
  )
}

# lineagespot graph generation -------------------------------------------------

# Barplot generation ----------------------------------------------------

lineagespot_barplot = function(x,
                          from,
                          to,
                          percentageThreshold,
                          xAxisSelection,
                          samplesCluster = NULL) {
  # User input selection
  # date range
  x = x[Date >= from & Date <= to]
  
  # MIN percentage threshold for report
  x = x[Percentage >= percentageThreshold / 100]
  
  # 100*Perc
  x$Percentage = x$Percentage * 100
  
  # Barplot drawing
  if (as.integer(xAxisSelection) == 1 &
      (samplesCluster == "None" | is.null(samplesCluster))) {
    # 1. Time series/No clusters
    
    # plot
    gr_plotly = plot_ly(
      x,
      x = ~ Date,
      y = ~ Percentage,
      type = 'bar',
      color = ~ lineage,
      stroke = ~ lineage,
      strokes = "grey50"
    ) |>
      layout(
        yaxis = list(ticksuffix = "%", title = ""),
        xaxis = list(showticklabels = TRUE, title = ""),
        barmode = 'stack',
        hovermode = "x unified"
      )
    
  } else if (as.integer(xAxisSelection) == 2 &
             (samplesCluster == "None" | is.null(samplesCluster))) {
    # 2.1. Samples/No clusters
    
    # reorder by date
    if('Date' %in% colnames(x)){
      # plot
      gr_plotly = plot_ly(
        x,
        x = ~ reorder(sample, Date),
        y = ~ Percentage,
        type = 'bar',
        color = ~ lineage,
        stroke = ~ lineage,
        strokes = "grey50"
      ) |>
        layout(
          yaxis = list(ticksuffix = "%", title = ""),
          xaxis = list(showticklabels = FALSE, title = ""),
          barmode = 'stack',
          hovermode = "x unified"
        )
    } else {
      # plot
      gr_plotly = plot_ly(
        x,
        x = ~ sample,
        y = ~ Percentage,
        type = 'bar',
        color = ~ lineage,
        stroke = ~ lineage,
        strokes = "grey50"
      ) |>
        layout(
          yaxis = list(ticksuffix = "%", title = ""),
          xaxis = list(showticklabels = FALSE, title = ""),
          barmode = 'stack',
          hovermode = "x unified"
        )
    }
    
  } else {
    # 2.2. Samples/Clusters
    
    # reorder by date
    if('Date' %in% colnames(x)){
      # plot
      gr_plotly = x %>%
        mutate_at("lineage", as.factor) %>%
        group_by(get(samplesCluster)) %>%
        group_map(# .keep = TRUE,
          .f = ~ {
            tidyr::complete(data = .x, lineage) %>%
              plot_ly(
                data = .,
                x = ~ reorder(sample, Date),
                y = ~ Percentage,
                type = 'bar',
                color = ~ lineage,
                stroke = ~ lineage,
                strokes = "grey50",
                showlegend = FALSE, # (.y == "B"),
                legendgroup = ~ lineage
              ) |>
              layout(
                yaxis = list(ticksuffix = "%", title = ""),
                xaxis = list(showticklabels = FALSE, title = unique(.x[[samplesCluster]])),
                barmode = 'stack',
                hovermode = "x unified"
              )
          }) %>%
        subplot(nrows = 1,
                shareX = TRUE,
                shareY = TRUE,
                titleX = TRUE)
    } else {
      # plot
      gr_plotly = x %>%
        mutate_at("lineage", as.factor) %>%
        group_by(get(samplesCluster)) %>%
        group_map(
          .keep = TRUE,
          .f = ~ {
            tidyr::complete(data = .x, lineage) %>%
              plot_ly(
                data = .,
                x = ~ sample,
                y = ~ Percentage,
                type = 'bar',
                color = ~ lineage,
                stroke = ~ lineage,
                strokes = "grey50",
                showlegend = FALSE, # (.y == "B"),
                legendgroup = ~ lineage
              ) |>
              layout(
                yaxis = list(ticksuffix = "%", title = ""),
                xaxis = list(showticklabels = FALSE, title = unique(.x[[samplesCluster]])),
                barmode = 'stack',
                hovermode = "x unified"
              )
          }) %>%
        subplot(nrows = 1,
                shareX = TRUE,
                shareY = TRUE,
                titleX = TRUE)
    }
    
  }
  
  # return result
  return(gr_plotly)
}
