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

freyja_barplot = function(x,
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
      color = ~ Variant,
      stroke = ~ Variant,
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
        x = ~ reorder(Sample, Date),
        y = ~ Percentage,
        type = 'bar',
        color = ~ Variant,
        stroke = ~ Variant,
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
        x = ~ Sample,
        y = ~ Percentage,
        type = 'bar',
        color = ~ Variant,
        stroke = ~ Variant,
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
        mutate_at("Variant", as.factor) %>%
        group_by(get(samplesCluster)) %>%
        group_map(# .keep = TRUE,
          .f = ~ {
            tidyr::complete(data = .x, Variant) %>%
              plot_ly(
                data = .,
                x = ~ reorder(Sample, Date),
                y = ~ Percentage,
                type = 'bar',
                color = ~ Variant,
                stroke = ~ Variant,
                strokes = "grey50",
                showlegend = FALSE, # (.y == "B"),
                legendgroup = ~ Variant
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
        mutate_at("Variant", as.factor) %>%
        group_by(get(samplesCluster)) %>%
        group_map(
          .keep = TRUE,
          .f = ~ {
            tidyr::complete(data = .x, Variant) %>%
              plot_ly(
                data = .,
                x = ~ Sample,
                y = ~ Percentage,
                type = 'bar',
                color = ~ Variant,
                stroke = ~ Variant,
                strokes = "grey50",
                showlegend = FALSE, # (.y == "B"),
                legendgroup = ~ Variant
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
