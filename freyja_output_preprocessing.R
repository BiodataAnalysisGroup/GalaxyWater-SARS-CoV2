# load libraries --------------------
library(data.table)
library(stringr)
library(ggplot2)
library(plotly)
library(paletteer)
# for colors https://r-charts.com/color-palettes/

# Freyja pipeline results -------------------------------------------------

freyja_results_preprocessing = function(x){
  
  # preprocess Freyja results -----------------------------------------------

  # load file
  freyja_res = fread(x)
  
  # # rename first column
  # colnames(freyja_res)[1] = "Sample"
  
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
  
  # Generate artificial metadata (to be removed) ----------------------------
  
  # # generate pseudo dates
  # pseudo_dates = seq(
  #   from = as.Date("2020-01-01"),
  #   by = "day",
  #   length.out = freyja_summarized$Sample |> unique() |> length()
  # )
  # 
  # # add pseudo dates to the results
  # freyja_summarized$pseudo_date = pseudo_dates[ freyja_summarized$Sample |> rleid() ]
  # freyja_lineages$pseudo_date = pseudo_dates[ freyja_lineages$Sample |> rleid() ]
  # 
  # # generate pseudo clusters
  # pseudo_clusters = c(
  #   rep("ClusterA", 100),
  #   rep("ClusterB", 100),
  #   rep( "ClusterC", length(unique(freyja_summarized$Sample))-200 )
  # )
  # # and mix
  # pseudo_clusters = sample(x = pseudo_clusters, size = length(pseudo_clusters), replace = FALSE)
  # 
  # # add pseudo dates to the results
  # freyja_summarized$pseudo_clusters = pseudo_clusters[ freyja_summarized$Sample |> rleid() ]
  # freyja_lineages$pseudo_clusters = pseudo_clusters[ freyja_lineages$Sample |> rleid() ]
  
  # Index for additional (unknown) columns
  idx = !colnames(freyja_res) %in% c("Sample",
                                     "summarized",
                                     "lineages",
                                     "abundances",
                                     "resid",
                                     "coverage")
  
  # match with Sample column (IDs) and add to summarized and lineages DTs
  for(i in seq(ncol(freyja_res[,idx, with = FALSE]))){
    # column name
    colName = colnames(freyja_res[,idx, with = FALSE])[i]
    # column values
    colValues = freyja_res[,idx, with = FALSE][,i, with = FALSE]
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

freyja_results_graph = function(
    x,    
    from, 
    to, 
    percentageThreshold,
    xAxisSelection,
    samplesCluster = NULL
  ){ 
  
  # User input selection
  # date range
  x$summarized = x$summarized[Date >= from & Date <= to]
  
  # MIN percentage threshold for report
  x$summarized = x$summarized[Percentage >= percentageThreshold/100]
  
  # barplot: Date vs Perc ---------------------------------------------------

  # tooltip label 
  # https://stackoverflow.com/questions/34605919/formatting-mouse-over-labels-in-plotly-when-using-ggplotly
  x$summarized$tooltip_label = paste0(
    x$summarized$Variant, ": ", round(100 * x$summarized$Percentage, digits = 3), "%"
  )
  # summarize and generate reactive label
  x$summarized = x$summarized[, by = Sample, tooltip_label := paste(tooltip_label, collapse = "\n")]
  
  # Graph drawing objective
  if( as.integer(xAxisSelection) == 1 & (samplesCluster == "None" | is.null(samplesCluster) )){
    # 1. Time series
    
    # update tooltip label
    x$summarized$tooltip_label = paste(
      x$summarized$tooltip_label,
      paste0("Sample: ", x$summarized$Sample),
      sep = '\n'
    )
    
    # static ggplot2
    gr = ggplot(data = x$summarized, aes(text = tooltip_label)) +
      
      geom_bar(
        aes(x = Date, y = Percentage, fill = Variant),
        position = "fill", # position_dodge2(preserve = "single"), # "dodge", #
        stat = "identity"
      ) +
      
      scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
      # scale_x_date(
      #   expand = c(0, 0), 
      #   date_breaks = paste0("1 ", scaleDate)
      #   ) +
      
      scale_fill_manual(values = paletteer_d("ggsci::default_igv")) +
      
      theme_minimal() +
      
      theme(
        panel.grid = element_blank(),
        
        axis.ticks = element_line(linewidth = .3),
        axis.text.x = element_text(
          size = 8,
          angle = 45,
          hjust = 1
        ),
        
        plot.margin = margin(20, 20, 20, 20)
      ) + 
      
      labs( x = "", y = "" )
    
    # interactive plotly
    gr_plotly = ggplotly(gr, tooltip = c("text", "x"), dynamicTicks = TRUE) 
    
  } else if( as.integer(xAxisSelection) == 2 & (samplesCluster == "None" | is.null(samplesCluster) )) {
    
    # 2.1. Samples/No clusters

    # update tooltip label
    x$summarized$tooltip_label = paste(
      x$summarized$tooltip_label,
      paste0("Date: ", x$summarized$Date),
      sep = '\n'
    )
    
    # static ggplot2
    gr = ggplot(data = x$summarized, aes(text = tooltip_label)) +
      
      geom_bar(
        aes(x = reorder(Sample, Date), y = Percentage, fill = Variant),
        position = "fill",
        stat = "identity"
      ) +
      
      scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
      # scale_x_date(expand = c(0, 0), date_breaks = paste0("1 ", scaleDate)) +
      
      scale_fill_manual(values = paletteer_d("ggsci::default_igv")) +
      
      theme_minimal() +
      
      theme(
        panel.grid = element_blank(),
        axis.ticks.y = element_line(linewidth = .3),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(20, 20, 20, 20)
      ) + 
      
      labs( x = "", y = "" ) # +
      # facet_grid(as.formula(paste0(".~ ", samplesCluster)), scales = "free_x")
    
    # interactive plotly
    gr_plotly = ggplotly(
      gr, tooltip = c("text", "x"), dynamicTicks = F
    )
    
  } else {
    # 2.2. Samples/Clusters
    
    # update tooltip label
    x$summarized$tooltip_label = paste(
      x$summarized$tooltip_label,
      paste0("Date: ", x$summarized$Date),
      sep = '\n'
    )
    
    # static ggplot2
    gr = ggplot(data = x$summarized, aes(text = tooltip_label)) +
      
      geom_bar(
        aes(x = reorder(Sample, Date), y = Percentage, fill = Variant),
        position = "fill",
        stat = "identity"
      ) +
      
      scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
      # scale_x_date(expand = c(0, 0), date_breaks = paste0("1 ", scaleDate)) +
      
      scale_fill_manual(values = paletteer_d("ggsci::default_igv")) +
      
      theme_minimal() +
      
      theme(
        panel.grid = element_blank(),
        axis.ticks.y = element_line(linewidth = .3),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(20, 20, 20, 20)
      ) + 
      
      labs( x = "", y = "" ) + 
      facet_grid(as.formula(paste0(".~ ", samplesCluster)), scales = "free_x")
    
    # interactive plotly
    gr_plotly = ggplotly(
      gr, tooltip = c("text", "x"), dynamicTicks = F
    )
  }
  
  # lineplot: Date vs Perc ---------------------------------------------------
  # lineplot: Date vs Perc ---------------------------------------------------
  
  # return result
  return(gr_plotly)
}
