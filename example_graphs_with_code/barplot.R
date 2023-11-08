

rm(list = ls())
gc()

# load libraries --------------------
library(data.table)
library(stringr)

library(ggplot2)
library(plotly)

# library(dplyr)
# library(scales)

# Load and process Freyja example results ---------------------------------

df = "us-freyja-tsv.tsv" |> fread()

colnames(df)[1] = "Sample"


df$Sample = str_split(df$Sample, "\\.", simplify = TRUE)[, 1]

list1 = df$summarized |>
    str_remove_all("\\[|\\]") |>
    str_split("\\)\\,\\ \\(") |>
    lapply(str_remove_all, "\\(|\\)") |>
    lapply(function(x) {
        data.table(
            "Variant" = str_split(x, "\\,", simplify = TRUE)[, 1] |> str_remove_all("\\'"),
            "Prop"    = str_split(x, "\\,", simplify = TRUE)[, 2] |> as.numeric()
        )
    })

names(list1) = df$Sample



list2 = mapply(function(x, y) {
    
    data.table(
        "lineage"   = x,
        "abundance" = y |> as.numeric()
    )
    
}, str_split(df$lineages, "\\ "), str_split(df$abundances, "\\ "), SIMPLIFY = FALSE)

names(list2) = df$Sample



list1 = list1 |> rbindlist(idcol = "Sample")
list2 = list2 |> rbindlist(idcol = "Sample")




# artificical metadata -----------------------------


list1$Sample |> unique() |> length()

pseudo_dates = seq(
    from = as.Date("2020-01-01"), by = "day",
    length.out = list1$Sample |> unique() |> length()
)


list1$pseudo_date = pseudo_dates[ list1$Sample |> rleid() ]
list2$pseudo_date = pseudo_dates[ list2$Sample |> rleid() ]


# graph ----------------------------------

# for colors https://r-charts.com/color-palettes/

library(paletteer)

list1 = list1[order(pseudo_date, -Prop), ]

# tooltip label 
# https://stackoverflow.com/questions/34605919/formatting-mouse-over-labels-in-plotly-when-using-ggplotly

list1$tooltip_label = paste0(
    list1$Variant, ": ", round(100 * list1$Prop, digits = 3), "%"
)
    
list1 = list1[, by = Sample, tooltip_label := paste(tooltip_label, collapse = "\n")]
    
    
gr = ggplot(data = list1, aes(text = tooltip_label)) +
    
    geom_bar(aes(x = pseudo_date, y = Prop, fill = Variant), 
             position = "fill", stat = "identity") +
    
    scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
    scale_x_date(expand = c(0, 0), date_breaks = "1 month") +
    
    scale_fill_manual(values = paletteer_d("ggsci::default_igv")) +
    
    theme_minimal() +
    
    theme(
        
        panel.grid = element_blank(),
        
        axis.ticks = element_line(linewidth = .3),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        
        plot.margin = margin(20, 20, 20, 20)
    )


ggsave(
    plot = gr, filename = "barplot.jpeg",
    width = 12, height = 6, units = "in"
)


gr_plotly = ggplotly(
    gr, tooltip = c("text", "x"), dynamicTicks = TRUE
)

htmlwidgets::saveWidget(
    widget = gr_plotly,
    file = "barplot.html", 
    selfcontained = TRUE 
)
