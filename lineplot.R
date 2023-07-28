

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


gr = ggplot(data = list1, aes(x = pseudo_date, y = Prop)) +
    
    geom_point(
        aes(fill = Variant), 
        shape = 21, color = "white", size = 2, stroke = .1
    ) +
    
    geom_smooth(aes(color = Variant, fill = Variant)) +
    
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1), labels = scales::percent) +
    scale_x_date(expand = c(0, 0), date_breaks = "1 month") +
    
    scale_color_manual(values = paletteer_d("ggsci::planetexpress_futurama")) +
    scale_fill_manual(values = paletteer_d("ggsci::planetexpress_futurama")) +
    
    theme_minimal() +
    
    theme(
        panel.grid.major = element_line(linewidth = .3, linetype = "dashed", color = "grey75"),
        panel.grid.minor = element_blank(),
        
        axis.ticks = element_line(linewidth = .3),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        
        plot.margin = margin(20, 20, 20, 20)
    )



ggsave(
    plot = gr, filename = "lineplot.jpeg",
    width = 12, height = 6, units = "in"
)


gr_plotly = ggplotly(
    gr, dynamicTicks = TRUE, tooltip = c("x", "color")
)

htmlwidgets::saveWidget(
    widget = gr_plotly,
    file = "lineplot.html", 
    selfcontained = TRUE 
)


gr2 = ggplot(data = list2, aes(x = pseudo_date, y = abundance)) +
    
    geom_point(
        aes(fill = lineage), 
        shape = 21, color = "white", size = 1.2, stroke = .1
    ) +
    
    # geom_smooth(aes(color = lineage, fill = lineage)) +
    
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1), labels = scales::percent) +
    scale_x_date(expand = c(0, 0), date_breaks = "1 month") +
    
    scale_color_manual(values = paletteer_c("viridis::plasma", 353)) +
    scale_fill_manual(values = paletteer_c("viridis::plasma", 353)) +
    
    theme_minimal() +
    
    theme(
        panel.grid.major = element_line(linewidth = .3, linetype = "dashed", color = "grey75"),
        panel.grid.minor = element_blank(),
        
        axis.ticks = element_line(linewidth = .3),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        
        plot.margin = margin(20, 20, 20, 20)
    )



gr_plotly = ggplotly(
    gr2, dynamicTicks = TRUE, tooltip = c("x", "fill")
)

htmlwidgets::saveWidget(
    widget = gr_plotly,
    file = "dotplot.html", 
    selfcontained = TRUE 
)















