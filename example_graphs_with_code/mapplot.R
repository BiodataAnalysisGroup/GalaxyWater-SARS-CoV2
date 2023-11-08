

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

library(maps)
library(sf)

library(paletteer)

pseudo_cities = us.cities[ list1$Sample |> rleid(), ] |> setDT()


list1 = list1 |> cbind(pseudo_cities)


world1 <- st_as_sf(map('world', plot = FALSE, fill = TRUE))
states <- st_as_sf(map('state', plot = FALSE, fill = TRUE))

gr = ggplot() + 
    
    geom_sf(data = world1, fill = "grey96", color = NA) +
    geom_sf(data = states, fill = NA, color = "grey10", linewidth = .15) +
    
    geom_point(data = list1, 
               aes(x = long, y = lat, fill = Variant, size = Prop, text = name),
               shape = 21, stroke = .1, color = "white", alpha = .5) +
    
    scale_fill_manual(
        values = paletteer_d("ggsci::planetexpress_futurama"),
        guide = guide_legend(override.aes = list(size = 3))
    ) +
    
    scale_size_continuous(
        range = c(1, 5),
        limits = c(0, 1),
        breaks = c(.25, .5, .75, 1),
        labels = scales::percent,
        guide = guide_legend(
            override.aes = list(color = "black", linewidth = .25)
        )
    ) +
    
    facet_wrap(vars(Variant), ncol = 2) +
    
    coord_sf(xlim = c(-130, -65), ylim = c(25, 50)) +
    
    theme_void() +
    
    theme(
        strip.text = element_text(margin = margin(b = 5), face = "bold"),
        panel.spacing = unit(2, "lines"),
        
        plot.margin = margin(10, 10, 10, 10)
    )


ggsave(
    plot = gr, filename = "mapplot.jpeg", dpi = 600,
    width = 10, height = 8, units = "in"
)



gr_plotly = ggplotly(
    gr, dynamicTicks = TRUE
)

htmlwidgets::saveWidget(
    widget = gr_plotly,
    file = "mapplot.html", 
    selfcontained = TRUE 
)






