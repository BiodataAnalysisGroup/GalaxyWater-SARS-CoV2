# load libraries
library(data.table)
library(ggplot2)
library(plotly)
library(dplyr)
library(stringr)
library(scales)

# Load and process Freyja example results ---------------------------------
freyja_results_us = fread("us-freyja-tsv.tsv")
colnames(freyja_results_us)

# # arbitrary min threshold (%) for summarization of covid lineages
# min_threshold 

# # split lineages and abundances
# covid_lineages = str_split(
#   string = freyja_results_us$lineages,
#   pattern = " "
# )
# covid_abundances = str_split(
#   string = freyja_results_us$abundances,
#   pattern = " "
# )
# all(unlist(lapply(d1, length)) == unlist(lapply(d2, length)))
# # TRUE

# add pseudo dates (7 pseudo-weeks) to the Freyja example dataset:
pseudo_dates = seq(
  from = as.Date("2020-03-01"),
  to = as.Date("2020-04-30"),
  by = "week"
)[1:7] %>% as.character

# ~50 samples per pseudo-week
pseudo_dates_split = split(
  x = 1:nrow(freyja_results_us),
  ceiling(seq_along(1:nrow(freyja_results_us))/50)
)

# add pseudo dates 
freyja_results_us$date = ""
for(i in seq(length(pseudo_dates_split))){
  freyja_results_us$date[pseudo_dates_split[[i]]] = pseudo_dates[i]
}
freyja_results_us$date = as.Date(freyja_results_us$date)

# emtpy list
freyja_results_us_long = list()
# iterate
for(i in seq(nrow(freyja_results_us))){
  # split summarized lineages and abundances
  # split summarized observation
  firstSplit = freyja_results_us$summarized[i] %>%
    sub(pattern = "^[[]", replacement = "") %>%
    sub(pattern = "[]]$", replacement = "") %>%
    str_split(pattern = "[)],") %>%
    unlist()
  # split further and convert to data.frame
  # extract lineage
  sum_lineage = str_extract(string = firstSplit, pattern = "'.*?'") %>%
    gsub(pattern = "[']", replacement = "")
  # extract abundance
  sum_abundance = str_extract(string = firstSplit, pattern = ", .*?$") %>%
    gsub(pattern = "[, ]", replacement = "") %>%
    gsub(pattern = "[)]", replacement = "") %>%
    as.numeric()
  # form data.frame and save in list
  freyja_results_us_long[[i]] = data.frame(
    sample = freyja_results_us$V1[i],
    lineage = sum_lineage,
    abundance = sum_abundance, # turn to %
    resid = freyja_results_us$resid[i],
    coverage = freyja_results_us$coverage[i],
    date = freyja_results_us$date[i]
    )
}
# convert to data.frame
freyja_results_us_long = do.call(
  rbind, freyja_results_us_long
)

# Stacked bar chart for SARS-CoV-2 variant lineage (%) per sample ---------

# Example plot of the first 100 samples
freyja_results_us_long_subset = subset(
  freyja_results_us_long, sample %in% unique(freyja_results_us_long$sample)[1:100]
)
# stacked bar graph
tiff(
  filename = "stacked_bar_graph_Freyja_us_100_samples.tif",
  units = "in",
  width = 12,
  height = 8,
  compression = "lzw",
  res = 300
)
ggplot(
  data = freyja_results_us_long_subset,
  aes(
    x = sample,
    y = abundance,
    fill = lineage
  )
)+
  geom_bar(stat="identity") +
  labs(
    title = "Freyja US example - 100 samples",
    x = "", 
    y = "SARS-CoV-2 variant lineage abundance (%)"
    )+
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 8),
    axis.title.y = element_text(margin = margin(r = 10, t = 0, l = 0, b = 0))
    ) + 
  scale_y_continuous(labels = percent) + 
  guides(fill=guide_legend(title="Lineage\n(summarized)"))
dev.off()


# Stacked bar chart for SARS-CoV-2 variant lineage (%) per date -----------

# convert to data.table
freyja_results_us_long_dt = setDT(freyja_results_us_long)

# Total sum of abundance by date
total_abundance = freyja_results_us_long_dt[, .(abundance = sum(abundance)), by = date]
# Sum of abundance by date and lineage
lineage_abundance = freyja_results_us_long_dt[, .(abundance = sum(abundance)), by = .(lineage, date)]
# Copy data.table for updates
lineage_abundance_2 = lineage_abundance
# Replace abundance values (abundance % of lineage per samples of specific date)
for(i in seq(length(total_abundance$date))){
  for(j in lineage_abundance$lineage){
    # index
    idx = which(lineage_abundance$date == total_abundance$date[i] & lineage_abundance$lineage == j)
    # update data.table
    lineage_abundance_2$abundance[idx] = lineage_abundance$abundance[idx]/total_abundance$abundance[i]
  }
}

# stacked bar graph
tiff(
  filename = "stacked_bar_graph_Freyja_us_pseudo_dates.tif",
  units = "in",
  width = 12,
  height = 8,
  compression = "lzw",
  res = 300
)
# change to English
Sys.setlocale("LC_ALL", "English")
# plot
ggplot(
  data = lineage_abundance_2,
  aes(
    x = date,
    y = abundance,
    fill = lineage
  )
)+
  geom_bar(stat="identity") +
  labs(
    title = "Freyja US example - pseudo-dates",
    x = "", 
    y = "SARS-CoV-2 variant lineage abundance (%)"
  )+
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1, size = 8),
    axis.title.y = element_text(margin = margin(r = 10, t = 0, l = 0, b = 0))
  ) + 
  scale_y_continuous(labels = percent) + 
  guides(fill=guide_legend(title="Lineage\n(summarized)"))
dev.off()
