library(data.table)

d1 = fread("cojac-freyja-output-examples/ca-freyja-tsv.tsv")

any(duplicated(d1$V1))
# FALSE

pseudo_dates = seq(
  from = as.Date("2020-01-01"), 
  by = "day",
  length.out = d1 |> nrow()
)

pseudo_clusters = c(
  rep('A', d1 |> nrow()/3),
  rep('B', d1 |> nrow()/3),
  rep('C', d1 |> nrow()/3)
)
pseudo_clusters = sample(x = pseudo_clusters, size = length(pseudo_clusters), replace = F)

# pseudo_metadata
pseudo_metadata = data.frame(
  Sample = d1$V1,
  Date = pseudo_dates,
  pseudo_clusters = pseudo_clusters
)

fwrite("cojac-freyja-output-examples/ca-freyja-tsv.tsv")

# save as TSV
fwrite(
  pseudo_metadata,
  file = "cojac-freyja-output-examples/ca-freyja-tsv-pseudometadata.tsv",
  row.names = F,
  col.names = T,
  quote = F,
  sep = '\t'
)
