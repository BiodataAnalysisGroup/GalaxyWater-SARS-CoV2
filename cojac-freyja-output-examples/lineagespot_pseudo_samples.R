library(lineagespot)
library(data.table)
# 
results <- lineagespot(vcf_folder = system.file("extdata", "vcf-files", 
                                                package = "lineagespot"),
                       
                       gff3_path = system.file("extdata", 
                                               "NC_045512.2_annot.gff3", 
                                               package = "lineagespot"),
                       
                       ref_folder = system.file("extdata", "ref", 
                                                package = "lineagespot"))
#
dim(results$variants.table)
View(results$variants.table)

# 
dim(results$lineage.hits)
View(results$lineage.hits)

#
dim(results$lineage.report)
View(results$lineage.report)

# 
fwrite(
  results$lineage.report,
  file = 'lineagespot_testing_data.tsv',
  row.names = F,
  col.names = T,
  quote = F,
  sep = '\t'
)
# 
pseudo_dates = seq(
  from = as.Date("2020-01-01"), by = "month",
  length.out = results$lineage.report$sample |> unique() |> length()
)
#
df1 = results$lineage.report
df1 = df1[order(df1$sample),]
df1$pseudo_date = pseudo_dates[ df1$sample |> rleid() ]
df1$pseudo_clusters = LETTERS[1:3][ df1$sample |> rleid() ]
# 
fwrite(
  df1[,c('sample','pseudo_date', 'pseudo_clusters')],
  file = 'lineagespot_testing_metadata.tsv',
  row.names = F,
  col.names = T,
  quote = F,
  sep = '\t'
)


