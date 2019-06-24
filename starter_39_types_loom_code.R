library(dplyr)
library(loomR)
library(here)


lfile <- connect(filename = here("l6_r4.agg.loom"), mode = "r+")
lfile

lfile[["col_graphs"]]

exp_matrix <- as_tibble(t(as.matrix(lfile[["matrix"]][,])))
colnames(exp_matrix) <- lfile$col.attrs$TaxonomyRank4[]
exp_matrix %<>% mutate(Gene = lfile$row.attrs$Gene[]) %>% select(Gene, everything())
exp_matrix
duplicated(exp_matrix$Gene)
#matrix has duplicate gene symbols -maybe just average them for now

