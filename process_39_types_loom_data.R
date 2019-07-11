library(dplyr)
library(loomR)
library(here)
library(stringr)
library(purrr)
library(tidyr)
library(magrittr)
library(homologene)

lfile <- connect(filename = here("l6_r4.agg.loom"), mode = "r+")
lfile

lfile[["col_graphs"]]

exp_matrix <- as_tibble(t(as.matrix(lfile[["matrix"]][,])))
colnames(exp_matrix) <- lfile$col.attrs$TaxonomyRank4[]
exp_matrix %<>% mutate(Gene = lfile$row.attrs$Gene[]) %>% select(Gene, everything())
exp_matrix
sum(duplicated(exp_matrix$Gene))
lfile$close_all()

#matrix has duplicate gene symbols -maybe just average them for now
exp_matrix %<>% group_by(Gene) %>% summarise_all(mean)

exp_matrix %<>% gather('cell_type', 'expression', -Gene) 
exp_matrix %<>% mutate(log1Expression=log(1+expression))

# zscore across genes
exp_matrix %<>% 
  group_by(Gene) %>% 
  mutate(log1ExpZ = (log1Expression - mean(log1Expression)) / sd(log1Expression))

dim(exp_matrix)
exp_matrix %<>% filter(!is.na(log1ExpZ))
dim(exp_matrix)

exp_matrix %<>% select(-expression, -log1Expression)

# should we avg the gene expression here for duplicate genes?
#linnarsson %<>% group_by(Gene, cluster_id) %>% summarize(log1ExpZ = mean(log1ExpZ))

exp_matrix %<>% group_by(cell_type) %<>%
  mutate(log1ExpZRank = rank(log1ExpZ)) %<>%
  select(-log1ExpZ)

mouseMatrix_39_celltypes <- exp_matrix %>% spread(cell_type, log1ExpZRank)
mouseMatrix_39_celltypes %<>% as.data.frame()
rownames(mouseMatrix_39_celltypes) <- mouseMatrix_39_celltypes$Gene
mouseMatrix_39_celltypes$Gene <- NULL

# we also want a matrix with only human reachable genes
unique_genes_all <- rownames(mouseMatrix_39_celltypes)
unique_genes_human_reachable <- mouse2human(unique_genes_all)$mouseGene

humanReachableMatrix_39_celltypes <- exp_matrix %>% filter(Gene %in% unique_genes_human_reachable)
humanReachableMatrix_39_celltypes %<>% mutate(log1ExpZRank = rank(log1ExpZRank)) #rerank after filtering
humanReachableMatrix_39_celltypes %<>% spread(cell_type, log1ExpZRank) 
humanReachableMatrix_39_celltypes %<>% as.data.frame()
rownames(humanReachableMatrix_39_celltypes) <- humanReachableMatrix_39_celltypes$Gene
humanReachableMatrix_39_celltypes$Gene <- NULL


cluster_id <- colnames(mouseMatrix_39_celltypes)

cluster_info <- cluster_id %>% 
  str_to_lower() %>% 
  str_replace_all('-', ' ') %>% 
  str_replace_all('  ', ' ') %>%
  str_replace_all(' ', '-') %>% 
  map_chr(~ paste0('r4_', .))

cell_type_descriptions_39 <- as_tibble(cbind(cluster_id, cluster_info))


save(mouseMatrix_39_celltypes, humanReachableMatrix_39_celltypes, cell_type_descriptions_39, file='processed_39_cellypes.Rdata')
