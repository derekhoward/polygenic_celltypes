library(readr)
library(dplyr)
library(magrittr)
getwd()

table <- read_tsv('./l5_all.agg.tab', col_names = FALSE)
descriptions <- table[29, 9:273]

table %<>% filter(!is.na(X1) | X8 == 'ClusterName')
colnames(table) <- c(unlist(table[2,1:7]), unlist(table[1, 8:ncol(table)]))
table <- table[-c(1,2),]

# columns to drop
dropcols <- c('Accession', '_LogCV', '_LogMean', '_Selected', '_Total', '_Valid', 'ClusterName', NA)
table %<>% select(-one_of(dropcols))

#to change all but the 'Gene' column name into numeric type
table[,2:length(colnames(table))] %<>% lapply(function(x) as.numeric(as.character(x)))

cell_type_info <- rbind(colnames(table[2:266]), descriptions)
cell_type_info <- t(cell_type_info)
colnames(cell_type_info) <- c('cluster_id', 'description')
write.csv(cell_type_info, file = "celltype_descriptions.csv", row.names=FALSE)
#now melt
linnarsson <- as_tibble(reshape2::melt(table, id.vars='Gene', variable.name='cluster_id', value.name='expression'))
linnarsson %<>% mutate(log1Expression=log(1+expression))

# zscore across genes
linnarsson %<>% 
  group_by(Gene) %>% 
  mutate(log1ExpZ = (log1Expression - mean(log1Expression)) / sd(log1Expression))

print(dim(linnarsson))
linnarsson %<>% filter(!is.na(log1ExpZ))       
print(dim(linnarsson))
linnarsson %<>% select(-expression, -log1Expression)

#some genes are duplicated - just average them
linnarsson %<>% group_by(Gene, cluster_id) %>% summarize(
  log1ExpZ = mean(log1ExpZ)
)

linnarsson %<>% group_by(cluster_id)
linnarsson %<>% mutate(log1ExpZRank = rank(log1ExpZ)) %>% select(-log1ExpZ)

linnarssonMatrix <- dcast(linnarsson, Gene ~ cluster_id, value.var = "log1ExpZRank")
rownames(linnarssonMatrix) <- linnarssonMatrix$Gene
linnarssonMatrix$Gene <- NULL

save(linnarssonMatrix, file='processed_zeisel.Rdata')


