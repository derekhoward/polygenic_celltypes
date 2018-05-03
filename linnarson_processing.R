library(readr)
library(dplyr)
library(magrittr)

table <- read_tsv('./l5_all.agg.tab', col_names = FALSE)

table %<>% filter(!is.na(X1) | X8 == 'ClusterName')
colnames(table) <- c(unlist(table[2,1:7]), unlist(table[1, 8:ncol(table)]))
table <- table[-c(1,2),]

# columns to drop
dropcols <- c('Accession', '_LogCV', '_LogMean', '_Selected', '_Total', '_Valid', 'ClusterName', NA)

table %<>% select(-one_of(dropcols))

#to change all but the 'Gene' column name into numeric type
table[,2:length(colnames(table))] %<>% lapply(function(x) as.numeric(as.character(x)))

#now melt
linnarsson <- reshape2::melt(table, id.vars='Gene', variable.name='cluster_id', value.name='expression')
linnarsson <- mutate(linnarsson, log1Expression=log(1+expression))

# zscore across genes
linnarsson %<>% 
  group_by(Gene) %>% 
  mutate(log1ExpZ = (log1Expression - mean(log1Expression)) / sd(log1Expression))

print(dim(linnarsson))
linnarsson %<>% filter(!is.na(log1ExpZ))       
print(dim(linnarsson))
linnarsson %<>% select(-expression, -log1Expression)
save(linnarsson, file='processed_zeisel.Rdata')


