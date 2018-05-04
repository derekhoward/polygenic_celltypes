library(ggplot2)
library(grid)
source("./AUCFunction.R")
library(reshape2)
library(magrittr)
library(dplyr)
library(foreach)
library(doMC)
registerDoMC(cores=4)

#todo - set rank to the proper indices
#set colors to dark blue to dark red with white at AUC = 0.5

load('./processed_zeisel.Rdata', verbose=TRUE)
descriptions <- read_csv("celltype_descriptions.csv")
cores <- 1

cleaned_gene_list <- c("Gpr151")
forIndices <- linnarssonMatrix[,1, drop=F]
forIndices$Gene <- rownames(linnarssonMatrix)
forIndices %<>% mutate(isTargetGene = Gene %in% cleaned_gene_list)
targetIndices <- forIndices$isTargetGene

wilcoxTests <- foreach(oneCol=iter(linnarssonMatrix, by='col'), .combine=rbind) %dopar% {
  data.frame(auc = auroc_analytic(oneCol, as.numeric(targetIndices)), 
             pValue=wilcox.test(oneCol[targetIndices], oneCol[!targetIndices], conf.int = F)$p.value)
}
wilcoxTests$cluster_id <- colnames(linnarssonMatrix)

#add descriptions
wilcoxTests <- inner_join(descriptions, wilcoxTests)
wilcoxTests %<>% arrange(-auc)
wilcoxTests %<>% mutate(adjusted_P = p.adjust(pValue))

wilcoxTests %<>% sample_n(size = nrow(wilcoxTests))
wilcoxTests$rank <- 1:265 #to be based on the true ranking/indices
wilcoxTests$dummyY <- 1 

#plot with all the extras    
(rasterPlot <- ggplot(wilcoxTests, aes(x = rank, y = dummyY)) +
    geom_tile(aes(fill = auc)) +
    coord_cartesian(expand=F) 
    )

#plot AUC with no margins
(rasterPlot <- ggplot(pValues, aes(x = rank, y = dummyY)) +
    geom_tile(aes(fill = pValue)) )

gt = ggplotGrob(rasterPlot)

gt = gtable::gtable_filter(gt, "panel")  

# Draw it
grid.newpage()
grid.draw(gt)

