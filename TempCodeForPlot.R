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

load('./processed_zeisel.Rdata', verbose=TRUE)
descriptions <- read_csv("celltype_descriptions.csv")
orders <- read_csv("png_tags.csv")
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

forPlot <- wilcoxTests %>% select(cluster_id, auc)
forPlot <- inner_join(forPlot, orders)
forPlot$dummyY <- 1 

#plot with all the extras    

(rasterPlot <- ggplot(forPlot, aes(x = index_in_png, y = dummyY)) +
    geom_tile(aes(fill = auc)) +
    coord_cartesian(expand=F) +
    scale_fill_gradientn(colours = c("darkblue", "white","darkred"), values = c(0, .5, 1), space = "Lab",
                           na.value = "grey50", guide = "colourbar", limits=c(0,1))
  
    )


gt = ggplotGrob(rasterPlot)

gt = gtable::gtable_filter(gt, "panel")  

# Draw it
grid.newpage()
grid.draw(gt)
plot(gt)

library(png)
library(cowplot)
img <- readPNG("./dendrogram-01.png")

g <- rasterGrob(img, interpolate=TRUE) 

treeImage <- ggplot() + 
  geom_blank() + 
  annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
gt_tree <- ggplotGrob(treeImage)
gt_tree <- gtable::gtable_filter(gt_tree, "panel")  
plot(gt_tree)

plot_grid(gt_tree, gt, ncol=1, axis="rlt", rel_heights = c(0.95,.05), align="vh")
