library(stringr)
library(homologene)

process_input_genes <- function(input_genes) {
  processed_genes <- unlist(str_split(trimws(input_genes), "[, \n]+") )
  return(processed_genes)
}

convert_genes <- function(input_genes) {
  mouse_genes <- human2mouse(input_genes)
  return(unique(mouse_genes$mouseGene))
}


genes <- c('CADM2', 'ZNF704', 'NCAM1', 'RABEP2', 'ATP2A1')
