#devtools::install_github("hadley/multidplyr")
library(shiny)
library(dplyr)
library(magrittr)
library(multidplyr)
source("./string_processing.R")
source("./AUCFunction.R")
load('./processed_zeisel.Rdata')
linnarsson %<>% group_by(cluster_id)
print("Done loading data")

unique_genes <- unique(linnarsson$Gene)
# Define UI for dataset viewer app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Polygenic celltypes in mouse"),
  
  # Sidebar layout with a input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      # Input: Selector for choosing dataset ----
      textAreaInput(inputId = "genelist",
                    label = "Input your gene list:",
                    value = 'Calca',
                    rows=5),
      actionButton("submit", "Submit")
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Verbatim text for data summary ----
      verbatimTextOutput("summary"),
      
      # Output: HTML table with requested number of observations ----
      dataTableOutput("view")
      
    )
  )
)

# Define server logic process and output top celltypes ----
server <- function(input, output) {
  
  observe({
    if (input$submit > 0) {
      cleaned_gene_list <- isolate(process_input_genes(input$genelist))
      #find the celltypes which express input genes
      linnarsson %<>% mutate(isTargetGene = Gene %in% cleaned_gene_list)
      
      start <- Sys.time()
      #do AUROC with gene list
      wilcoxTests <- linnarsson %>% partition(cluster_id) %>% summarize(
        pValue = wilcox.test(log1ExpZ ~ isTargetGene, correct=F, conf.int=F)$p.value, 
        ) %>% collect()
      print(paste0("Wilcox time taken:", Sys.time() - start))
      aucs <- linnarsson %>% summarize(
        auc = auroc_analytic(rank(log1ExpZ), as.numeric(isTargetGene))
      ) 
      print(paste0("AUCs time taken:", Sys.time() - start))
      wilcoxTests <- inner_join(wilcoxTests, aucs)

      wilcoxTests %<>% arrange(-auc)
      
      output$summary <- renderPrint({
        #count of intersection of submitted genes with total gene list
        cat(paste("Time taken:", round(Sys.time() - start), "seconds"))
        cat(paste("\nGenes found in data:",sum(cleaned_gene_list %in% unique_genes), " of ", length(cleaned_gene_list)))
      })
      
      output$view <- renderDataTable({
        wilcoxTests %<>% mutate(cluster_id = sprintf('<a href="http://mousebrain.org/doku.php?id=clusters:%s" target="_blank">%s</a>', cluster_id, cluster_id))
        wilcoxTests
      }, escape = FALSE)
      }
  })
}

shinyApp(ui, server)
