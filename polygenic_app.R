library(shiny)
library(dplyr)
library(magrittr)
source("./string_processing.R")
source("./AUCFunction.R")
load('./processed_zeisel.Rdata')
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
      tableOutput("view")
      
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
      
      #do AUROC with gene list
      wilcoxTests <- linnarsson %>% group_by(cluster_id) %>% summarize(
        pValue = wilcox.test(log1ExpZ ~ isTargetGene)$p.value, 
        auc = auroc_analytic(rank(log1ExpZ), as.numeric(isTargetGene)))
      
      wilcoxTests %<>% arrange(-auc)
      
      output$summary <- renderPrint({
        #count of intersection of submitted genes with total gene list
        print(sum(cleaned_gene_list %in% unique_genes))
      })
      
      # Show the first top observations ----
      #output$view <- renderTable({
      #  head(wilcoxTests)
      output$view <- DT::dataTableOutput({
       head(wilcoxTests)
      })
      }
  })
}

shinyApp(ui, server)
