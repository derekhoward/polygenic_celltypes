library(readr)
library(reshape2)
library(doMC)
library(shiny)
library(dplyr)
library(magrittr)
library(shinyjs)

source("./string_processing.R")
source("./AUCFunction.R")
load('./processed_zeisel.Rdata', verbose=TRUE)
descriptions <- read_csv("celltype_descriptions.csv")
cores <- 1

if( Sys.info()['nodename'] == "RES-C02RF0T2.local" ) { 
  cores <- 4
}


registerDoMC(cores=cores)
print("Done loading data")

unique_genes <- rownames(linnarssonMatrix)

ui <- fluidPage(
  shinyjs::useShinyjs(),  
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
      selectInput('species', 'Species:',
                  choices=c('Mouse', 'Human')),
      actionButton("submit", "Submit")
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      div(id = "main",
      p("This tool is made possible by data from:"),
      a("Molecular architecture of the mouse nervous system", href="https://doi.org/10.1101/294918"),
      p(" by Amit Zeisel, Hannah Hochgerner, Peter Lonnerberg, Anna Johnsson, Fatima Memic, Job van der Zwan, Martin Haring, Emelie Braun, Lars Borm, Gioele La Manno, Simone Codeluppi, Alessandro Furlan, Nathan Skene, Kenneth D Harris, Jens Hjerling Leffler, Ernest Arenas, Patrik Ernfors, Ulrika Marklund, and Sten Linnarsson."),
      br(),
      br()),
      
      
      # Output: Verbatim text for data summary ----
      verbatimTextOutput("summary"),
      
      # Output: HTML table with requested number of observations ----
      dataTableOutput("view")
      
    )
  )
)

# Define server logic process and output top celltypes ----
server <- function(input, output) {
  output$summary <- renderPrint({
    cat(paste("cores set to", cores))
    cat("\nResults will load here when complete")
    cat("\n")
    print(gc())
    print(Sys.info()['nodename'])
  })
  
  observe({
    if (input$submit > 0) {
      shinyjs::hide("main")
      shinyjs::disable("submit") 
      start <- Sys.time()
      cleaned_gene_list <- isolate(process_input_genes(input$genelist))
      if (input$species == 'Human') {
        cleaned_gene_list <- convert_genes(cleaned_gene_list)
      }
      
      print(paste0("Before time taken:", Sys.time() - start))

      #for indices - use dplyr for ease
      forIndices <- linnarssonMatrix[,1, drop=F]
      forIndices$Gene <- rownames(linnarssonMatrix)
      forIndices %<>% mutate(isTargetGene = Gene %in% cleaned_gene_list)
      targetIndices <- forIndices$isTargetGene
      
      

      wilcoxTests <- foreach(oneCol=iter(linnarssonMatrix, by='col'), .combine=rbind) %dopar% {
        data.frame(auc = auroc_analytic(oneCol, as.numeric(targetIndices)), 
                   pValue=wilcox.test(oneCol[targetIndices], oneCol[!targetIndices], conf.int = F)$p.value)
      }
      wilcoxTests$cluster_id <- colnames(linnarssonMatrix)

      print(paste0("Wilcox time taken:", Sys.time() - start))

      #add descriptions
      wilcoxTests <- inner_join(descriptions, wilcoxTests)
      wilcoxTests %<>% arrange(-auc)
      wilcoxTests %<>% mutate(adjusted_P = p.adjust(pValue))
      
      output$summary <- renderPrint({
        #count of intersection of submitted genes with total gene list
        cat(paste("Time taken:", round(Sys.time() - start), "seconds"))
        cat(paste("\nGenes found in data:",sum(cleaned_gene_list %in% unique_genes), "of", length(cleaned_gene_list)))
      })
      
      output$view <- renderDataTable({
        wilcoxTests %<>% mutate(cluster_id = sprintf('<a href="http://mousebrain.org/doku.php?id=clusters:%s" target="_blank">%s</a>', cluster_id, cluster_id))
        wilcoxTests %<>% mutate(pValue = signif(pValue, digits=3), auc = signif(auc, digits=3))
        wilcoxTests
      }, escape = FALSE)
    }
    shinyjs::enable("submit")
    
  })
}

shinyApp(ui, server)
