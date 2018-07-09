#todo
#filter genome for those reachable from human genes when using human genes as input
#have a better dendrogram
#one-sided p for enrichment, not depletion
#switch reactivePlot to renderPlot

library(ggplot2)
library(grid)
library(png)
library(cowplot)

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
orders <- read_csv("png_tags.csv")

cores <- 1

if( Sys.info()['nodename'] == "RES-C02RF0T2.local" ) { 
  cores <- 4
}


registerDoMC(cores=cores)
print("Done loading data")

unique_genes <- rownames(linnarssonMatrix)

ui <- fluidPage(
  shinyjs::useShinyjs(),
  tags$head(includeHTML("google-analytics.html")),
  # App title ----
  titlePanel("Polygenic celltype tester for the mouse nervous system"),
  
  # Sidebar layout with a input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      # Input: Selector for choosing dataset ----
      textAreaInput(inputId = "genelist",
                    label = "Input your gene list:",
                    value = 'Mag\nMobp\nMog\nMbp\nOmg',
                    rows=7),
      selectInput('species', 'Species:',
                  choices=c('Mouse', 'Human')),
      actionButton("submit", "Submit")
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      div(id = "main",
      p("This tool is made possible by data and the dendrogram from:"),
      a("Molecular architecture of the mouse nervous system", href="https://doi.org/10.1101/294918"),
      p(" by Amit Zeisel, Hannah Hochgerner, Peter Lonnerberg, Anna Johnsson, Fatima Memic, Job van der Zwan, Martin Haring, Emelie Braun, Lars Borm, Gioele La Manno, Simone Codeluppi, Alessandro Furlan, Nathan Skene, Kenneth D Harris, Jens Hjerling Leffler, Ernest Arenas, Patrik Ernfors, Ulrika Marklund, and Sten Linnarsson."),
      br(),
      p("This tool was made by Derek Howard, Navona Calarco and Leon French during BrainHack 2018 Toronto."),
      a("Source code", href="https://github.com/derekhoward/polygenic_celltypes"),
      br(),
      br()),
      
      
      # Output: Verbatim text for data summary ----
      verbatimTextOutput("summary"),
      br(),
      plotOutput("plot"),#,  width = "100%"),
      br(),
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
  
  observeEvent(input$submit, {
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

    output$summary <- renderPrint({
      #count of intersection of submitted genes with total gene list
      cat(paste("Time taken:", round(Sys.time() - start), "seconds"))
      cat(paste("\nGenes found in data:",sum(cleaned_gene_list %in% unique_genes), "of", length(cleaned_gene_list)))
    })
    
    output$view <- renderDataTable({
      wilcoxTests %<>% mutate(cluster_id = sprintf('<a href="http://mousebrain.org/celltypes/%s.html" target="_blank">%s</a>', cluster_id, cluster_id))
      wilcoxTests %<>% mutate(pValue = signif(pValue, digits=3), auc = signif(auc, digits=3), adjusted_P = signif(p.adjust(pValue), digits=3))
      wilcoxTests
    }, escape = FALSE)
    
    #Figure
    output$plot <- reactivePlot(function() {
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

      img <- readPNG("./dendrogram-01.png")

      g <- rasterGrob(img, interpolate=TRUE)

      treeImage <- ggplot() +
        geom_blank() +
        annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
      gt_tree <- ggplotGrob(treeImage)
      gt_tree <- gtable::gtable_filter(gt_tree, "panel")
      plot(gt_tree)

      final_p <- plot_grid(gt_tree, gt, ncol=1, axis="rltb", rel_heights = c(0.95,.05), align="v")
      return(final_p)
    }, height=188, width = 922)
    
    shinyjs::enable("submit")
  }
  
  
  )
}

shinyApp(ui, server)
