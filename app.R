#todo
#filter genome for those reachable from human genes when using human genes as input
#have a better dendrogram
#one-sided p for enrichment, not depletion
#switch reactivePlot to renderPlot

library(ggplot2)
library(grid)
library(png)
library(cowplot)
library(homologene)

library(readr)
library(reshape2)
library(doMC)
library(shiny)
library(dplyr)
library(magrittr)
library(tidyr)
library(shinyjs)

source("./string_processing.R")
source("./AUCFunction.R")

# load full dataset:
# --> linnarssonMatrixMouse // linnarssonMatrixHumanReachable
load('./processed_zeisel.Rdata', verbose = TRUE)
# load descriptions for full dataset
descriptions <- read_csv("celltype_descriptions.csv")
orders <- read_csv("png_tags.csv")

# load 39 celltypes dataset
# --> mouseMatrix_39_celltypes // humanReachableMatrix_39_celltypes
# also loads celltype descriptions --> cell_type_descriptions_39
load('processed_39_cellypes.Rdata', verbose = TRUE)

unique_genes_all <- rownames(linnarssonMatrixMouse)
unique_genes_human_reachable <-
  mouse2human(unique_genes_all)$mouseGene

cores <- 1

if( Sys.info()['nodename'] == "RES-C02RF0T2.local" ) { 
  cores <- 4
}


registerDoMC(cores=cores)


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
      selectInput(
        inputId = 'dataset',
        label = 'Dataset:',
        choices = c('Full dataset (265 celltypes)',
                    '39 celltypes')
      ),
      textAreaInput(
        inputId = "genelist",
        label = "Input your gene list:",
        value = 'Mag\nMobp\nMog\nMbp\nOmg',
        rows = 7
      ),
      selectInput('species', 'Species:',
                  choices = c('Mouse', 'Human')),
      actionButton("submit", "Submit"),
      br(),
      br(),
      downloadButton(outputId = "download_data", label = "Download results as .csv"),
      br(),
      br(),
      div(
        id = "main",
        p("This tool is made possible by data and the dendrogram from:"),
        a("Molecular architecture of the mouse nervous system", href = "http://dx.doi.org/10.1016/j.cell.2018.06.021"),
        p(
          " by Amit Zeisel, Hannah Hochgerner, Peter Lonnerberg, Anna Johnsson, Fatima Memic, Job van der Zwan, Martin Haring, Emelie Braun, Lars Borm, Gioele La Manno, Simone Codeluppi, Alessandro Furlan, Nathan Skene, Kenneth D Harris, Jens Hjerling Leffler, Ernest Arenas, Patrik Ernfors, Ulrika Marklund, and Sten Linnarsson."
        ),
        br(),
        p(
          "This tool was made by Derek Howard, Navona Calarco and Leon French during BrainHack 2018 Toronto."
        ),
        a("Source code", href = "https://github.com/derekhoward/polygenic_celltypes")
      )
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      # Output: Verbatim text for data summary ----
      verbatimTextOutput("summary"),
      br(),
      plotOutput("plot"),
      #,  width = "100%"),
      br(),
      # Output: HTML table with requested number of observations ----
      dataTableOutput("view")
      
    )
  )
)

server <- function(input, output) {
  output$summary <- renderPrint({
    shinyjs::disable("download_data")
    cat("\nResults will load here when complete")
    cat("\n")
    print(gc())
    print(Sys.info()['nodename'])
  })
  
  observeEvent(input$submit, {
    #eventReactive(input$submit, {
    shinyjs::hide("main")
    shinyjs::disable("submit")
    shinyjs::disable("download_data")
    start <- Sys.time()
    # should check to ensure that user inputs gene symbols of correct species
    # this doesn't work yet, try alternative to validate() that works with observeEvent()
    cleaned_gene_list <- isolate(if (length(process_input_genes(input$genelist)) == 0) {
      return()
    } else {
      process_input_genes(input$genelist)
    })
    
    if (input$species == 'Human') {
      cleaned_gene_list <- convert_genes(cleaned_gene_list)
      
      if (input$dataset == '39 celltypes') {
        linnarssonMatrix <- humanReachableMatrix_39_celltypes
        unique_genes <- rownames(linnarssonMatrix)
        descriptions <- cell_type_descriptions_39
        paste('Loaded 39 celltypes data')
      } else {
        # load standard human reachable dataset
        linnarssonMatrix <- linnarssonMatrixHumanReachable
        unique_genes <- rownames(linnarssonMatrix)
        descriptions <- descriptions
        paste('Loaded default 265 celltypes data')
      }
      
    } else {
      # species is mouse!
      if (input$dataset == '39 celltypes') {
        #load('processed_39_cellypes.Rdata', verbose = TRUE)
        linnarssonMatrix <- mouseMatrix_39_celltypes
        unique_genes <- rownames(linnarssonMatrix)
        descriptions <- cell_type_descriptions_39
        paste('Loaded 39 celltypes data')
      } else {
        # dataset is 265 celltypes for mouse
        linnarssonMatrix <- linnarssonMatrixMouse
        unique_genes <- unique_genes_all
      }
    }
    
    print(paste0("Before time taken:", Sys.time() - start))
    
    
    #for indices - use dplyr for ease
    forIndices <- linnarssonMatrix[, 1, drop = F]
    forIndices$Gene <- rownames(linnarssonMatrix)
    forIndices %<>% mutate(isTargetGene = Gene %in% cleaned_gene_list)
    targetIndices <- forIndices$isTargetGene
    
    wilcoxTests <-
      foreach(oneCol = iter(linnarssonMatrix, by = 'col'),
              .combine = rbind) %dopar% {
                data.frame(
                  AUROC = auroc_analytic(oneCol, as.numeric(targetIndices)),
                  pValue = wilcox.test(oneCol[targetIndices], oneCol[!targetIndices], conf.int = F)$p.value
                )
              }
    wilcoxTests$cluster_id <- colnames(linnarssonMatrix)
    
    
    print(paste0("Wilcox time taken:", Sys.time() - start))
    
    #add descriptions to table
    wilcoxTests <- inner_join(descriptions, wilcoxTests)
    wilcoxTests %<>% arrange(-AUROC)
    wilcoxTests %<>% mutate(rank=rank(-AUROC)) 
    print(wilcoxTests)
    
    output$summary <- renderPrint({
      #count of intersection of submitted genes with total gene list
      cat(paste("Time taken:", round(Sys.time() - start), "seconds"))
      cat(paste(
        "\nGenes found in data:",
        sum(cleaned_gene_list %in% unique_genes),
        "of",
        length(cleaned_gene_list)
      ))
      cat(paste(
        "\nUsing background set of",
        length(unique_genes),
        "genes"
      ))
    })
    
    output$view <- renderDataTable({
      if (input$dataset != '39 celltypes') {
        wilcoxTests %<>% mutate(
          cluster_id = sprintf(
            '<a href="http://mousebrain.org/celltypes/%s.html" target="_blank">%s</a>',
            cluster_id,
            cluster_id
          )
        )
      } else {
        wilcoxTests %<>% mutate(
          # modified this column such that there is only 'cluster_id' which links to the information linnarson lab website
          #cluster_info = sprintf(
          cluster_id = sprintf(
            '<a href="http://mousebrain.org/taxonomy/%s.html" target="_blank">%s</a>',
            cluster_info,
            cluster_id
          )
        )
        # drop cluster_info so you can reorder columns in same fashion regardless of using 39 celltypes or all of them
        wilcoxTests %<>% select(-cluster_info)
      }
      wilcoxTests %<>% mutate(
        pValue = signif(pValue, digits = 3),
        AUROC = signif(AUROC, digits = 3),
        adjusted_P = signif(p.adjust(pValue), digits = 3)
      )
      wilcoxTests %<>% select(rank, everything())#cluster_id, AUROC, pValue, adjusted_P)
      wilcoxTests
    }, escape = FALSE)
    
    #Figure
    if (input$dataset != '39 celltypes') {
      output$plot <- reactivePlot(function() {
        forPlot <- wilcoxTests %>% select(cluster_id, AUROC)
        forPlot <- inner_join(forPlot, orders)
        forPlot$dummyY <- 1
        
        #plot with all the extras
        
        (
          rasterPlot <- ggplot(forPlot, aes(x = index_in_png, y = dummyY)) +
            geom_tile(aes(fill = AUROC)) +
            coord_cartesian(expand = F) +
            scale_fill_gradientn(
              colours = c("darkblue", "white", "darkred"),
              values = c(0, .5, 1),
              space = "Lab",
              na.value = "grey50",
              guide = "colourbar",
              limits = c(0, 1)
            )
          
        )
        gt = ggplotGrob(rasterPlot)
        gt = gtable::gtable_filter(gt, "panel")
        
        img <- readPNG("./dendrogram-01.png")
        
        g <- rasterGrob(img, interpolate = TRUE)
        
        treeImage <- ggplot() +
          geom_blank() +
          annotation_custom(
            g,
            xmin = -Inf,
            xmax = Inf,
            ymin = -Inf,
            ymax = Inf
          )
        gt_tree <- ggplotGrob(treeImage)
        gt_tree <- gtable::gtable_filter(gt_tree, "panel")
        plot(gt_tree)
        
        final_p <-
          plot_grid(
            gt_tree,
            gt,
            ncol = 1,
            axis = "rltb",
            rel_heights = c(0.95, .05),
            align = "v"
          )
        return(final_p)
      }, height = 188, width = 922)
    }
    else {
      output$plot <- NULL
    }
    
    shinyjs::enable("submit")
    shinyjs::enable("download_data")
    
    output$download_data <-
      downloadHandler(
        filename = "polygenic_cell_types_AUC_results.csv",
        content = function(file) {
          write_csv(wilcoxTests, file)
        }
      )
    
  })
}

shinyApp(ui, server)
