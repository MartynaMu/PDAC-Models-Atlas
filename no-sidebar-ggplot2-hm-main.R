library(shiny)
library(bslib)
library(plotly)
library(PCAtools)
# UI============================================================================
ui <- fluidPage(
  titlePanel("PDAC In Vitro Models Atlas"),
  mainPanel(
    tabsetPanel(
      tabPanel("PCAplots",
               plotlyOutput("screeplot", height = 350, width = "50%"),
               selectInput("x_pc", "X axis PC", choices = colnames(p$rotated), selected = "PC1"),selectInput("y_pc", "Y axis PC", choices = colnames(p$rotated), selected = "PC2"),
               plotlyOutput("PCAplot", width = "50%")
      ),
      tabPanel(
        "PC and classification correlation",
        plotOutput("eigen", width = "50%")
      ),
      tabPanel(
        "Heatmap",
        fluidRow(
          column(3,
                 textInput("protein_ids", "Enter protein IDs (comma-separated):", value = example_ids)
          ),
          column(3,
                 fileInput("protein_file", "Or upload protein IDs file (.txt):", accept = c(".txt"))
          )
        ),
        actionButton("plot_btn", "Plot Heatmap"),
        fluidRow(
          column(6,
            plotlyOutput("hm_plot_prot")
          ),
          column(6,
            plotlyOutput("hm_plot_rna")
          )
        ),
        textOutput("not_found")
      )
    )
  )
)
# SERVER========================================================================
server <- function(input, output, session) {
  
  prot.means <- read.csv("../RR_proteomics/data/output/prot-means.csv") %>% column_to_rownames("Genes")
  rna.means <- read.csv("../RR_proteomics/data/output/rna-means.csv") %>% column_to_rownames("Genes")
  
  output$PCAplot <- renderPlotly({
    
    x <- input$x_pc
    y <- input$y_pc
    
    plot <- biplot(p,
                   lab = NULL,
                   colkey = metadata_colors$Model.detailed,
                   max.overlaps = FALSE,
                   colby = "Model.detailed",
                   shape = "Cell.line",
                   title = "Global proteome PCA, no filter",
                   hline = 0,
                   vline = 0,
                   legendPosition = "right",
                   x = x,
                   y = y)+ 
      aes(text=colnames(ds))
    ggplotly(plot, tooltip = "text") %>%
      layout(legend = list(traceorder="grouped"))
  })
  
  output$screeplot <- renderPlotly({
    screeplot <- screeplot(p, components = paste0("PC", 1:10), title = "Variance explained, PC1-10")
    ggplotly(screeplot, tooltip = "y")
  })
  
  output$eigen <- renderPlot({
    eigencorplot(p, 
                 metavars = colnames(metadata),
                 rotLabX = 45,
                 col = c('blue2', 'black', 'red2'),
                 colCorval = 'white',
                 main = "Global proteome PC correlation to sample classification",
                 posColKey = "top",
                 colFrame = "white",
                 cexMain = 1.7,
                 cexCorval = 1.1,
                 fontCorval = 2,
                 components = paste0("PC", c(1:10)))
  })
  
  selected_proteins <- eventReactive(input$plot_btn, {
    # Read protein IDs from file if uploaded
    file_ids <- NULL
    if (!is.null(input$protein_file)) {
      file_content <- readLines(input$protein_file$datapath)
      # Assume one protein ID per line, trim whitespace
      file_ids <- trimws(file_content)
    }
    
    # If file IDs exist and not empty, use them; else use text input
    ids <- if (!is.null(file_ids) && length(file_ids) > 0) {
      file_ids
    } else {
      trimws(unlist(strsplit(input$protein_ids, ",")))
    }
    
    found <- ids[ids %in% rownames(prot.means)]
    not_found <- ids[!ids %in% rownames(prot.means)]
    list(found = found, not_found = not_found)
  })
  
  
  output$not_found <- renderText({
    nf <- selected_proteins()$not_found
    if (length(nf) > 0) {
      paste("Not found:", paste(nf, collapse = ", "))
    }
  })
  
  ## HEATMAPS-------------------------------------------------------------------
  
  p_ann1 <- (ggplot(df.ann, aes(x=Group, y = .1, fill=Cell.line)) + 
               geom_tile(color="black")+
               scale_fill_brewer(palette = "Greys")) %>%
    ggplotly() %>% 
    layout(yaxis = list(showticklabels = FALSE,
                        ticks = ""))
  p_ann2 <- (ggplot(df.ann, aes(x=Group, y = .1, fill=Model.detailed)) + 
               geom_tile(color="black")+
               scale_fill_brewer(palette = "Paired",direction = -1)) %>%
    ggplotly() %>% 
    layout(yaxis = list(showticklabels = FALSE,
                        ticks = ""))
  ### PROT HM----------------------------------------------------------------------
  
  output$hm_plot_prot <- renderPlotly({
    req(length(selected_proteins()$found) > 0)
    mat <- prot.means[selected_proteins()$found, , drop = FALSE]
    
    n_proteins <- nrow(mat)
    
    mat.L <- mat %>% 
      scale() %>%
      as.data.frame() %>%
      rownames_to_column("Genes") %>%
      pivot_longer(2:16,
                   names_to = "Group",
                   values_to = "Mean.Intensity")
    
    p_hm <- (ggplot(mat.L, aes(x=Group, y= Genes, fill=Mean.Intensity)) + 
      geom_tile()+
      scale_fill_distiller(palette = "RdBu")) %>%
      ggplotly()
    
    ann_height = 25
    hm_height = n_proteins * 15
    ann_prop <- ann_height / (ann_height + hm_height)
    hm_prop <- 1 - 2*ann_prop
    plot_height <- 2*ann_height + hm_height
    
    subplot(p_ann1,
            p_ann2,
            p_hm, 
            nrows=3, 
            shareX = TRUE,
            heights = c(ann_prop, ann_prop, hm_prop), 
            margin = 0) %>%
      layout(height = plot_height)
    
  })
  ### RNA HM--------------------------------------------------------------------
  output$hm_plot_rna <- renderPlotly({
    req(length(selected_proteins()$found) > 0)
    mat <- rna.means[selected_proteins()$found, , drop = FALSE]
    
    n_proteins <- nrow(mat)
    
    mat.L <- mat %>% 
      scale() %>%
      as.data.frame() %>%
      rownames_to_column("Genes") %>%
      pivot_longer(2:16,
                   names_to = "Group",
                   values_to = "Mean.Intensity")
    
    p_hm <- (ggplot(mat.L, aes(x=Group, y= Genes, fill=Mean.Intensity)) + 
               geom_tile()+
               scale_fill_distiller(palette = "RdBu")) %>%
      ggplotly()
    
    ann_height = 25
    hm_height = n_proteins * 15
    ann_prop <- ann_height / (ann_height + hm_height)
    hm_prop <- 1 - 2*ann_prop
    plot_height <- 2*ann_height + hm_height
    
    subplot(p_ann1,
            p_ann2,
            p_hm, 
            nrows=3, 
            shareX = TRUE,
            heights = c(ann_prop, ann_prop, hm_prop), 
            margin = 0) %>%
      layout(height = plot_height)
    
  })
  
}
# Run===========================================================================
shinyApp(ui = ui, server = server)

