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
          plotlyOutput("heatmap_plot", width = "50%"),
          textOutput("not_found")
        )
      )
  )
)
# SERVER========================================================================
server <- function(input, output, session) {
  
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
  
  output$heatmap_plot <- renderPlotly({
    req(length(selected_proteins()$found) > 0)
    mat <- prot.means[selected_proteins()$found, , drop = FALSE]
    n_proteins <- nrow(mat)
    
    plot_height = 300 + (n_proteins * 20)
    
    heatmaply::heatmaply(mat, main = "Protein Heatmap",
                         colors = colorRampPalette(c("blue", "white", "red"))(256),
                         col_side_colors = ann, 
                         Rowv = FALSE, 
                         Colv = FALSE,
                         plot_method = "plotly",
                         scale="row",
                         colorbar_len = .5,
                         side_color_colorbar_len = .5,
                         height = plot_height)
  })

}
# Run===========================================================================
shinyApp(ui = ui, server = server)

