library(shiny)
library(tidyverse)
library(bslib)
library(plotly)
library(PCAtools)
# UI============================================================================
ui <- fluidPage(
  titlePanel("PDAC Models Atlas"),
  mainPanel(
    tabsetPanel(
      tabPanel("PCAplots",
               plotlyOutput("screeplot", height = 350, width = "50%"),
               plotOutput("eigen", width = "50%"),
               selectInput("x_pc", "X axis PC", choices = colnames(p$rotated), selected = "PC1"),selectInput("y_pc", "Y axis PC", choices = colnames(p$rotated), selected = "PC2"),
               plotlyOutput("PCAplot", width = "50%"),
               plotlyOutput("loadplot", width = "50%")
      ),
      tabPanel("Explore",
               fluidRow(
                 column(3,
                        textInput("protein_ids", "Enter protein IDs (comma-separated):", value = example_ids)
                 ),
                 column(3,
                        fileInput("protein_file", "Or upload protein IDs file (.txt):", accept = c(".txt"))
                 )
               ),
               actionButton("plot_btn", "Plot Heatmap"),
               tags$div(
                 style = "width: 1200px; overflow-x: scroll;",
                 plotlyOutput("dotplot", height = "600px")
               ),
               textOutput("not_found")
      ),
      tabPanel(
        "Heatmap",
        plotlyOutput("hm_plot"),

      )
    )
  )
)
# SERVER========================================================================
server <- function(input, output, session) {
  
  prot.means <- read.csv("../RR_proteomics/data/output/prot-means.csv") %>% column_to_rownames("Genes")
  rna.means <- read.csv("../RR_proteomics/data/output/rna-means.csv") %>% column_to_rownames("Genes")
  colnames(rna.means) <- colnames(rna.means) %>% str_sub(start = 5)
  
  
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
  
  output$loadplot <- renderPlotly({
    
    x <- input$x_pc |> paste0("$")
    y <- input$y_pc |> paste0("$")
    
    (p$loadings %>% 
        select(matches(c(x,y))) %>% 
        rownames_to_column("Genes") %>%
        pivot_longer(2:3, 
                     names_to = "PC",
                     values_to = "Component.loading") %>%
        filter(abs(Component.loading) > 0.01) %>%
        ggplot(aes(x=Component.loading,
                   y=PC,
                   fill=Component.loading,
                   text=Genes))+
        geom_point(size=7)+
        geom_hline(yintercept = 0, linetype=2)+
        geom_vline(xintercept = 0, linetype=2)+
        scale_fill_distiller(palette = "RdBu")+
        labs(title=paste0(input$x_pc, ", ", input$y_pc, " loadings,\nabsolute value > 0.01"))
    ) %>% 
      ggplotly(tooltip = c("x", "text"))
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
  
  output$hm_plot <- renderPlotly({
    
    mat.L <- mat %>% 
      scale() %>%
      as.data.frame() %>%
      rownames_to_column("Genes") %>%
      pivot_longer(2:16,
                   names_to = "Group",
                   values_to = "Mean.Intensity")
    
    p_hm_prot <- (ggplot(mat.L, aes(x=Group, y= Genes, fill=Mean.Intensity)) + 
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
            
            nrows=3, 
            shareX = TRUE,
            heights = c(ann_prop, ann_prop, hm_prop),
            margin = 0) %>%
      layout(height = plot_height,
             xaxis = list(title = "",
                          tickangle = 45),
             xaxis2 = list(title = "",
                           tickangle = 45),
             yaxis = list(title = ""),
             yaxis2 = list(title = ""))
    
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
  
  output$dotplot <- renderPlotly({
    
    req(length(selected_proteins()$found) > 0)
    mat <- int.cor.L %>% 
      filter(Genes %in% selected_proteins()$found)
    
    n_proteins <- mat %>%
      select(Genes) %>%
      n_distinct()
    
    hm_width = n_proteins * 20
    plot_width = hm_width + 150
    limit <- mat$Intensity.z.score %>%
      range() %>%
      abs() %>%
      max()
    
    (mat %>%
       ggplot(aes(x=Genes,
                  y=Group,
                  fill=Intensity.z.score,
                  size=abs(Cell.line.corr),
                  text=paste0(
                    Genes,"<br>",
                    Group, "<br>",
                    "Log2 Intensity: ", round(Intensity,2), "<br>",
                    "Cell-line protein-mRNA \nPearson correlation: ", round(Cell.line.corr,2)
                  )
                  )
              )+
        geom_point()+
        geom_hline(yintercept = c(5.5,10.5), linetype=2)+
        scale_fill_distiller(palette = "RdBu",
                             limit=c(-limit,limit))+
        labs(x=NULL, 
             y=NULL,
             fill = "Intensity \nz-score")+
        theme(axis.text.x = element_text(angle=45))+
        annotate("rect", 
                 xmin=0,
                 xmax=-1, 
                 ymin=c(0,5.5,10.5), 
                 ymax=c(5.5,10.5,16), 
                 fill = RColorBrewer::brewer.pal("Greys", n=3), 
                 color="black")) |>
      ggplotly(width = plot_width,
               tooltip="text") |>
      layout(legend=list(
        x = 0, y = 1, xanchor = "left", yanchor = "top"))
    
  })
  
}
# Run===========================================================================
shinyApp(ui = ui, server = server)

