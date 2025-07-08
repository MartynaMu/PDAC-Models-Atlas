library(shiny)
library(bslib)
library(plotly)
library(PCAtools)
# UI============================================================================
ui <- fluidPage(
  titlePanel("PDAC In Vitro Models Atlas"),
  sidebarLayout(
    sidebarPanel(plotlyOutput("screeplot", height = 350, width = "100%"),
                 selectInput("x_pc", "X axis PC", choices = colnames(p$rotated), selected = "PC1"),
                 selectInput("y_pc", "Y axis PC", choices = colnames(p$rotated), selected = "PC2")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("PCAplots",
                 splitLayout(cellWidths = c("50%", "50%"),
                             plotlyOutput("PCAplot")
                 )
        ),
        tabPanel(
          "PC and classification correlation",
          plotOutput("eigen", width = "50%")
        )
      )
    )
  )
)
# SERVER========================================================================
server <- function(input, output) {
  
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
}
# Run===========================================================================
shinyApp(ui = ui, server = server)

