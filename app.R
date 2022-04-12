library(shiny)
library(ggplot2)
library(colourpicker)
library(DT)


ui <- fluidPage(
  titlePanel("BF591 Assignment 7"),
  p("To use this application, download the CSV deseq_res.csv from the data directory of this app's repository."),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file",
                label = "Load differential expression results",
                placeholder = "deseq_res.csv",
                multiple = FALSE,
                accept = ".csv"),
      p('A volcano plot can be generated with "log2 fold-change" on the x-axis and "p-adjusted" on the y-axis.'),
      radioButtons("x_name", "Choose the column for the x-axis",
                   choices = c(baseMean = "baseMean",
                               log2FoldChange = "log2FoldChange",
                               lfcSE = "lfcSE",
                               stat = "stat",
                               pvalue = "pvalue",
                               padj = "padj"),
                   selected = "log2FoldChange"),
      radioButtons("y_name", "Choose the column for the y-axis",
                   choices = c(baseMean = "baseMean",
                               log2FoldChange = "log2FoldChange",
                               lfcSE = "lfcSE",
                               stat = "stat",
                               pvalue = "pvalue",
                               padj = "padj"),
                   selected = "padj"),
      colourInput("base", "Base point color", value = "midnightblue"),
      colourInput("highlight", "Highlight point color", value = "gold"),
      sliderInput("slider", "Select the magnitude of the p adjusted coloring:",
                  min = -300, max = 0, value = -150),
      submitButton(text = "Plot", icon = icon("dog", lib = "font-awesome"), width = '100%')
    ),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Plot", plotOutput("volcano")),
                  tabPanel("Table", DT::dataTableOutput("table"))
      )
    )
  )
)

server <- function(input, output, session) {
  
  # Load data
  load_data <- reactive({
    validate(
      need(input$file, "\nPlease select a file"),
    )
    datafile <- input$file
    return(read.csv(datafile$datapath))
  })
  
  # Create the volcano plot
  volcano_plot <-
    function(dataf, x_name, y_name, slider, highlight, base) {
      ggplot(dataf) +
        geom_point(aes(x=get(x_name), y=-log10(get(y_name)), color = padj < 10**slider)) +
        scale_color_manual(name = paste('padj < 1e', as.character(slider), sep =""), values = setNames(c(highlight, base), c(T, F))) +
        theme_bw() +
        theme(legend.position="bottom") +
        xlab(as.character(x_name)) +
        ylab(as.character(y_name))
    }
  
  # Draw and format the table
  draw_table <- function(dataf, slider) {
    df <- dataf[which(dataf$padj < 10**slider),]
    df$baseMean <- round(df$baseMean, 1)
    df$log2FoldChange <- round(df$log2FoldChange, 2)
    df$lfcSE <- round(df$lfcSE, 2)
    df$stat <- round(df$stat, 2) 
    df$pvalue <- formatC(df$pvalue, digits = 3)
    df$padj <- formatC(df$padj, digits = 3)
    df
  }
  
  # Return volcano plot output
  output$volcano <- renderPlot({volcano_plot(load_data(), input$x_name, input$y_name, input$slider, input$highlight, input$base)})
  
  # Return output table
  output$table <- DT::renderDataTable({draw_table(load_data(), input$slider)}, rownames=FALSE)
}

# Run the application
shinyApp(ui = ui, server = server)