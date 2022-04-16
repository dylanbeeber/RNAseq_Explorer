library(shiny)
library(tidyverse)
library(colourpicker)
library(DT)
library(bigPint)

ui <- fluidPage(
  navbarPage("RNAseq Explorer",
    tabPanel("Sample Explorer",
      sidebarPanel(
        p("Select a data set:"),
        selectInput("dataset", "Data Set",
                    c("Soybean" = "se_soybean_cn_sub")),
        submitButton(text = "Plot", icon = icon("dog", lib = "font-awesome"), width = '100%')
      ),
      mainPanel(
        p("hello"),
        DT::dataTableOutput(data.frame(c(1, 2, 3), c(2, 4, 3)))
      )
    ),
    tabPanel("Counts Explorer",
    ),
    tabPanel("Differential Expression",
    ),
    tabPanel("Other",
    )
  )
)

server <- function(input, output) {
  
}

# Run the application
shinyApp(ui = ui, server = server)