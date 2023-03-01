library(shiny)
library(tidyverse)
library(colourpicker)
library(DT)
library(bigPint)
library(shinyjs)
library(SummarizedExperiment)
library(ggbeeswarm)
library(gplots)

ui <- fluidPage(
  shinyjs::useShinyjs(),
  tags$head(tags$style(HTML(".shiny-split-layout > div {overflow: visible}"))),
  navbarPage("RNAseq Explorer",
    tabPanel("Sample Explorer",
      sidebarPanel(
        p("Select a data set:"),
        selectInput("dataset", "Data Set",
                    c("Choose one" = "",
                      "Soybean" = "se_soybean_cn_sub")),
        submitButton(text = "Select", width = '100%'),
        p(),
        downloadButton("DL_counts", label = "Download Counts Data", width = '100%'),
        p(),
        downloadButton("DL_DE", label = "Download DE Data", width = '100%'),
      width = 3),
      mainPanel(
        DT::dataTableOutput("metadata_table")
      )
    ),
    tabPanel("Counts Explorer",
      sidebarPanel(
        # Remove file input for now, as data source is determined by first tab
        #fileInput("counts_data",
          #label = "Load Counts Matrix",
          #multiple = FALSE,
          #accept = ".csv"),
        p("Filter gene based on percentile of variance:"),
        sliderInput("counts_filter", "Counts Filter",
                    min = 0, max = 1, value = .5),
        p("Filter genes with at most N samples of nontrivial expression:"),
        sliderInput("num_samples", "Number of samples with expression greater than 1",
                    min = 0, max = 9, value = 1),
        submitButton(text = "Submit", width = '100%'),
        width = 3),
      mainPanel(
        tabsetPanel(type = "tabs",
                    tabPanel("Summary Table", tableOutput("filtered_table"),
                             p("Table displaying the results of filtering based on
                               the sidebar sliders")),
                    tabPanel("Scatter Plot", plotOutput("mean_var_scatter"),
                             p("Plot displaying the ranked mean vs. the variance of
                               filtered expression data")),
                    tabPanel("Heatmap", p(), plotOutput("heatmap"),
                             p("Clustered heatmap displaying the groupings of samples
                               based on filtered gene expression")),
                    tabPanel("PCA Plot",
                             fluidRow(
                               splitLayout(
                                  selectInput("x_axis", "X Axis",
                                             c("PC1" = "PC1",
                                               "PC2" = "PC2",
                                               "PC3" = "PC3",
                                               "PC4" = "PC4",
                                               "PC5" = "PC5",
                                               "PC6" = "PC6",
                                               "PC7" = "PC7",
                                               "PC8" = "PC8",
                                               "PC9" = "PC9"),
                                             selected = "PC1"),
                                  selectInput("y_axis", "Y Axis",
                                             c("PC1" = "PC1",
                                               "PC2" = "PC2",
                                               "PC3" = "PC3",
                                               "PC4" = "PC4",
                                               "PC5" = "PC5",
                                               "PC6" = "PC6",
                                               "PC7" = "PC7",
                                               "PC8" = "PC8",
                                               "PC9" = "PC9"),
                                             selected = "PC2"),
                               )
                             ),
                             plotOutput("pca_plot"),
                             p("Principal component analysis of filtered gene expression
                               data"))
        )
      )
    ),
    tabPanel("Differential Expression",
      sidebarPanel(
        # Remove file input for now, as data source is determined by first tab
        #fileInput("de_data",
          #label = "Load Differential Expression Matrix",
          #multiple = FALSE,
          #accept = ".csv"),
        p("First select a comparison:"),
        selectInput("comp_select", "Comparison",
                    c("Early vs. Middle" = "S1_S2",
                      "Middle vs. Late" = "S2_S3",
                      "Early vs. Late" = "S1_S3")),
        radioButtons("x_name", "Choose the column for the x-axis",
                     choices = c(log2FoldChange = ".logFC",
                                 logCPM = ".logCPM",
                                 likelihoodratio = ".LR",
                                 pvalue = ".PValue",
                                 falsediscoveryrate = ".FDR"),
                     selected = ".logFC"),
        radioButtons("y_name", "Choose the column for the y-axis",
                     choices = c(log2FoldChange = ".logFC",
                                 logCPM = ".logCPM",
                                 likelihoodratio = ".LR",
                                 pvalue = ".PValue",
                                 falsediscoveryrate = ".FDR"),
                     selected = ".PValue"),
        colourpicker::colourInput("base", "Base point color", value = "midnightblue"),
        colourpicker::colourInput("highlight", "Highlight point color", value = "gold"),
        sliderInput("slider", "Select the magnitude of the p adjusted coloring:",
                    min = -10, max = 0, value = -5),
        submitButton(text = "Submit", width = '100%'),
        width = 3),
      mainPanel(
        tabsetPanel(type = "tabs",
                    tabPanel("DE Table", DT::dataTableOutput("filtered_de_data"),
                             p("Table displaying the filtered results of differential
                               expression analysis")),
                    tabPanel("Volcano Plot", plotOutput("volcano_plot"),
                             p("Plot displaying the results of differential expression
                               analysis")),
        )
      )
    ),
    tabPanel("Single Gene Visualization",
      sidebarPanel(
        # Remove file input for now, as data source is determined by first tab
        #fileInput("sg_counts_data",
                  #label = "Load Counts Matrix",
                  #multiple = FALSE,
                  #accept = ".csv"),
        uiOutput("moreControls"),
        selectInput("plot_type", "Plot Type",
                    c("Bar Plot" = "bar_plot",
                      "Boxplot" = "boxplot",
                      "Violin Plot" = "violin_plot",
                      "Beeswarm Plot" = "beeswarm_plot")),
        submitButton(text = "Submit", width = '100%'),
        width = 3),
      mainPanel(plotOutput("single_gene_plot")),
    )
  )
)

server <- function(input, output) {
  
  # Code for sample explorer
  
  # Create function to take summarized experiment as text and output
  # a metadata tibble
  sample_info <- function(data) {
    # Create dictionary for later use in the function
    encoding_dict <- c("1"="early", "2"="middle", "3"="late")
    
    # Create metadata tibble
    metadata <- tibble(colnames(get(data)))
    # Rename samples column
    metadata <- dplyr::rename(metadata, samples = colnames(metadata)[1]) %>%
      # Add column for replicates
      mutate(replicate = str_extract(samples, "(?<=\\.).")) %>%
      # Add column for the growth stage
      mutate(stage = encoding_dict[str_extract(samples, "(?<=S).")])
    # Convert stage column to factor
    metadata$stage <- as.factor(metadata$stage)
    metadata
  }
  
  # Metadata table output
  output$metadata_table <- DT::renderDataTable({load_sample_data()}, rownames=FALSE)
  
  # Download counts data
  output$DL_counts <- downloadHandler(
    # Set name of file
    filename = function() {
      paste(input$dataset, '_counts.csv', sep='')
    },
    content = function(con) {
      # Extract counts data and write to .csv
      write.csv(assay(get(input$dataset)), con)
    }
  )
  
  # Download DE data
  output$DL_DE <- downloadHandler(
    # Set name of file
    filename = function() {
      paste(input$dataset, '_de.csv', sep='')
    },
    content = function(con) {
      # Extract DE data and write to .csv
      write.csv(rowData(get(input$dataset)), con)
    }
  )
  
  
  
  # Code for counts explorer
  
  # Load sample data
  load_sample_data <- reactive({
    # Ensure dataset is selected before loading table
    validate(
      need(input$dataset, "\nPlease select a dataset"),
    )
    dataset <- input$dataset
    # Format sample data for display
    return(sample_info(dataset))
  })
  
  # Load counts data
  load_counts_data <- reactive({
    validate(
      need(input$dataset, "\nPlease select dataset in first tab"),
    )
    return(read.csv(paste("data/", input$dataset, "_counts.csv",  sep = "", collapse = "")))
  })
  
  # Function to remove genes with certain number of rows less than one
  custom_filter <- function(data, percentile, n_samples){
    # Remove gene list from data
    selected_data <- dplyr::select(data, -X)
    # Get variance counts
    variance_counts <- apply(selected_data, var, MARGIN = 1)
    # Get a value for the specified percentile
    percentile_val <- quantile(variance_counts, percentile)
    # Create boolean list for genes over that percentile
    check <- variance_counts >= percentile_val
    # Filter down list of data based on percentiles
    filtered <- dplyr::filter(data, check == TRUE)
    # Again remove gene list from filtered data
    selected_data <- dplyr::select(filtered, -X)
    # Create function for checking if rows are less than one
    # Note I use one here instead of zero because of the way normalization
    # is performed
    parameter_check <- function(x){
      sum((x) > 1)
    }
    # Create boolean list for genes that meet the criteria
    check <- apply(selected_data, parameter_check, MARGIN = 1) >= n_samples
    # Produce final list based on booleans
    dplyr::filter(filtered, check == TRUE)
  }
  
  # Make table of filter related outputs
  filter_table_results <- function(data, percentile, n_samples) {
    num_samples <- ncol(data)-1
    num_genes <- nrow(data)
    # Filter data with custom filter
    filtered_data <- custom_filter(data, percentile, n_samples)
    num_passing_genes <- nrow(filtered_data)
    num_excluded_genes <- num_genes - num_passing_genes
    # Create data frame
    first_col <- c("Number of samples", "Number of genes", "Number of passing genes", "Number of excluded genes")
    second_col <- c(num_samples, num_genes, 
                    num_passing_genes, num_excluded_genes)
    second_col <- round(second_col, 0)
    second_col[3] <- paste(second_col[3], " (", round(num_passing_genes/num_genes * 100, 1), "%)", sep = "")
    second_col[4] <- paste(second_col[4], " (", round(num_excluded_genes/num_genes * 100, 1), "%)", sep = "")
    results <- data.frame(first_col, second_col)
    colnames(results) <- c("Attributes", "Results")
    results
  }
  
  # Create a table of filtered summary results
  output$filtered_table <- renderTable({validate(
                                        need(input$dataset, "\nPlease select a dataset in first tab"))
                                        filter_table_results(load_counts_data(),
                                                             input$counts_filter,
                                                             input$num_samples)})
  
  # Create scatter plot
  create_mean_var_scatter <- function(data, percentile, n_samples) {
    # Filter data
    filtered_data <- custom_filter(data, percentile, n_samples)
    # Remove gene column
    selected_data <- dplyr::select(data, -X)
    # Calculate mean counts
    mean_counts <- apply(selected_data, mean, MARGIN = 1)
    # Calculate variance
    variance_counts <- apply(selected_data, var, MARGIN = 1)
    # Rank mean values
    ranked_means <- rank(mean_counts)
    # Check which data has been filtered
    check <- data$X %in% filtered_data$X
    # Create tibble for the variance/mean data
    temp <- tibble(mean_counts, variance_counts, ranked_means, check)
    # Plot
    p <- ggplot(temp, aes(x = ranked_means, y = variance_counts, col = check)) + 
      geom_point() +
      scale_color_manual(name = "Filtered", values = c("gray", "black"), c(T, F)) +
      geom_smooth(method='gam', formula = y ~ s(x, bs = "cs"), color = "red") +
      labs(title = "Rank vs. Variance",
           x = "Rank(Mean)",
           y = "Variance") +
      scale_y_continuous(trans = 'log10')
    p
  }
  
  # Create a scatter plot
  output$mean_var_scatter <- renderPlot({create_mean_var_scatter(load_counts_data(),
                                                               input$counts_filter,
                                                               input$num_samples)})
  
  # Function for creating a heatmap
  create_heatmap <- function(data, percentile, n_samples) {
    # Filter data
    filtered <- custom_filter(data, percentile, n_samples)
    # Convert filtered data to matrix
    matrix_data <- as.matrix(dplyr::select(filtered, -X))
    # Add gene names to matrix
    rownames(matrix_data) <- filtered$X
    # Create heatmap
    heatmap.2(matrix_data, trace = "none", ylab="Genes",
              ColSideColors=rep(c("green","blue","purple"), each=3),
              xlab="Samples", labRow = FALSE, key = TRUE)
  }
  
  # Create the heatmap
  output$heatmap <- renderPlot({create_heatmap(load_counts_data(),
                                               input$counts_filter,
                                               input$num_samples)})
  
  # Create principal analysis plot
  plot_pca <- function(data, percentile, n_samples, x_axis, y_axis) {
    # Filter the data
    filtered_data <- custom_filter(data, percentile, n_samples)
    # Remove gene column
    selected_data <- select(filtered_data, -X)
    # Transpose data and perform PCA
    PCA <- prcomp(t(selected_data))
    # Extract principal components
    PCs <- as_tibble(PCA$x)
    # Get number from PC
    x_val <- substr(x_axis, 3, 4)
    y_val <- substr(y_axis, 3, 4)
    # Get explained variance
    xaxis_exp_var <- round(summary(PCA)$importance[2, as.numeric(x_val)] * 100)
    yaxis_exp_var <- round(summary(PCA)$importance[2, as.numeric(y_val)] * 100)
    # Manually create metadata
    stages <- c("early", "early", "early", "middle", "middle", "middle",
                "late", "late", "late")
    # Plot
    p <- ggplot(PCs, aes(x = get(x_axis), y = get(y_axis), color = stages)) +
      geom_point() + 
      labs(title = "Principal Component Analysis",
           x = paste(x_axis, ": ", xaxis_exp_var, "% variance", sep = ""),
           y = paste(y_axis, ": ", yaxis_exp_var, "% variance", sep = "")) +
      theme_bw()
    p
  }
  
  output$pca_plot <- renderPlot({plot_pca(load_counts_data(),
                                          input$counts_filter,
                                          input$num_samples,
                                          input$x_axis,
                                          input$y_axis)})
  
  
  
  # Code for differential expression analysis
  
  # Load de data
  load_de_data <- reactive({
    validate(
      need(input$dataset, "\nPlease select dataset in first tab"),
    )
    return(read.csv(paste("data/", input$dataset, "_de.csv",  sep = "", collapse = "")))
  })
  
  # Subset de data based on dropdown menu
  subset_de_table <- function(data, var) {
    name_combiner <- function(ex) {
      paste(var, ex, sep = "")
    }
    names_vect <- c(".logFC", ".logCPM", ".LR", ".PValue", ".FDR")
    result_names <- as.character(sapply(names_vect, name_combiner))
    result_names <- append(result_names, "ID", 0)
    dplyr::select(data, result_names)
  }
  
  # Draw table with filtered results
  draw_table <- function(dataf, slider, comparison) {
    dataf <- subset_de_table(dataf, comparison)
    pval <- paste(comparison, ".PValue", sep = "")
    df <- dataf[which(dataf[5] < 10**slider),]
    df[2] <- round(df[2], 3)
    df[3] <- round(df[3], 3)
    df[4] <- round(df[4], 3)
    df[6] <- round(df[6], 3)
    df
  }
  
  output$filtered_de_data <- DT::renderDataTable({draw_table(load_de_data(),
                                                            input$slider,
                                                            input$comp_select)})
  
  # Create the volcano plot
  volcano_plot <-
    function(dataf, x_name, y_name, slider, highlight, base, comparison) {
      # Subset data
      dataf <- subset_de_table(dataf, comparison)
      # Get column names
      x_name <- paste(comparison, x_name, sep = "")
      y_name <- paste(comparison, y_name, sep = "")
      pval <- paste(comparison, ".PValue", sep = "")
      # Plot data
      ggplot(dataf) +
        geom_point(aes(x=get(x_name), y=-log10(get(y_name)), color = get(pval) < 10**slider)) +
        scale_color_manual(name = paste('pvalue < 1e', as.character(slider), sep =""), values = setNames(c(highlight, base), c(T, F))) +
        theme_bw() +
        theme(legend.position="bottom") +
        xlab(as.character(x_name)) +
        ylab(as.character(y_name))
    }
  
  # Return volcano plot output
  output$volcano_plot <- renderPlot({volcano_plot(load_de_data(),
                                             input$x_name, input$y_name,
                                             input$slider, input$highlight,
                                             input$base, input$comp_select)})
  
  
  
  # Load sg counts data
  load_sg_counts_data <- reactive({
    validate(
      need(input$dataset, "\nPlease select a dataset in first tab"),
    )
    return(read.csv(paste("data/", input$dataset, "_counts.csv",  sep = "", collapse = "")))
  })
  
  # Add UI Element
  output$moreControls <- renderUI({
    validate(
      need(input$dataset, "")
    )
    selectInput("gene", "Select a gene of interest",
                choices = load_sg_counts_data()[1])
  })
  
  # Obtain single gene data for plotting
  get_plot_data <- function (data, gene) {
    # Extract row based on gene of interest
    ext_data <- as.numeric(data[data$X==gene,2:10])
    # Manually create stages
    stages <- factor(c("early", "early", "early", "middle", "middle", "middle", "late", "late", "late"))
    stages <- factor(stages, levels=c("early", "middle", "late"))
    # Create temporary data frame for plotting
    plot_data <- data.frame(ext_data, stages, colnames(data[2:10]))
    # Fix colname for samples
    colnames(plot_data)[3] <- "Samples"
    plot_data
  }
  
  plot_se_data <- function (data, gene, plot) {
    # Extract data for the gene of interest
    plot_data <- get_plot_data(data, gene)
    if (plot == "bar_plot") {
      # For bar plot
      p <- ggplot(data=plot_data, aes(x=Samples, y=ext_data, fill=stages)) +
        geom_bar(stat='identity') +
        theme_bw() +
        labs(title = paste("Expression of ", gene, sep = ""),
             y = "Gene Expression")
      p
    } else if (plot == "boxplot") {
      # For boxplot
      p <- ggplot(data=plot_data, aes(x=stages, y=ext_data, fill=stages)) +
        geom_boxplot() +
        theme_bw() +
        labs(title = paste("Expression of ", gene, sep = ""),
             y = "Gene Expression",
             x = "Stages")
      p
    } else if (plot == "violin_plot") {
      # For violin plot
      p <- ggplot(data=plot_data, aes(x=stages, y=ext_data, fill=stages)) +
        geom_violin() +
        theme_bw() +
        labs(title = paste("Expression of ", gene, sep = ""),
             y = "Gene Expression",
             x = "Stages")
      p
    } else if (plot == "beeswarm_plot") {
      # For beeswarm plot
      p<-ggplot(data=plot_data, aes(x=stages, y=ext_data, col=stages)) +
        geom_beeswarm(groupOnX = TRUE) +
        theme_bw() +
        labs(title = paste("Expression of ", gene, sep = ""),
             y = "Gene Expression",
             x = "Stages")
      p
    }
  }
  
  # Output plot for single gene examination
  output$single_gene_plot <- renderPlot({plot_se_data(load_sg_counts_data(),
                                                      input$gene,
                                                      input$plot_type)
  })
  
  output$text <- renderText(input$dataset)
}

# Run the application
shinyApp(ui = ui, server = server)