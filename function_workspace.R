library(tidyverse)
library(gplots)

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
    mutate(replicate = as.integer(str_extract(samples, "(?<=\\.)."))) %>%
    # Add column for the growth stage
    mutate(stage = encoding_dict[str_extract(samples, "(?<=S).")])
  # Convert stage column to factor
  metadata$stage <- as.factor(metadata$stage)
  metadata
}

test <- read.csv("se_soybean_cn_sub_counts.csv")

parameter_check <- function(x){
  sum((x) > 1)
}
check <- apply(test[,2:10], parameter_check, MARGIN = 1) >= 1
filtered <- dplyr::filter(test, check == TRUE)


selected_data <- dplyr::select(test, -X)
# Calculate mean counts
mean_counts <- apply(selected_data, mean, MARGIN = 1)
# Calculate variance
variance_counts <- apply(selected_data, var, MARGIN = 1)

percentile_val <- quantile(variance_counts, .5)

check <- variance_counts >= percentile_val
filtered <- dplyr::filter(test, check == TRUE)

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

create_mean_var_scatter <- function(data, percentile, n_samples) {
  # Filter data
  filtered_data <- custom_filter(data, percentile, n_samples)
  # Remove gene column
  selected_data <- dplyr::select(filtered_data, -X)
  # Calculate mean counts
  mean_counts <- apply(selected_data, mean, MARGIN = 1)
  # Calculate variance
  variance_counts <- apply(selected_data, var, MARGIN = 1)
  # Rank mean values
  ranked_means <- rank(mean_counts)
  # Create tibble for the variance/mean data
  temp <- tibble(mean_counts, variance_counts,ranked_means)
  # Plot
  p <- ggplot(temp, aes(x = ranked_means, y = variance_counts)) + 
    geom_point() +
    geom_smooth(method='gam', formula = y ~ s(x, bs = "cs")) +
    labs(title = "Rank vs. Variance",
         x = "Rank(Mean)",
         y = "Variance") +
    scale_y_continuous(trans = 'log10')
  print(p)
}

plot_pca <- function(data, percentile, n_samples) {
  # Filter the data
  filtered_data <- custom_filter(data, percentile, n_samples)
  # Remove gene column
  selected_data <- select(filtered_data, -X)
  # Transpose data and perform PCA
  PCA <- prcomp(t(selected_data))
  # Extract principal components
  PCs <- as_tibble(PCA$x)
  # Get explained variance
  PC1_exp_var <- round(summary(PCA)$importance[2, 1] * 100)
  PC2_exp_var <- round(summary(PCA)$importance[2, 2] * 100)
  stages <- c("early", "early", "early", "middle", "middle", "middle",
    "late", "late", "late")
  # Plot
  p <- ggplot(PCs, aes(x = PC1, y = PC2, color = stages)) +
    geom_point() + 
    labs(title = "Principal Component Analysis",
         x = paste("PC1: ", PC1_exp_var, "% variance", sep = ""),
         y = paste("PC2: ", PC2_exp_var, "% variance", sep = ""))
  p
}

var <- "S1_S2"
name_combiner <- function(ex) {
  paste(var, ex, sep = "")
}

result_names <- sapply(ex, name_combiner)
result_names[1]

subset_de_table <- function(data, var) {
  name_combiner <- function(ex) {
    paste(var, ex, sep = "")
  }
  names_vect <- c(".logFC", ".logCPM", ".LR", ".PValue", ".FDR")
  result_names <- as.character(sapply(names_vect, name_combiner))
  result_names <- append(result_names, "ID", 0)
  select(data, result_names)
}

hello <- subset_de_table("S1_S2")
select(test_de, hello)

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

heatmap.2(test_matrix, trace = "none", ylab="Genes",
          ColSideColors=rep(c("green","blue","purple"), each=3),
          xlab="Samples", labRow = FALSE, key = TRUE)


heatmap.2(test_matrix, trace = "none")

placeholder = geom_bar(stat='identity')

gene = "hello"


get_plot_data <- function (data, gene) {
  # Extract row based on gene of interest
  ext_data <- as.numeric(data[data$X==gene,2:10])
  # Manually create stages
  stages <- factor(c("early", "early", "early", "middle", "middle", "middle", "late", "late", "late"))
  # Create temporary data frame for plotting
  plot_data <- data.frame(ext_data, stages, colnames(test[2:10]))
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
      labs(title = paste("Expression of ", gene, sep = ""),
           y = "Gene Expression")
    p
  } elif (plot == "boxplot") {
    # For boxplot
    p <- ggplot(data=plot_data, aes(x=stages, y=ext_data, fill=stages)) +
      geom_boxplot() +
      labs(title = paste("Expression of ", gene, sep = ""),
           y = "Gene Expression",
           x = "Stages")
    p
  } elif (plot == "violin_plot") {
    # For violin plot
    p <- ggplot(data=plot_data, aes(x=stages, y=ext_data, fill=stages)) +
      geom_violin() +
      labs(title = paste("Expression of ", gene, sep = ""),
           y = "Gene Expression",
           x = "Stages")
    p
  } elif (plot == "beeswarm_plot") {
    # For beeswarm plot
    p<-ggplot(data=plot_data, aes(x=stages, y=ext_data, col=stages)) +
      geom_beeswarm() +
      labs(title = paste("Expression of ", gene, sep = ""),
           y = "Gene Expression",
           x = "Stages")
    p
  }
}

stages <- factor(c("early", "early", "early", "middle", "middle", "middle", "late", "late", "late"))
stages <- factor(stages, levels=c("early", "middle", "late"))


selected_data <- dplyr::select(test, -X)
mean_counts <- apply(selected_data, mean, MARGIN = 1)
variance_counts <- apply(selected_data, var, MARGIN = 1)
ranked_means <- rank(mean_counts)
check <- test$X %in% new_data$X
temp <- tibble(mean_counts, variance_counts, ranked_means, check)

filtered_data <- custom_filter(test, .5, 5)
excluded_data <- dplyr::anti_join(test, new_data)

check <- test$X %in% new_data$X

