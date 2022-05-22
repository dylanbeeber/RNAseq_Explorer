# RNA Seq Explorer

This RNA seq explorer Shiny App serves as an example of the different ways to visualize and interperet RNA seq data in RShiny. In this app, we examine data from the soybean data set provided by BigPint Bioconductor package. To run this app, simply ensure that all package dependencies are met (listed at the top of the script), and run the script. 

In the app, first select the soybean data set in the initial tab. Download the counts matrix and differential expression data set using the download buttons within the page. These will download as .csv files that will serve as the input for the subsequent tabs. The second tab is the counts matrix explorer that takes the counts matrix file as input. The third tab allows for visualization of the differential expression data from the differential expression (de) file. Finally, the fourt tabs allows for analysis of patterns of expression in individual genes, this tab uses the counts matrix file as input.
