args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("A file path must be provided.")
}

library(ComplexHeatmap)
library(JAMS)

# Load your SummarizedExperiment object
# This is just a placeholder, adjust according to your actual data
ExpObj <- readRDS(args[1])

# Call the function to plot the heatmap
output_file <- "heatmap.png"
png(filename = output_file)
plot_relabund_heatmap(ExpObj = ExpObj)
dev.off()

cat(output_file)
