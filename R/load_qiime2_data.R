#' load qiime2 files into metagenomeSeq object
#' 
#' 
#' @param otu_table a relative abundance otu file created by qiime2
#' @param tax_table a taxonomy file created by qiime2
#' @param metadata_table a metatable file
#' @param total total number of desired reads per sample
#' @export
#' @examples
#' load_qiime2_data("otu.tsv", "taxonomy.tsv", "metadata.tsv")

load_qiime2_data <- function(otu_table, tax_table, metadata_table, total = 100000) {
  counts <- round(read.table(otu_table, header=TRUE, skip=1, 
                             comment.char = '%', 
                             row.names = 1, sep="\t") * total)
  tax <- fix_qiime2_taxonomy(tax_table)
  metadata <- read.table(metadata_table, header=TRUE, 
                         comment.char = '%', row.names = 1, sep="\t")
  counts <- counts[,row.names(metadata)]
  counts <- counts[order(rowSums(counts), decreasing = TRUE),]
  tax <- tax[row.names(counts),]
  mrexp <- newMRexperiment(counts, AnnotatedDataFrame(metadata), 
                         AnnotatedDataFrame(tax))
  attr(mrexp, "analysis") <- "LKT"
  return(mrexp)
}
