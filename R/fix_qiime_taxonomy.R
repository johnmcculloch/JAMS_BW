#' Convert qiime2 taxonomy table to one compatible with JAMS
#'
#'
#' @param taxtsv a "taxonomy.tsv" file created by qiime2.
#' @export
#' @examples
#' fix_qiime2_taxonomy("taxonomy.tsv")

fix_qiime2_taxonomy <- function(taxtsv) {
  tax <- read.table(
    taxtsv,
    header = TRUE,
    sep = "\t",
    row.names = 1,
    stringsAsFactors = FALSE
  )
  taxtable <- sapply(tax$Taxon, function(x) {
    x <- str_remove_all(x, fixed("["))
    x <- str_remove_all(x, fixed("]"))
    x <- str_replace(x, "D_0__", "d__")
    x <- str_replace(x, "D_1__", "p__")
    x <- str_replace(x, "D_2__", "c__")
    x <- str_replace(x, "D_3__", "o__")
    x <- str_replace(x, "D_4__", "f__")
    x <- str_replace(x, "D_5__", "g__")
    x <- str_replace(x, "D_6__", "s__")
    levs = str_split_fixed(x, fixed(";"), n = 7)
    levs
  })
  taxtable <- t(taxtable)
  rownames(taxtable) <- row.names(tax)
  colnames(taxtable) <- c("Domain",
                          "Phylum",
                          "Class",
                          "Order",
                          "Family",
                          "Genus",
                          "Species")
  taxtable[taxtable == ""] <- "Unclassified"
  taxtable[grepl("uncultured", taxtable)] <- "Unclassified"
  taxtable <- cbind(taxtable,
                    LKT = infer_LKT(data.frame(taxtable,
                                               stringsAsFactors = FALSE))$LKT)
  return(data.frame(taxtable))
}
