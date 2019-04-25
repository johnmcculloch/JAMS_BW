#' Creates Venn Diagram (2, 3, or 4 way) from metagenomeSeq object and returns list of features broken down from diagram
#'
#'
#' @param mgseqobj a metagenomeSeq object.
#' @param category string giving pData field to group samples by.
#' @param glomby string giving taxonomic level (from tax_table) to plot at.
#' @param samplesToKeep vector with samples to plot.
#' @param featuresToKeep vector with features to plot.
#' @param mgSeqnorm whether to normalize counts.
#' @param featmaxatleastPPM only count features with at least a certain PPM at max value.
#' @param featcutoff only count features with at least a certain ppm in a certain percentage of samples.
#' @param title plot title.
#' @param colours use list of colors.
#' @param labels to overwrite metagenomeseq names
#' @param pdf write plot as a pdf to a certain file.
#' @export
#' @examples
#' function plot_Venn(mgseqobj = NULL,category = NULL, glomby = NULL, samplesToKeep = NULL, featuresToKeep = NULL, mgSeqnorm = FALSE, featmaxatleastPPM = 0, featcutoff = c(0, 0), aggFun = rowMeans, title = "", colours= NULL, pdf = NULL, ...)

plot_Venn <- function (mgseqobj = NULL,category = NULL, glomby = NULL, samplesToKeep = NULL, featuresToKeep = NULL, mgSeqnorm = FALSE, featmaxatleastPPM = 0, featcutoff = c(0, 0),  aggFun = rowMeans, title = "", colours= NULL, labels = NULL, pdf = NULL, ...) {
  obj <- mgseqobj
  if (!(is.null(glomby))) {
    obj <- aggTax(obj,
                  lvl = glomby,
                  out = "MRexperiment",
                  norm = FALSE)
  }
  if (!(is.null(category))) {
    if (!is.factor(obj[[category]])) {
      category <- factor(obj[[category]])
    }
    obj <- aggSamp(obj, category, aggfun = aggFun)
  }
  if (ncol(obj) < 2 || ncol(obj) > 4) {
    print("Can only make Venn diagram of 2, 3, 4 categories/samples.")
    print(str_c("You have ", ncol(obj), "."))
    return()
  }
  obj <- filter_experiment(
    mgseqobj = obj,
    featmaxatleastPPM = featmaxatleastPPM,
    featcutoff = featcutoff,
    samplesToKeep = samplesToKeep,
    asPA = TRUE,
    mgSeqnorm = mgSeqnorm
  )
  vennMatrix <- data.frame(MRcounts(obj))
  results = list()
  if(!is.null(labels)) {
    colnames(obj) <- labels
  }
  if (!is.null(pdf)) {
    pdf(pdf, bg = "white")
  }
  if (ncol(obj) == 2) {
    results[["10"]] <- subset(vennMatrix, vennMatrix[, 1] == 1 &
                                vennMatrix[, 2] == 0)
    results[["01"]] <-
      subset(vennMatrix, vennMatrix[, 1] == 0 & vennMatrix[, 2] == 1)
    results[["11"]] <-
      subset(vennMatrix, vennMatrix[, 1] == 1 & vennMatrix[, 2] == 1)
    counts <- lapply(results, function(x) nrow(x)) %>% unlist()
    colorfulVennPlot::plotVenn2d(counts, colnames(obj),
                                 Colors = colours, Title = title)
  } else if (ncol(obj) == 3) {
    for (i in 0:1) {
      for (j in 0:1) {
        for (k in 0:1) {
          key = str_c(i,j,k)
          results[[key]] <- subset(vennMatrix, vennMatrix[, 1] == i & vennMatrix[, 2] == j &
                                     vennMatrix[, 3] == k)
        }
      }
    }
    counts <- lapply(results, function(x) nrow(x)) %>% unlist()
    colorfulVennPlot::plotVenn3d(counts, colnames(obj), Colors = colours,
                                 Title = title)
  } else if (ncol(obj) == 4) {
    for (i in 0:1) {
      for(j in 0:1) {
        for(k in 0:1) {
          for(l in 0:1) {
            key = str_c(i,j,k,l)
            results[[key]] <- subset(vennMatrix,vennMatrix[, 1] == i & vennMatrix[, 2] == j &
                                       vennMatrix[, 3] == k & vennMatrix[, 4] == l)
          }
        }
      }
    }
    counts <- lapply(results, function(x) nrow(x)) %>% unlist()
    colorfulVennPlot::plotVenn4d(counts, colnames(obj), Colors = colours,
                                 Title = title)
  }
  if (!is.null(pdf)) {
    dev.off()
  }
  return(results)
}
