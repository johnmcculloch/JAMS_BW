#' Creates bar graph of relative taxonomic abundance from SummarizedExperiment object
#'
#'
#' @param mgseqobj a SummarizedExperiment object.
#' @param glomby string giving taxonomic level (from tax_table) to plot at.
#' @param samplesToKeep vector with samples to plot.
#' @param featuresToKeep vector with features to plot.
#' @param units unit to normalize data to -- default percentage (100).
#' @param threshold minimum value to include as separate entry (below is Misc_Low_Abundance).
#' @param title plot title
#' @param category colData category to use for labels
#' @param cat_pos position of category line (change if ugly in plot)
#' @param cat_text_size text size of category labels
#' @param border_color color of bar borders
#' @param cdict use JAMS colour dictionary
#' @param colors use named list of taxon colors
#' @export
#' @examples
#' plot_bar_graph(mrexp1, glomby = "Phylum", category = "Group")

plot_bar_graph<-function(mgseqobj=NULL, glomby=NULL, samplesToKeep=NULL, featuresToKeep=NULL, units=100, threshold = 2, title="", category=NULL, cat_pos = -2, cat_text_size=3, border_color="white", cdict=NULL, colors=NULL, grid = TRUE, ...) {

    require(reshape2)

    #Get appropriate object to work with
    obj<-mgseqobj
    tax_name <- "Taxon"

    #Exclude samples and features if specified
    if(!(is.null(samplesToKeep))){
        obj<-obj[, samplesToKeep]
    }

    if(!(is.null(featuresToKeep))){
        obj<-obj[featuresToKeep, ]
    }

    #Aggregate if required
    if(!(is.null(glomby))){
        obj <- agglomerate_features(obj, glomby)
        tax_name <- glomby
    }

    #Define analysis type
    analysis<-attr(obj, "analysis")

    new_counts <- apply(assay(obj), 2, function(x) units * x/(sum(x) + 0.01))
    if (threshold > 0) {
      new_counts <- new_counts[rowMeans(new_counts) > threshold, ]
      misc_counts <- units - colSums(new_counts)
      new_counts <- rbind(new_counts, Misc_Low_Abundance = misc_counts)
    }
    mdf <- melt(new_counts)
    colnames(mdf) <- c(tax_name,"Sample","Abundance")
    fill_sum <- aggregate(as.formula(paste("Abundance ~", tax_name)),
                          mdf, sum)
    mdf[, tax_name] <- factor(mdf[, tax_name],
                              levels =
                                fill_sum[order(-fill_sum$Abundance),
                                                         tax_name])
    mdf$Sample <- factor(mdf$Sample, levels = colnames(new_counts))
    p <- ggplot(mdf, aes_string(x = "Sample", y = "Abundance",
                                fill = tax_name))
    p <- p + geom_bar(stat = "identity", position = "stack",
                      color = border_color, size = 0)
    p <- p + theme(axis.title.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())
    p <- p + theme(plot.margin = unit(c(0, 1, 1, 0.5), "lines"))
    p <- p + theme(strip.background = element_rect(fill = "white",
                                                   colour = NA))

    p <- p + theme(panel.background = element_rect(fill = "white",
                                                   colour = "white"))
    p <- p + theme(plot.margin = margin(c(1, 1, 3, 1), unit = "lines"))
    p <- p + theme(axis.line = element_blank())
    p <- p + ggtitle(title)
    if(!(is.null(cdict))) {
      ct<-cdict[[colourby]]
      groupcols<-setNames(as.character(ct$Colour), as.character(ct$Name))
      p <- p + scale_color_manual(values = groupcols)
    }
    if (!is.null(colors)) {
      p <- p + scale_fill_manual(values = colors)
    }
    if (!is.null(category)) {
      cat_values <- as.factor(colData(obj)[,category])
      for(level in levels(cat_values)) {
        startx <- min(which(cat_values==level))
        endx <- max(which(cat_values==level))
        p <- p + annotate("segment", x=startx,
                          xend=endx, y = cat_pos,
                          yend = cat_pos)
        p <- p + annotate("text", label = level,
                   x = startx + (endx - startx)/2,
                   y = 2*cat_pos, cex = cat_text_size)
      }
    }
    return(p)
}
