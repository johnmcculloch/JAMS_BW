#' plot_genome_taxonomy_donut(opt = opt, contigsdata = NULL)
#'
#' JAMSalpha function
#' @export

plot_genome_taxonomy_donut <- function(opt = opt, assemblystats = NULL, taxlevel = "LKT", isolatename = NULL, percentage_threshold = 99, ntop = 8, brewer_palette = 3){

    if (!is.null(opt)){
        assemblystats <- opt$assemblystats
        isolatename <- opt$prefix
    }
    genomedf <- subset(assemblystats, TaxLevel == taxlevel)
    rownames(genomedf) <- genomedf$Taxon

    genomesizeinfo <- paste(isolatename, "total genome size:", paste(sum(genomedf$ContigSum), "bp"))
    taxlvlnam <- switch(taxlevel, "Species" = "Species", "IS1" = "Strain (infraspecies)", "LKT" = "Contig Last Known Taxon")

    lvlmsg <- paste("", taxlvlnam)
    genomedf <- genomedf[ , c("Taxon", "ContigSum")]

    #Shamelessly based on https://r-graph-gallery.com/128-ring-or-donut-plot.html

    #Get genome length proportions
    genomedf$Percentage <- round(((genomedf$ContigSum / sum(genomedf$ContigSum)) * 100), 2)
    genomedf$CumPercentage <- cumsum(genomedf$Percentage)

    if (percentage_threshold <= min(genomedf$CumPercentage)){
        percentage_threshold <- min(genomedf$CumPercentage)
    }

    if (length(which(genomedf$CumPercentage <= percentage_threshold)) > ntop){
        rowsToKeep <- 1:ntop
    } else {
        rowsToKeep <- which(genomedf$CumPercentage <= percentage_threshold)
    }

    remainder <- genomedf[!(rownames(genomedf) %in% rowsToKeep), ]
    if (nrow(remainder) > 0){
        genomedf <- genomedf[rowsToKeep, ]
        genomedf[(nrow(genomedf) + 1), ] <- c(paste0("Other (n=", nrow(remainder), ")"), colSums(remainder[,2:3]), 100)
        genomedf$Percentage <- as.numeric(genomedf$Percentage)
        genomedf$CumPercentage <- as.numeric(genomedf$CumPercentage)
    }
    genomedf$ymax <- genomedf$CumPercentage
    genomedf$ymin <- c(0, head(genomedf$ymax, n=-1))

    genomedf$labelPosition <- (genomedf$ymax + genomedf$ymin) / 2
    genomedf$label <- paste(genomedf$Taxon, paste0(genomedf$Percentage, "%"), sep = "\n")

    genome_donut <- ggplot(genomedf, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = Taxon))
    genome_donut <- genome_donut + geom_rect()
    #genome_donut <- genome_donut + geom_label_repel(x = 3.5, aes(label = label, y = labelPosition), size = 2, min.segment.length = 0, box.padding = 3, segment.size = 0.5, segment.color = "black", segment.curvature = -0.05, segment.ncp = 1, segment.inflect = TRUE, max.overlaps = Inf)
    genome_donut <- genome_donut + geom_label_repel(x = 3.5, aes(label = label, y = labelPosition), size = 2, min.segment.length = 0, box.padding = 3, segment.size = 0.5, segment.color = "black", max.overlaps = Inf)

    genome_donut <- genome_donut + scale_fill_brewer(palette = brewer_palette) + scale_color_brewer(palette = brewer_palette)
    genome_donut <- genome_donut + coord_polar(theta = "y") + xlim(c(-1, 4))
    genome_donut <- genome_donut + theme_void() + theme(plot.title = element_text(size = 12),legend.position = "none", plot.background = element_rect(colour = "black", fill = NA, size = 1))+ labs(title = lvlmsg, caption = NULL)
    #genome_donut <- genome_donut + geom_text(x = -0.5, y = 0, size = 3.5, label = taxlevel)
    #genome_donut <- genome_donut + ggtitle(title = taxlevel)

    return(genome_donut)
}
