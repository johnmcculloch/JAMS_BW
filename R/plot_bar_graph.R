#' Creates bar graph of relative taxonomic abundance from SummarizedExperiment object
#'
#'
#' @param ExpObj a SummarizedExperiment object.
#' @param glomby string giving taxonomic level (from tax_table) to plot at.
#' @param samplesToKeep vector with samples to plot.
#' @param featuresToKeep vector with features to plot.
#' @param total_units unit to normalize data to -- default percentage (100).
#' @param threshold minimum value to include as separate entry (below is Misc_Low_Abundance).
#' @param groupby colData category to use for labels
#' @param cat_pos position of category line (change if ugly in plot)
#' @param cat_text_size text size of category labels
#' @param legend_text_size text size of legend labels
#' @param border_color color of bar borders
#' @param cdict use JAMS colour dictionary
#' @param feature_cdict use named list of taxon colors
#' @param addtit optional string to add more text to title
#' @export
#' @examples
#' plot_bar_graph(ExpObj = mrexp1, glomby = "Phylum", groupby = "Group")

plot_bar_graph <- function(ExpObj = NULL, glomby = NULL, samplesToKeep = NULL, featuresToKeep = NULL, groupby = NULL, colourby = NULL, subsetby = NULL, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL, PPM_normalize_to_bases_sequenced = FALSE, ignoreunclassified = FALSE, total_units = 100, threshold = 2, cat_pos = -2, cat_text_size = 3, legend_text_size = 4, border_color = "white", cdict = NULL, feature_cdict = NULL, addtit = NULL, grid = TRUE, class_to_ignore = "N_A", ...) {

    require(reshape2)
    require(Polychrome)

    valid_vars <- c(groupby, subsetby)[which(!is.na(c(groupby, subsetby)))]

    #Vet experiment object
    obj <- ExpObjVetting(ExpObj = ExpObj, samplesToKeep = samplesToKeep, featuresToKeep = featuresToKeep, glomby = glomby, variables_to_fix = valid_vars, class_to_ignore = class_to_ignore)

    analysis <- metadata(obj)$analysis
    if (!is.null(glomby)){
        analysisname <- glomby
    } else {
        analysisname <- analysis
    }

    presetlist <- declare_filtering_presets(analysis = analysis, applyfilters = applyfilters, featcutoff = featcutoff, GenomeCompletenessCutoff = GenomeCompletenessCutoff, PctFromCtgscutoff = PctFromCtgscutoff)

    if (!(is.null(subsetby))){
        subset_points <- sort(unique(colData(obj)[, which(colnames(colData(obj)) == subsetby)]))
    } else {
        subset_points <- "none"
    }

    #Create list vector to hold plots
    gvec <- list()
    plotnum <- 1
    bartitbase <- paste("Relative abundances of", analysisname)

    for (sp in 1:length(subset_points)){

        if (!(is.null(subsetby))){
            sp_samplesToKeep <- rownames(colData(obj))[which(colData(obj)[ , subsetby] == subset_points[sp])]
            flog.info(paste("Plotting within", subset_points[sp]))
            subsetname <- subset_points[sp]
            bartit <- paste(bartitbase, "within", subset_points[sp])
        } else {
            sp_samplesToKeep <- rownames(colData(obj))
            subsetname <- "no_sub"
            bartit <- bartitbase
        }
        bartit <- paste(c(bartit, presetlist$filtermsg), collapse = "\n")

        currobj <- filter_experiment(ExpObj = obj, featcutoff = presetlist$featcutoff, samplesToKeep = sp_samplesToKeep, featuresToKeep = featuresToKeep, asPPM = TRUE, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, GenomeCompletenessCutoff = presetlist$GenomeCompletenessCutoff, PctFromCtgscutoff = presetlist$PctFromCtgscutoff)

        currpt <- as.data.frame(colData(currobj))

        if (PPM_normalize_to_bases_sequenced == TRUE){
            bartit <- paste(c(bartit, "Normalized to total number of bases sequenced in sample"), collapse = "\n")
        } else {
            bartit <- paste(c(bartit, "Normalized to number of bases for analysis in sample"), collapse = "\n")
        }

        #Get counts matrix
        countmat <- as.matrix(assays(currobj)$BaseCounts)

        #Protect against rows with empty data
        rowsToKeep <- which(rowSums(countmat) > 0 & rownames(countmat) != "")
        countmat <- countmat[rowsToKeep, ]

        if (ignoreunclassified == TRUE){
            dunno <- c(paste(analysis, "none", sep = "_"), "LKT__d__Unclassified", "LKT__Unclassified")
            rowsToKeep <- which(!(rownames(countmat) %in% dunno))
            countmat <- countmat[rowsToKeep, ]
        }

        #Re-normalize to percentage
        for (colm in 1:ncol(countmat)){
            countmat[ , colm] <- (countmat[ , colm] / sum(countmat[ , colm])) * total_units
        }

        if (threshold > 0) {
            new_counts <- countmat[rowMeans(countmat) > threshold, ]
            misc_counts <- total_units - colSums(new_counts)
            new_counts <- rbind(new_counts, Misc_Low_Abundance = misc_counts)
        } else {
            new_counts <- countmat
        }

        mdf <- reshape2::melt(new_counts)
        colnames(mdf) <- c(analysisname, "Sample", "Abundance")

        fill_sum <- aggregate(as.formula(paste("Abundance ~", analysisname)), mdf, sum)
        mdf[ , analysisname] <- factor(mdf[ , analysisname], levels = fill_sum[order(-fill_sum$Abundance), analysisname])

        mdf$Sample <- factor(mdf$Sample, levels = colnames(new_counts))

        p <- ggplot(mdf, aes_string(x = "Sample", y = "Abundance", fill = analysisname))
        p <- p + geom_bar(stat = "identity", position = "stack", color = border_color, size = 0)
        if(!is.null(addtit)) {
          bartit <- str_c(bartit, " ", addtit)
        }
        p <- p + ggtitle(bartit)
        p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
        p <- p + theme(plot.margin = unit(c(0, 1, 1, 0.5), "lines"))
        p <- p + theme(strip.background = element_rect(fill = "white", colour = NA))
        p <- p + theme(panel.background = element_rect(fill = "white", colour = "white"))
        p <- p + theme(plot.margin = margin(c(1, 1, 3, 1), unit = "lines"))
        p <- p + theme(axis.line = element_blank())
        p <- p + theme(legend.position = "bottom")
        p <- p + theme(legend.text = element_text(colour = "black", size = legend_text_size, face="italic"))
        p <- p + theme(legend.title = element_text(colour = "black", size = legend_text_size))
        p <- p + theme(title = element_text(size = 8))

        if(!(is.null(cdict))) {
            ct <- cdict[[colourby]]
            groupcols <- setNames(as.character(ct$Colour), as.character(ct$Name))
            p <- p + scale_color_manual(values = groupcols)
        }

        if (is.null(feature_cdict)){
            feature_cdict <- createPalette(nrow(new_counts),  c("#ff0000", "#00ff00", "#0000ff"))
            names(feature_cdict) <- rownames(new_counts)
        }

        p <- p + scale_fill_manual(values = feature_cdict)

        if (!is.null(groupby)) {
            cat_values <- as.factor(currpt[ , groupby])
            for (level in levels(cat_values)) {
                startx <- min(which(cat_values==level))
                endx <- max(which(cat_values==level))
                p <- p + annotate("segment", x=startx, xend=endx, y = cat_pos, yend = cat_pos)
                p <- p + annotate("text", label = level, x = startx + (endx - startx) / 2, y = 2 * cat_pos, cex = cat_text_size)
            }
        }
        gvec[[plotnum]] <- p
        plotnum <- plotnum + 1
    }

    return(gvec)
}
