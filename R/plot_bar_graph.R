#' Creates bar graph of relative taxonomic abundance from SummarizedExperiment object.

#' @export

plot_bar_graph <- function(ExpObj = NULL, glomby = NULL, absolute = FALSE, samplesToKeep = NULL, featuresToKeep = NULL, groupby = NULL, label_samples = FALSE, colourby = NULL, subsetby = NULL, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, PPM_normalize_to_bases_sequenced = FALSE, total_units = 100, threshold = 2, cat_pos = -2, cat_text_size = 3, legend_text_size = 4, border_color = "white", cdict = NULL, feature_cdict = NULL, addtit = NULL, grid = TRUE, ignoreunclassified = FALSE, class_to_ignore = "N_A", ...) {

    require(reshape2)
    require(Polychrome)

    valid_vars <- c(groupby, subsetby)[which(!is.na(c(groupby, subsetby)))]

    transp <- TRUE

    #Vet experiment object
    obj <- ExpObjVetting(ExpObj = ExpObj, samplesToKeep = samplesToKeep, featuresToKeep = featuresToKeep, glomby = glomby, variables_to_fix = valid_vars, class_to_ignore = class_to_ignore)
    if (!is.null(groupby)) {
      if (!is.null(obj[[groupby]])) {
        obj <- obj[,order(obj[[groupby]])]
      } else {
        flog.error(str_c(groupby, " is not a valid column of your metadata"))
        groupby <- NULL
      }

    }
    analysis <- metadata(obj)$analysis
    if (!is.null(glomby)){
        analysisname <- glomby
    } else {
        analysisname <- analysis
    }

    presetlist <- declare_filtering_presets(analysis = analysis, applyfilters = applyfilters, featcutoff = featcutoff, GenomeCompletenessCutoff = GenomeCompletenessCutoff)

    if (!(is.null(subsetby))){
        subset_points <- sort(unique(colData(obj)[, which(colnames(colData(obj)) == subsetby)]))
    } else {
        subset_points <- "none"
    }

    #Create list vector to hold plots
    gvec <- list()
    plotnum <- 1
    if (absolute) {
      bartitbase <- paste("Absolute abundances of", analysisname)
    } else {
      bartitbase <- paste("Relative abundances of", analysisname)
    }

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

        currobj <- filter_experiment(SEobj = obj, featcutoff = presetlist$featcutoff, samplesToKeep = sp_samplesToKeep, featuresToKeep = featuresToKeep, normalization = "relabund", PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, GenomeCompletenessCutoff = presetlist$GenomeCompletenessCutoff)

        currpt <- as.data.frame(colData(currobj))

        if (PPM_normalize_to_bases_sequenced == TRUE){
            bartit <- paste(c(bartit, "Normalized to total number of bases sequenced in sample"), collapse = "\n")
        } else if (!absolute) {
            bartit <- paste(c(bartit, "Normalized to number of bases for analysis in sample"), collapse = "\n")
        }

        #Get counts matrix
        if (absolute){
            countmat <- as.matrix(assays(currobj)$BaseCounts)
        } else {
            countmat <- as.matrix(assays(currobj)$PPM)
        }

        #Protect against rows with empty data
        rowsToKeep <- which(rowSums(countmat) > 0 & rownames(countmat) != "")
        countmat <- countmat[rowsToKeep, ]

        if (ignoreunclassified == TRUE){
            dunno <- c(paste(analysis, "none", sep = "_"), "LKT__d__Unclassified", "LKT__Unclassified")
            rowsToKeep <- which(!(rownames(countmat) %in% dunno))
            countmat <- countmat[rowsToKeep, ]
        }

        #Re-normalize to percentage
        if (!absolute) {
          for (colm in 1:ncol(countmat)){
              countmat[ , colm] <- (countmat[ , colm] / sum(countmat[ , colm])) * total_units
          }
        }
        if (threshold > 0) {
            new_counts <- countmat[rowMeans(countmat) > threshold,, drop=FALSE ]
            if (absolute) {
              misc_counts <- colSums(countmat) - colSums(new_counts)
            } else {
              misc_counts <- total_units - colSums(new_counts)
            }
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
        p <- p + theme(strip.background=element_rect(fill = "white", colour = NA))
        p <- p + theme(panel.background = element_rect(fill = "white", colour = "white"))
        p <- p + theme(plot.margin = margin(c(1, 1, 3, 1), unit = "lines"))
        p <- p + theme(axis.line = element_blank())
        p <- p + theme(legend.position = "bottom")
        p <- p + theme(legend.text = element_text(colour = "black", size = legend_text_size, face="italic"))
        p <- p + theme(legend.title = element_text(colour = "black", size = legend_text_size))
        p <- p + theme(title = element_text(size = 8))
        if (label_samples) {
          p <- p +  theme(axis.text.x = element_text(angle=90))
        }
        if(!(is.null(cdict))) {
            ct <- cdict[[colourby]]
            groupcols <- setNames(as.character(ct$Colour), as.character(ct$Name))
            p <- p + scale_color_manual(values = groupcols)
        }
        if (is.null(feature_cdict)){
            feature_cdict <- createPalette(nrow(new_counts),  c("#ff0000", "#00ff00", "#0000ff"))
            names(feature_cdict) <- rownames(new_counts)
        }
	feature_cdict <- feature_cdict[rownames(new_counts)] # remove entries not on graph

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
        return(p)
        gvec[[plotnum]] <- p
        plotnum <- plotnum + 1
    }

    return(gvec)
}
