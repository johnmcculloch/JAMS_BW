#' Creates bar graph of relative taxonomic abundance from SummarizedExperiment object
#'
#' @param ExpObj JAMS-style SummarizedExperiment object

#' @param glomby String giving the taxonomic level at which to agglomerate counts. This argument should only be used with taxonomic SummarizedExperiment objects. When NULL (the default), there is no agglomeration

#' @param absolute Requires a logical value. If set to true, will not normalize to 100% (or whatever) but will show absolute abundance instead. Default is FALSE.

#' @param samplesToKeep Vector with sample names to keep. If NULL, all samples within the SummarizedExperiment object are kept. Default is NULL.

#' @param featuresToKeep Vector with feature names to keep. If NULL, all features within the SummarizedExperiment object are kept. Default is NULL. Please note that when agglomerating features with the glomby argument (see above), feature names passed to featuresToKeep must be post-agglomeration feature names. For example, if glomby="Family", featuresToKeep must be family names, such as "f__Enterobacteriaceae", etc.

#' @param groupby String specifying the metadata variable name. colData category to use for labels.

#' @param label_samples Requires a logical value. If set to TRUE, will include labels for sample names on graph. Default is FALSE, and will not include sample names.

#' @param colourby String specifying the metadata variable name for colouring in samples. If NULL, all samples will be black. Default is NULL.

#' @param subsetby String specifying the metadata variable name for subsetting samples. If passed, multiple plots will be drawn, one plot for samples within each different class contained within the variable.  If NULL, data is not subset. Default is NULL.

#' @param applyfilters Optional string specifying filtration setting "combos", used as a shorthand for setting the featcutoff, GenomeCompletenessCutoff, minl2fc and minabscorrcoeff arguments in JAMS plotting functions. If NULL, none of these arguments are set if not specified. Permissible values for applyfilters are "light", "moderate" or "stringent". The actual values vary whether the SummarizedExperiment object is taxonomical (LKT) or not. For a taxonomical SummarizedExperiment object, using "light" will set featcutoff=c(50, 5), GenomeCompletenessCutoff=c(5, 5), minl2fc=1, minabscorrcoeff=0.4; using "moderate" will set featcutoff=c(250, 15), GenomeCompletenessCutoff=c(10, 5), minl2fc=1, minabscorrcoeff=0.6; and using "stringent" will set featcutoff=c(2000, 15), GenomeCompletenessCutoff=c(30, 10), minl2fc=2, minabscorrcoeff=0.8. For non-taxonomical (i.e. functional) SummarizedExperiment objects, using "light" will set featcutoff=c(0, 0), minl2fc=1, minabscorrcoeff=0.4; using "moderate" will set featcutoff=c(5, 5), minl2fc=1, minabscorrcoeff=0.6; and using "stringent" will set featcutoff=c(50, 15), minl2fc=2.5, minabscorrcoeff=0.8. When using applyfilters, one can still set one or more of featcutoff, GenomeCompletenessCutoff, minl2fc and minabscorrcoeff, which will then take the user set value in lieu of those set by the applyfilters shorthand. Default is light.

#' @param featcutoff Requires a numeric vector of length 2 for specifying how to filter out features by relative abundance. The first value of the vector specifies the minimum relative abundance in Parts per Million (PPM) and the second value is the percentage of samples which must have at least that relative abundance. Thus, passing c(250, 10) to featcutoff would filter out any feature which does not have at least 250 PPM (= 0.025 percent) of relative abundance in at least 10 percent of all samples being plot. Please note that when using the subsetby option (q.v.) to automatically plot multiple plots of sample subsets, the featcutoff parameters are applied within the subset. The default is c(0, 0), meaning no feature is filtered. If NULL is passed, then the value defaults to c(0, 0). See also applyfilters for a shorthand way of applying multiple filtration settings.

#' @param GenomeCompletenessCutoff Requires a numeric vector of length 2 for specifying how to filter out features by genome completeness. This is, of course, only applicble for taxonomic shotgun SummarizedExperiment objects. When passed to non-taxonomic shotgun SummarizedExperiment objects, GenomeCompletenessCutoff will be ignored. The first value of the vector specifies the minimum genome completeness in percentage  and the second value is the percentage of samples which must have at least that genome completeness. Thus, passing c(50, 5) to GenomeCompletenessCutoff would filter out any taxonomic feature which does not have at least 50 percent of genome completeness in at least 5 percent of all samples being plot. Please note that when using the subsetby option (q.v.) to automatically plot multiple plots of sample subsets, the GenomeCompletenessCutoff parameters are applied within the subset. The default is c(0, 0), meaning no feature is filtered. If NULL is passed, then the value defaults to c(0, 0). See also applyfilters for a shorthand way of applying multiple filtration settings.

#' @param PPM_normalize_to_bases_sequenced Requires a logical value. Non-filtered JAMS feature counts tables (the BaseCounts assay within SummarizedExperiment objects) always includes unclassified taxonomical features (for taxonomical SummarizedExperiment objects) or unknown/unattributed functional features (for non-taxonomical SummarizedExperiment objects), so the relative abundance for each feature (see normalization) will be calculated in Parts per Million (PPM) by dividing the number of bases covering each feature by the sum of each sample column **previous to any filtration**. Relative abundances are thus representative of the entirety of the genomic content for taxonomical objects, whereas for non-taxonomical objects, strictly speaking, it is the abundance of each feature relative to only the coding regions present in the metagenome, even if these are annotationally unatributed. In other words, intergenic regions are not taken into account. In order to relative-abundance-normalize a **non-taxonomical** SummarizedExperiment object with the total genomic sequencing content, including non-coding regions, set PPM_normalize_to_bases_sequenced = TRUE. Default is FALSE.

#' @param total_units Numerical value specifying the unit to normalize data to. Default percentage is 100.

#' @param threshold Numberical value specifying the minimum value to include as separate entry (below is Misc_Low_Abundance). Default is 2.

#' @param cat_pos Numerical value that specifies the position of the category line (change if ugly in plot, otherwise mess with at your own risk). Default is -2.

#' @param cat_text_size Numerical value setting the text size of category labels. Default is 3.

#' @param legend_text_size Numerical value setting the text size of legend labels. Default is 4.

#' @param border_color String that specifies the color of bar borders. Default is white.

#' @param feature_cdict Vector containing a dictionary of colours with a corresponding taxon. Allows you to manually define the colours for taxonomy. Default is NULL.

#' @param addtit Optional string with text to append to main title. Default is NULL.

#' @param grid Requires a logical value. If set to FALSE, background will be one solid color within the plot, rather than include a grid behind the plot. Default is TRUE, meaning that the background will display a grid.

#' @param ignoreunclassified Requires a logical value. If set to TRUE, for taxonomical SummarizedExperiment objects, the feature "LKT__Unclassified" will be omitted from being shown. In the case of non-taxonomical SummarizedExperiment objects, the completely unannotated features will be omitted. For example, for an ECNumber SummarizedExperiment object, genes *without* an Enzyme Commission Number annotation (feature "EC_none") will not be shown. Statistics are, however, computed taking the completely unclassifed feature into account, so p-values will not change.

#' @param class_to_ignore String or vector specifying any classes which should lead to samples being excluded from the comparison within the variable passed to compareby. Default is N_A. This means that within any metadata variable passed to compareby containing the "N_A" string within that specific variable, the sample will be dropped from that comparison.

#' @examples plot_bar_graph(ExpObj = expvec$LKT, glomby = "Phylum", groupby = "Group")

#' @export
plot_bar_graph <- function(ExpObj = NULL, glomby = NULL, absolute = FALSE, samplesToKeep = NULL, featuresToKeep = NULL, groupby = NULL, label_samples = FALSE, colourby = NULL, subsetby = NULL, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, PPM_normalize_to_bases_sequenced = FALSE, total_units = 100, threshold = 2, cat_pos = -2, cat_text_size = 3, legend_text_size = 4, border_color = "white", cdict = NULL, feature_cdict = NULL, addtit = NULL, grid = TRUE, ignoreunclassified = FALSE, class_to_ignore = "N_A", ...) {

    require(reshape2)
    require(Polychrome)

    valid_vars <- c(groupby, subsetby)[which(!is.na(c(groupby, subsetby)))]

    #Hardwire PctFromCtgscutoff, as this should never be used without filtering because of huge amounts of false positives when evaluating taxonomic information from unassembled reads. The use of classifying unassembled reads is deprecated in JAMS and the default is to NOT classify unassembled reads, so this is usually not an issue.
    PctFromCtgscutoff <- c(70, 50)
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

    presetlist <- declare_filtering_presets(analysis = analysis, applyfilters = applyfilters, featcutoff = featcutoff, GenomeCompletenessCutoff = GenomeCompletenessCutoff, PctFromCtgscutoff = PctFromCtgscutoff)

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

        if (absolute) {
          currobj <- obj
        } else {
            currobj <- filter_experiment(ExpObj = obj, featcutoff = presetlist$featcutoff, samplesToKeep = sp_samplesToKeep, featuresToKeep = featuresToKeep, asPPM = TRUE, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, GenomeCompletenessCutoff = presetlist$GenomeCompletenessCutoff, PctFromCtgscutoff = presetlist$PctFromCtgscutoff)
        }

        currpt <- as.data.frame(colData(currobj))

        if (PPM_normalize_to_bases_sequenced == TRUE){
            bartit <- paste(c(bartit, "Normalized to total number of bases sequenced in sample"), collapse = "\n")
        } else if (!absolute) {
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
