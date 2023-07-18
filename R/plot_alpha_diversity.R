#' Generates relative abundance plots per feature annotated by the metadata using as input a SummarizedExperiment object


#' @param ExpObj JAMS-style SummarizedExperiment object

#' @param measures String giving the alpha diversity measurements used to quantify the alpha diversity within each category. Default includes Observed, Chao1, Shannon, Simpson, InvSimpson, and GeneCount.

#' @param stratify_by_kingdoms Requires a logical value. If TRUE, will concatenate all of the taxonomical features to the kingdom level and create individual plots for each of the kingdoms for each measure specified. If FALSE, will forgo the kingdom concatenation. Default is TRUE.

#' @param glomby String giving the taxonomic level at which to agglomerate counts. This argument should only be used with taxonomic SummarizedExperiment objects. When NULL (the default), there is no agglomeration

#' @param samplesToKeep Vector with sample names to keep. If NULL, all samples within the SummarizedExperiment object are kept. Default is NULL.

#' @param featuresToKeep Vector with feature names to keep. If NULL, all features within the SummarizedExperiment object are kept. Default is NULL. Please note that when agglomerating features with the glomby argument (see above), feature names passed to featuresToKeep must be post-agglomeration feature names. For example, if glomby="Family", featuresToKeep must be family names, such as "f__Enterobacteriaceae", etc.

#' @param subsetby String specifying the metadata variable name for subsetting samples. If passed, multiple plots will be drawn, one plot for samples within each different class contained within the variable.  If NULL, data is not subset. Default is NULL.

#' @param compareby  String specifying the metadata variable name for grouping samples. This will define which metadata variable grouping to calculate PERMANOVA p-value. If not specified, and argument permanova is set to TRUE, (see permanova), the compareby argument will be set by colourby or shapeby. If these latter two are also NULL, and permanova is TRUE, permanova will be set to FALSE. Default is NULL.

#' @param compareby_order String or vector specifying the order in which to compare by, if this order is different than alphabetical order of the compareby parameter. Default is NULL.

#' @param colourby String specifying the metadata variable name for colouring the lines of the boxes for the samples. If NULL, all samples will be black. Default is NULL.

#' @param shapeby String specifying the metadata variable name for attributing shapes to samples. If NULL, all samples will be a round dot (pch = 19). Default is NULL. If there are more than 27 classes within the variable, samples will be attributed letters (A-Z, then a-z) automatically.

#' @param fillby String specifying the metadata variable with which to colour/fill in the boxes with. If NULL, all boxes will be filled in white. Default is NULL.

#' @param pairby .


#' @param connectby String specifying the metadata variable name for drawing a line connecting samples belonging to the same class. If NULL, samples are not connected. Default is NULL.


#' @param facetby .


#' @param wrap_facet .


#' @param overlay_boxplot Requires a logical value. If FALSE, will not overlay the boxplots ontop of one-another. If TRUE, the boxplots will all be plotted, one on top of another. Default is FALSE.

#' @param applyfilters Optional string specifying filtration setting "combos", used as a shorthand for setting the featcutoff, GenomeCompletenessCutoff, minl2fc and minabscorrcoeff arguments in JAMS plotting functions. If NULL, none of these arguments are set if not specified. Permissible values for applyfilters are "light", "moderate" or "stringent". The actual values vary whether the SummarizedExperiment object is taxonomical (LKT) or not. For a taxonomical SummarizedExperiment object, using "light" will set featcutoff=c(50, 5), GenomeCompletenessCutoff=c(5, 5), minl2fc=1, minabscorrcoeff=0.4; using "moderate" will set featcutoff=c(250, 15), GenomeCompletenessCutoff=c(10, 5), minl2fc=1, minabscorrcoeff=0.6; and using "stringent" will set featcutoff=c(2000, 15), GenomeCompletenessCutoff=c(30, 10), minl2fc=2, minabscorrcoeff=0.8. For non-taxonomical (i.e. functional) SummarizedExperiment objects, using "light" will set featcutoff=c(0, 0), minl2fc=1, minabscorrcoeff=0.4; using "moderate" will set featcutoff=c(5, 5), minl2fc=1, minabscorrcoeff=0.6; and using "stringent" will set featcutoff=c(50, 15), minl2fc=2.5, minabscorrcoeff=0.8. When using applyfilters, one can still set one or more of featcutoff, GenomeCompletenessCutoff, minl2fc and minabscorrcoeff, which will then take the user set value in lieu of those set by the applyfilters shorthand. Default is light.

#' @param featcutoff Requires a numeric vector of length 2 for specifying how to filter out features by relative abundance. The first value of the vector specifies the minimum relative abundance in Parts per Million (PPM) and the second value is the percentage of samples which must have at least that relative abundance. Thus, passing c(250, 10) to featcutoff would filter out any feature which does not have at least 250 PPM (= 0.025 percent) of relative abundance in at least 10 percent of all samples being plot. Please note that when using the subsetby option (q.v.) to automatically plot multiple plots of sample subsets, the featcutoff parameters are applied within the subset. The default is c(0, 0), meaning no feature is filtered. If NULL is passed, then the value defaults to c(0, 0). See also applyfilters for a shorthand way of applying multiple filtration settings.

#' @param GenomeCompletenessCutoff Requires a numeric vector of length 2 for specifying how to filter out features by genome completeness. This is, of course, only applicble for taxonomic shotgun SummarizedExperiment objects. When passed to non-taxonomic shotgun SummarizedExperiment objects, GenomeCompletenessCutoff will be ignored. The first value of the vector specifies the minimum genome completeness in percentage  and the second value is the percentage of samples which must have at least that genome completeness. Thus, passing c(50, 5) to GenomeCompletenessCutoff would filter out any taxonomic feature which does not have at least 50 percent of genome completeness in at least 5 percent of all samples being plot. Please note that when using the subsetby option (q.v.) to automatically plot multiple plots of sample subsets, the GenomeCompletenessCutoff parameters are applied within the subset. The default is c(0, 0), meaning no feature is filtered. If NULL is passed, then the value defaults to c(0, 0). See also applyfilters for a shorthand way of applying multiple filtration settings.

#' @param PPM_normalize_to_bases_sequenced Requires a logical value. Non-filtered JAMS feature counts tables (the BaseCounts assay within SummarizedExperiment objects) always includes unclassified taxonomical features (for taxonomical SummarizedExperiment objects) or unknown/unattributed functional features (for non-taxonomical SummarizedExperiment objects), so the relative abundance for each feature (see normalization) will be calculated in Parts per Million (PPM) by dividing the number of bases covering each feature by the sum of each sample column **previous to any filtration**. Relative abundances are thus representative of the entirety of the genomic content for taxonomical objects, whereas for non-taxonomical objects, strictly speaking, it is the abundance of each feature relative to only the coding regions present in the metagenome, even if these are annotationally unatributed. In other words, intergenic regions are not taken into account. In order to relative-abundance-normalize a **non-taxonomical** SummarizedExperiment object with the total genomic sequencing content, including non-coding regions, set PPM_normalize_to_bases_sequenced = TRUE. Default is FALSE.

#' @param addtit Optional string with text to append to heatmap main title. Default is NULL.

#' @param signiflabel String specifying the label used to determine if comparisions are significant or not. Deault is p.format.

#' @param max_pairwise_cats Numerical value specifying the maximum number of categories to be plot. Default is 4. This means that if you give a category that would require more than four different boxes, it will not be plot.

#' @param ignoreunclassified Requires a logical value. If set to TRUE, for taxonomical SummarizedExperiment objects, the feature "LKT__Unclassified" will be omitted from being shown. In the case of non-taxonomical SummarizedExperiment objects, the completely unannotated features will be omitted. For example, for an ECNumber SummarizedExperiment object, genes *without* an Enzyme Commission Number annotation (feature "EC_none") will not be shown. Statistics are, however, computed taking the completely unclassifed feature into account, so p-values will not change.

#' @param class_to_ignore String or vector specifying any classes which should lead to samples being excluded from the comparison within the variable passed to compareby. Default is N_A. This means that within any metadata variable passed to compareby containing the "N_A" string within that specific variable, the sample will be dropped from that comparison.

#' @export

plot_alpha_diversity <- function(ExpObj = NULL, measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson", "GeneCount"), stratify_by_kingdoms = TRUE, glomby = NULL, samplesToKeep = NULL, featuresToKeep = NULL, subsetby = NULL, compareby = NULL, compareby_order = NULL,colourby = NULL, shapeby = NULL, fillby = NULL, pairby = NULL, connectby = NULL, facetby = NULL, wrap_facet = FALSE, overlay_boxplot = FALSE, applyfilters = "light", featcutoff = NULL, GenomeCompletenessCutoff = NULL, PPM_normalize_to_bases_sequenced = FALSE, cdict = NULL, addtit = NULL, signiflabel = "p.format", max_pairwise_cats = 4, ignoreunclassified = TRUE, class_to_ignore = "N_A", returnstats = FALSE, ...){

    variables_to_fix <- c(compareby, subsetby, colourby, shapeby)

    #Hardwire PctFromCtgscutoff, as this should never be used without filtering because of huge amounts of false positives when evaluating taxonomic information from unassembled reads. The use of classifying unassembled reads is deprecated in JAMS and the default is to NOT classify unassembled reads, so this is usually not an issue.
    PctFromCtgscutoff <- c(70, 50)

    #Vet experiment object
    obj <- ExpObjVetting(ExpObj = ExpObj, samplesToKeep = samplesToKeep, featuresToKeep = NULL, glomby = glomby, variables_to_fix = variables_to_fix, class_to_ignore = class_to_ignore)

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

    #Initialize Graphics list
    gvec <- list()
    plotcount <- 1
    svec <- list()
    svn <- 1

    #subset by metadata column
    for (sp in 1:length(subset_points)){

        if (!(is.null(subsetby))){
            samplesToKeep <- rownames(colData(obj))[which(colData(obj)[ , subsetby] == subset_points[sp])]
            flog.info(paste("Plotting within", subset_points[sp]))
            subsetname <- subset_points[sp]
        } else {
            samplesToKeep <- rownames(colData(obj))
            subsetname <- "no_sub"
        }

        #See if there are enough samples and features to go ahead
        proceed <- TRUE
        curr_pt <- colData(obj)[samplesToKeep, ]

        if ((dim(curr_pt)[1] * dim(curr_pt)[2]) < 2){
            #There are less than 2 cells, a plot is meaningless.
            proceed <- FALSE
        }

        if (proceed){

            hmtypemsg <- "Alpha Diversity Plot"
            asPA <- FALSE
            hmasPA <- FALSE
            if (can_be_made_numeric(curr_pt[ , compareby])){
                stattype <- "spearman"
            } else {
                stattype <- "auto"
            }

            currobj <- filter_experiment(ExpObj = obj, featcutoff = presetlist$featcutoff, samplesToKeep = samplesToKeep, featuresToKeep = NULL, asPPM = TRUE, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, GenomeCompletenessCutoff = presetlist$GenomeCompletenessCutoff, PctFromCtgscutoff = presetlist$PctFromCtgscutoff)

        } else {

            flog.info("Unable to make plots with the current metadata for this comparison.")
            return(NULL)

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
            if (nrow(countmat) < 1){
                #abort, nothing is left over
                flog.info("None of the wanted features were found in SummarizedExperiment object when using the current filtration parameters.")
                return(NULL)
            }
        }

        if ("GeneCounts" %in% names(assays(currobj))){
            genecountmat <- as.matrix(assays(currobj)$GeneCounts)
            genecountmat <- genecountmat[rownames(countmat), ]
        } else {
            genecountmat <- NULL
        }

        if (all(c((analysis == "LKT"), (stratify_by_kingdoms == TRUE)))) {
            tt <- as.data.frame(rowData(currobj))
            tt <- tt[rownames(countmat), ]
            #Need at least 10 features within each Kingdom to measure alpha diversity
            table(tt$Kingdom)[which(table(tt$Kingdom) >= 10)]
            validanalyses <- names(table(tt$Kingdom)[which(table(tt$Kingdom) >= 10)])
        } else {
            validanalyses <- analysisname
        }

        #Loop though each valid analysis subsetting the current object to contain only features within that subset
        for (curranalysis in validanalyses){
            #Compose an appropriate title for the plot
            if (length(unique(subset_points)) > 1){
                maintit <- paste(hmtypemsg, curranalysis, paste("within", subset_points[sp]), sep = " | ")
                tablename <- paste("AlphaDiv", curranalysis, paste("within", subset_points[sp]), sep = "_")
            } else {
                maintit <- paste(hmtypemsg, curranalysis, sep = " | ")
                tablename <- paste("AlphaDiv", curranalysis, subsetname, sep = "_")
            }
            if (!is.null(addtit)) {
                maintit <- paste(addtit, maintit, sep = "\n")
            }

            #Define which features are going to be kept within the subset
            if (all(c((analysis == "LKT"), (stratify_by_kingdoms == TRUE)))){
                validfeats <- subset(tt, Kingdom == curranalysis)[ , analysisname]
                currcountmat <- countmat[validfeats, ]
            } else {
                currcountmat <- countmat
            }

            #calculate alpha diversity measures
            tmat <- t(currcountmat)
            alphadiv <- estimateR(tmat)
            for (meas in c("invsimpson", "simpson", "shannon")){
                alphadiv2 <- diversity(tmat, index = meas)
                alphadiv <- rbind(alphadiv, alphadiv2[colnames(alphadiv)])
                rownames(alphadiv)[nrow(alphadiv)] <- meas
            }
            alphadiv[c(3,5,6,7,8), ] <- round(alphadiv[c(3,5,6,7,8), ], 2)
            alphadiv[c(1,2,4), ] <- round(alphadiv[c(1,2,4), ], 0)
            alphadiv <- alphadiv[!(rownames(alphadiv) %in% c("se.chao1", "se.ACE")), ]
            rownames(alphadiv) <- unlist(sapply(1:nrow(alphadiv), function(x) { switch(rownames(alphadiv)[x], "S.obs" = "Observed", "S.chao1" = "Chao1", "S.ACE" = "ACE", "shannon" = "Shannon", "simpson" = "Simpson", "invsimpson" = "InvSimpson") }))

            if ("GeneCounts" %in% names(assays(currobj))){
                currgenecountmat <- genecountmat[rownames(currcountmat), ]
                genecountmatsum <- colSums(genecountmat)
                alphadiv <- rbind(alphadiv, genecountmatsum[colnames(alphadiv)])
                rownames(alphadiv)[nrow(alphadiv)] <- "GeneCounts"
            }

            alphadiv <- as.data.frame(t(alphadiv))
            alphadiv$Sample <- rownames(alphadiv)
            dat <- left_join(alphadiv, as.data.frame(curr_pt), by = "Sample")

            #Bank dat
            statsdf <- as.data.frame(dat)
            statsdf <- statsdf[ , which(colnames(statsdf) != "ACE")]
            statsdf$Sample_Set <- tablename
            measurescols <- colnames(statsdf)[!(colnames(statsdf) %in% c("Samples", "Sample_Set", colnames(curr_pt)))]
            othercols <- colnames(statsdf)[!(colnames(statsdf) %in% c("Samples", "Sample_Set", measurescols))]
            statsdf <- statsdf[ , c("Sample_Set", "Sample", measurescols, othercols) ]
            #if (!is.null(presetlist$filtermsg)){
            #    statsdf$Filtering <- as.character(presetlist$filtermsg)
            #}
            svec[[svn]] <- statsdf
            names(svec)[svn] <- tablename
            svn <- svn + 1

            dat$Compareby <- dat[ , compareby]

            #if there is an explicit order to compareby then set it to that
            if (!is.null(compareby_order)){
                dat$Compareby <- factor(dat$Compareby, levels = compareby_order)
            }

            if (!is.null(shapeby)){
                dat$Shape <- dat[ , shapeby]
            }

            if (!is.null(colourby)){
                dat$Colour <- dat[ , colourby]
            }

            if (!is.null(fillby)){
                dat$Fill <- dat[ , fillby]
            }

            if (!is.null(pairby)){
              dat$Pairby <- dat[ , pairby]
            }

            if (!is.null(connectby)){
                dat$Connect <- dat[ , connectby]
            }

            if (!is.null(facetby)){
                dat$Facetby <- dat[ , facetby]
            }

            measures_to_show <- measures[measures %in% colnames(alphadiv)]

            for (meas in measures_to_show){

                #Cheap trick, but it works
                currdat <- dat
                if (!is.null(pairby)) {
                  dat <- dat[order(dat$Pairby),]
                }
                colnames(currdat)[which(colnames(dat) == meas)] <- "AlphaMeas"
                #Start building a plot
                p <- ggplot(currdat, aes(x = Compareby, y = AlphaMeas))

                #If compareby is discrete, do boxplot, if numeric, do scatterplot
                if (is.numeric(currdat$Compareby)){
                    #Make a scatterplot
                    p <- p + geom_point()
                    p <- p + geom_smooth(method = lm, aes(group = 1), se = FALSE)
                    if (!(is.null(shapeby))){
                        p <- p + aes(shape = Shape)
                        numshapes <- length(unique(currdat$Shape))
                        p <- p + scale_shape_manual(values = 15:(numshapes + 15))
                    }
                    rotang <- 0
                    comparison <- cor.test(currdat$AlphaMeas, currdat$Compareby, method = "spearman")
                    correlstat <- paste0("spearman_corr_coeff=", round(comparison$estimate, 3))
                    overallp <- paste0("pval=", round(comparison$p.value, 4))
                    stattit <- paste(correlstat, overallp, sep = " | ")

                } else {

                    #Code for a boxplot
                    if(length(unique(currdat$Compareby)) < nrow(curr_pt)){
                        jitfact <- -( 0.3 / nrow(colData(currobj))) * (length(unique(currdat$Compareby))) + 0.25
                    } else {
                        jitfact <- 0
                    }

                    if (!overlay_boxplot){
                        p <- p + geom_boxplot(outlier.shape = NA)
                    }

                    if (!is.null(fillby)){
                        p <- p + aes(fill = Fill)
                        #if there is a colour dictionary, then use that
                        if (!(is.null(cdict))){
                            ct <- cdict[[fillby]]
                            groupcols <- setNames(as.character(ct$Colour), as.character(ct$Name))
                            p <- p + scale_fill_manual(values = groupcols)
                        }
                    }

                    if (!(is.null(shapeby))){
                        p <- p + geom_jitter(position = position_jitter(width = jitfact, height = 0.0), aes(shape = Shape))
                    } else {
                        p <- p + geom_jitter(position = position_jitter(width = jitfact, height = 0.0))
                    }

                    if (!(is.null(shapeby))){
                        p <- add_shape_to_plot_safely(p = p, shapevec = currdat$Shape, shapeby = shapeby, cdict = cdict)
                    }

                    if (overlay_boxplot){
                        p <- p + geom_boxplot(outlier.shape = NA)
                    }

                    rotang <- 0

                    if ((length(unique(currdat$Compareby)) > 1) && (length(unique(currdat$Compareby)) <= max_pairwise_cats)){
                        if (is.null(signiflabel)){
                            signiflabel <- "p.format"
                        }
                        #Add pval
                        my_comparisons <- combn(unique(as.character(currdat$Compareby)), m = 2, simplify = FALSE)

                        p <- p + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, paired = !is.null(pairby), label = signiflabel)
                    } else {
                        flog.warn("There are too many combinations to plot significance.")
                    }
                    rotang <- 90
                    stattit <- NULL
                }

                #Add bells and whistles
                if (!is.null(colourby)){
                    p <- p + aes(colour = Colour)

                    if (is.numeric(currdat$Colour)){
                        #If it is numeric, check that range is enough for a gradient
                        #if ((range(dat$colours)[2] - range(dat$colours)[1]) != 0){
                        #p <- p + scale_color_gradient(low = "blue", high = "red")
                        #} else {
                        #dat$colours <- as.character(dat$colours)
                        #groupcols <- setNames("black", unique(dat$colours))
                        #p <- p + scale_color_manual(values = groupcols)
                        #}
                        p <- p + scale_color_gradient(low = "blue", high = "red")
                    } else {
                        #if there is a colour dictionary, then use that
                        if (!(is.null(cdict))){
                            ct <- cdict[[colourby]]
                            groupcols <- setNames(as.character(ct$Colour), as.character(ct$Name))
                            p <- p + scale_color_manual(values = groupcols)
                        }
                    }
                }

                #Deal with titles and legends
                if (!is.null(facetby)){
                    if (wrap_facet){
                        p <- p + facet_wrap( ~ Facetby)
                    } else {
                        p <- p + facet_grid( ~ Facetby)
                    }
                }

                p <- p + theme_minimal()
                #Build plot title
                measureexpl <- switch(meas, "Observed" = "Number of Observed Features", "Chao1" = "Chao1 index", "ACE" = "Abundance Based Coverage Estimator", "Shannon" = "Shannon Index", "Simpson" = "Simpson Index", "InvSimpson" = "Inverse Simpson Index", "GeneCounts" = "Number of Individual Genes")

                msgs <- c(maintit, presetlist$filtermsg, measureexpl, stattit)
                plotit <- paste0(msgs, collapse = "\n")

                p <- p + ggtitle(plotit)

                if (!(is.null(colourby))){
                    p <- p + labs(colour = colourby)
                }

                if (!(is.null(shapeby))){
                    p <- p + labs(shape = shapeby)
                }

                p <- p + labs(x = compareby, y = measureexpl)
                p <- p + theme(axis.text.x = element_text(angle = rotang, size = rel(1), colour = "black"))
                p <- p + theme(plot.title = element_text(size = 10))

                gvec[[plotcount]] <- p
                names(gvec)[plotcount] <- paste(maintit, meas, sep = " | ")
                plotcount <- plotcount + 1

            }#End loop for plotting each measure
        }#End loop for each valid analysis
    }#End loop for each subset

    if (returnstats){
        gvec <- merge.list(gvec, svec)
    }

    return(gvec)

}
