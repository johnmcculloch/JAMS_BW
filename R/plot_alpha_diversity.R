#' plot_alpha_diversity(ExpObj = NULL, measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson", "GeneCount"), stratify_by_kingdoms = TRUE, glomby = NULL, samplesToKeep = NULL, featuresToKeep = NULL, subsetby = NULL, compareby = NULL, compareby_order = NULL, colourby = NULL, shapeby = NULL, fillby = NULL, facetby = NULL, wrap_facet = FALSE, overlay_boxplot = FALSE, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL, PPM_normalize_to_bases_sequenced = FALSE, cdict = NULL, addtit = NULL, signiflabel = "p.format", max_pairwise_cats = 4, ignoreunclassified = TRUE, class_to_ignore = "N_A", returnstats = FALSE, ...)
#'
#' Generates relative abundance plots per feature annotated by the metadata using as input a SummarizedExperiment object
#' @export

plot_alpha_diversity <- function(ExpObj = NULL, measures = c("Observed", "Chao1", "Shannon", "Simpson", "InvSimpson", "GeneCount"), stratify_by_kingdoms = TRUE, glomby = NULL, samplesToKeep = NULL, featuresToKeep = NULL, subsetby = NULL, compareby = NULL, compareby_order = NULL, colourby = NULL, shapeby = NULL, fillby = NULL, connectby = NULL, facetby = NULL, wrap_facet = FALSE, overlay_boxplot = FALSE, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL, PPM_normalize_to_bases_sequenced = FALSE, cdict = NULL, addtit = NULL, signiflabel = "p.format", max_pairwise_cats = 4, ignoreunclassified = TRUE, class_to_ignore = "N_A", returnstats = FALSE, ...){

    variables_to_fix <- c(compareby, subsetby, colourby, shapeby)

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

                        p <- p + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, paired = FALSE, label = signiflabel)
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
