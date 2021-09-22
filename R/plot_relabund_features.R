#' plot_relabund_features(ExpObj = NULL, glomby = NULL, samplesToKeep = NULL, featuresToKeep = NULL, aggregatefeatures = FALSE, aggregatefeatures_label = "Sum_of_wanted_features", subsetby = NULL, compareby = NULL, compareby_order = NULL, colourby = NULL, shapeby = NULL, fillby = NULL, facetby = NULL, wrap_facet = FALSE, overlay_boxplot = FALSE, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL, ntop = NULL, minabscorrcoeff = NULL, adjustpval = TRUE, padjmeth = "fdr", showonlypbelow = NULL, showonlypadjusted = FALSE, maxl2fc = NULL, minl2fc = NULL, addtit = NULL, PPM_normalize_to_bases_sequenced = FALSE, uselog = FALSE, statsonlog = FALSE, cdict = NULL, stratify_by_taxlevel = NULL, maxnumplots = NULL, signiflabel = "p.format", max_pairwise_cats = 4, dump_interpro_descriptions_to_plot = FALSE, numthreads = 1, nperm = 99, ignoreunclassified = TRUE, class_to_ignore = "N_A", maxnumtaxa = 20, horizontal = TRUE, plot_points_on_taxonomy = FALSE, use_heatmap_for_stratification = TRUE, return_plots = FALSE, rescale_axis_quantiles = NULL, ...)
#'
#' Generates relative abundance plots per feature annotated by the metadata using as input a SummarizedExperiment object
#' @export

plot_relabund_features <- function(ExpObj = NULL, glomby = NULL, samplesToKeep = NULL, featuresToKeep = NULL, aggregatefeatures = FALSE, aggregatefeatures_label = "Sum_of_wanted_features", subsetby = NULL, compareby = NULL, compareby_order = NULL, colourby = NULL, shapeby = NULL, fillby = NULL, facetby = NULL, wrap_facet = FALSE, overlay_boxplot = FALSE, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL, ntop = NULL, minabscorrcoeff = NULL, adjustpval = TRUE, padjmeth = "fdr", showonlypbelow = NULL, showonlypadjusted = FALSE, maxl2fc = NULL, minl2fc = NULL, addtit = NULL, PPM_normalize_to_bases_sequenced = FALSE, uselog = FALSE, statsonlog = FALSE, cdict = NULL, stratify_by_taxlevel = NULL, maxnumplots = NULL, signiflabel = "p.format", max_pairwise_cats = 4, dump_interpro_descriptions_to_plot = FALSE, numthreads = 1, nperm = 99, ignoreunclassified = TRUE, class_to_ignore = "N_A", maxnumtaxa = 20, horizontal = TRUE, plot_points_on_taxonomy = FALSE, use_heatmap_for_stratification = TRUE, return_plots = FALSE, rescale_axis_quantiles = NULL, ...){

    variables_to_fix <- c(compareby, subsetby, colourby, shapeby)

    #Vet experiment object
    obj <- ExpObjVetting(ExpObj = ExpObj, samplesToKeep = samplesToKeep, featuresToKeep = NULL, glomby = glomby, variables_to_fix = variables_to_fix, class_to_ignore = class_to_ignore)

    analysis <- metadata(obj)$analysis
    if (!is.null(glomby)){
        analysisname <- glomby
    } else {
        analysisname <- analysis
    }

    presetlist <- declare_filtering_presets(analysis = analysis, applyfilters = applyfilters, featcutoff = featcutoff, GenomeCompletenessCutoff = GenomeCompletenessCutoff, PctFromCtgscutoff = PctFromCtgscutoff, maxl2fc = maxl2fc, minl2fc = minl2fc, minabscorrcoeff = minabscorrcoeff)

    if (!(is.null(subsetby))){
        subset_points <- sort(unique(colData(obj)[, which(colnames(colData(obj)) == subsetby)]))
    } else {
        subset_points <- "none"
    }

    #Initialize Graphics list
    gvec <- list()
    plotcount <- 1

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

            hmtypemsg <- "Relative Abundance Plot"
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

        #There must be at least one feature requested in the object

        if (is.null(featuresToKeep)){

            wantedfeatures <- rownames(currobj)

        } else {

            #Just be sure
            featuresToKeep <- unique(featuresToKeep)
            wantedfeatures <- featuresToKeep[featuresToKeep %in% rownames(currobj)]

            if(length(wantedfeatures) < 1){
                #abort, nothing is left over
                flog.info("None of the wanted features were found in SummarizedExperiment object when using the current filtration parameters.")
                return(NULL)
            }

            if (length(wantedfeatures) < length(featuresToKeep)){
                #warn that some features were not found
                flog.warn(paste("Some of the wanted features were not found in SummarizedExperiment object when using the current filtration parameters. Only", paste0(length(wantedfeatures), "/", length(featuresToKeep)), "are still present."))
            }

        }

        if (length(wantedfeatures) > 0){

            #Compose an appropriate title for the plot
            if (length(unique(subset_points)) > 1){
                maintit <- paste(hmtypemsg, analysisname, paste("within", subset_points[sp]), sep = " | ")
            } else {
                maintit <- paste(hmtypemsg, analysisname, sep = " | ")
            }
            if (!is.null(addtit)) {
                maintit <- paste(addtit, maintit, sep = "\n")
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

            if ("GenomeCompleteness" %in% names(assays(currobj))){
                genomecompletenessdf <- as.matrix(assays(currobj)$GenomeCompleteness)
            } else {
                genomecompletenessdf <- NULL
            }

            if ("PctFromCtgs" %in% names(assays(currobj))){
                PctFromCtgsdf <- as.matrix(assays(currobj)$PctFromCtgs)
            } else {
                PctFromCtgsdf <- NULL
            }

            #Aggregate if appropriate
            if (aggregatefeatures == TRUE){
                #If aggregating features, then cull to wantedfeatures now and aggregate
                wantedfeatures <- wantedfeatures[wantedfeatures %in% rownames(countmat)]
                countmat <- countmat[wantedfeatures, ]
                originalwantedfeatures <- wantedfeatures
                #Count matrix should not be in log2 at this point
                aggcountmat <- colSums(countmat)
                aggcountmat <- t(as.matrix(aggcountmat))
                rownames(aggcountmat) <- aggregatefeatures_label
                countmat <- aggcountmat

                #Also aggregate Genome Completeness, if applicable
                if ("GenomeCompleteness" %in% names(assays(currobj))){
                    genomecompletenessdf <- genomecompletenessdf[wantedfeatures, ]
                    agggenomecompletenessdf <- colSums(genomecompletenessdf)
                    agggenomecompletenessdf <- t(as.matrix(agggenomecompletenessdf))
                    rownames(agggenomecompletenessdf) <- aggregatefeatures_label
                    genomecompletenessdf <- agggenomecompletenessdf
                }

                #Also get mean percentage from contigs, if applicable
                if ("PctFromCtgs" %in% names(assays(currobj))){
                    PctFromCtgsdf <- PctFromCtgsdf[wantedfeatures, ]
                    aggPctFromCtgsdf <- colMeans(PctFromCtgsdf)
                    aggPctFromCtgsdf <- t(as.matrix(aggPctFromCtgsdf))
                    rownames(aggPctFromCtgsdf) <- aggregatefeatures_label
                    PctFromCtgsdf <- aggPctFromCtgsdf
                }

                wantedfeatures <- aggregatefeatures_label
            }

            matrixSamples <- colnames(countmat)
            matrixRows <- rownames(countmat)

            #Calculate matrix stats and get new matrix.
            cl <- colData(currobj)[ , which(colnames(colData(currobj)) == compareby)]
            discretenames <- sort(unique(cl))

            matstats <- calculate_matrix_stats(countmatrix = countmat, uselog = FALSE, statsonlog = FALSE, stattype = stattype, classesvector = cl, invertbinaryorder = FALSE, numthreads = numthreads, nperm = nperm)

            ffeatmsg <- paste0("Number of features assessed = ", nrow(matstats))

            #Cull to only features wanted
            wantedfeatures <- wantedfeatures[wantedfeatures %in% rownames(matstats)]
            matstats <- matstats[wantedfeatures, ]
            #Reorder matrix by p-value
            matstats <- matstats[order(matstats$pval), ]
            topcats <- nrow(matstats)
            if (!(is.null(ntop))) {
                topcats <- min(topcats, ntop)
            }

            if (!is.null(showonlypbelow)){
                if (showonlypadjusted == TRUE) {
                    sigmeas <- paste("padj", padjmeth, sep = "_")
                } else {
                    sigmeas <- "pval"
                }
                rowcutoff <- which(matstats[ , sigmeas] < showonlypbelow)
            } else {
                rowcutoff <- 1:nrow(matstats)
            }

            #Limit number of features to requested number or number available
            rowcutoff <- rowcutoff[1:(min(topcats, length(rowcutoff)))]

            if (any(c(is.na(rowcutoff), (length(rowcutoff) == 0)))){
                #abort, nothing is left over
                flog.warn("None of the wanted features were found in SummarizedExperiment object when using the current p-value filtration parameters.")

                return(NULL)
            }

            matstats <- matstats[rowcutoff, ]

            #Filter by l2fc if applicable
            if (all(c((!is.null(presetlist$minl2fc)), ("absl2fc" %in% colnames(matstats))))){

                matstats <- subset(matstats, absl2fc >= presetlist$minl2fc)

                if (nrow(matstats) < 1){
                    #abort, nothing is left over
                    flog.warn("None of the wanted features were found in the SummarizedExperiment object when using the current log2 foldchange filtration parameters.")

                    return(NULL)
                }
            }

            #Filter by correlation coefficient, if applicable
            if (all(c((!is.null(presetlist$minabscorrcoeff)), ("abscorrel" %in% colnames(matstats))))){
                matstats <- subset(matstats, abscorrel >= presetlist$minabscorrcoeff)

                if (nrow(matstats) < 1){
                    #abort, nothing is left over
                    flog.warn("None of the wanted features were found in the SummarizedExperiment object when using the current absolute correlation coefficient filtration parameters.")

                    return(NULL)
                }
            }

            #Redefine countmat to include only features matching filtering criteria
            countmat <- countmat[rownames(matstats), ]

            if (class(countmat)[1] != "matrix"){
                countmat <- t(as.matrix(countmat))
                rownames(countmat) <- rownames(matstats)
            }

            if ("GenomeCompleteness" %in% names(assays(currobj))){
                #genomecompletenessdf <- as.matrix(assays(currobj)$GenomeCompleteness)
                genomecompletenessdf <- genomecompletenessdf[rownames(matstats), ]
                if (class(genomecompletenessdf)[1] != "matrix"){
                    genomecompletenessdf <- t(as.matrix(genomecompletenessdf))
                    rownames(genomecompletenessdf) <- rownames(matstats)
                }
            } else {
                genomecompletenessdf <- NULL
            }

            if ("PctFromCtgs" %in% names(assays(currobj))){
                #PctFromCtgsdf <- as.matrix(assays(currobj)$PctFromCtgs)
                PctFromCtgsdf <- PctFromCtgsdf[rownames(matstats), ]
                if (class(PctFromCtgsdf)[1] != "matrix"){
                    PctFromCtgsdf <- t(as.matrix(PctFromCtgsdf))
                    rownames(PctFromCtgsdf) <- rownames(matstats)
                }
            } else {
                PctFromCtgsdf <- NULL
            }

        } else {

            #abort, nothing is left over
            flog.warn("None of the wanted features were found in SummarizedExperiment object when using the current filtration parameters.")

            return(NULL)
        }

        #Now, for the tricky bit of stratifying by taxa
        if (!is.null(stratify_by_taxlevel)){
            if (stratify_by_taxlevel == TRUE){
                stratify_by_taxlevel <- "LKT"
            }

            if (stratify_by_taxlevel == FALSE){
                stratify_by_taxlevel <- NULL
            }

            if (aggregatefeatures == TRUE){
                featnamesforsubset <- originalwantedfeatures
            } else {
                featnamesforsubset <- rownames(countmat)
            }

            #See if current SummarizedExperiment object allows for stratification by taxa.
            if (all(c("allfeaturesbytaxa_index", "allfeaturesbytaxa_matrix") %in% names(metadata(currobj)))){
                taxsplit <- retrieve_features_by_taxa(FuncExpObj = currobj, wantedfeatures = featnamesforsubset, wantedsamples = colnames(countmat), asPPM = TRUE, PPMthreshold = 0)
                colnames(taxsplit)[which(colnames(taxsplit) == compareby)] <- "Compareby"
            } else {
                flog.warn("Current SummarizedExperiment object does not contain the necessary data for stratifying this function by taxonomy. Check your input.")
            }

            LKTcolumns <- colnames(taxsplit)[!(colnames(taxsplit) %in% unique(c(colnames(curr_pt), c("Sample", "Accession", "Compareby"))))]

            #Fix taxsplit, if aggregating
            if (aggregatefeatures == TRUE){
                aggtaxsplit <- NULL
                for (smpl in unique(taxsplit$Sample)){
                    taxsplitsmpl <- subset(taxsplit, Sample == smpl)
                    taxvalues <- as.matrix(taxsplitsmpl[ , LKTcolumns])
                    aggtaxvalues <- colSums(taxvalues)
                    aggtaxvalues <- t(as.matrix(aggtaxvalues))
                    rownames(aggtaxvalues) <- aggregatefeatures_label
                    aggtaxvalues <- as.data.frame(aggtaxvalues)
                    aggtaxsplitsmpl <- cbind((taxsplitsmpl[1, colnames(taxsplitsmpl)[!(colnames(taxsplitsmpl) %in% LKTcolumns)]]), aggtaxvalues)
                    aggtaxsplitsmpl[which(colnames(aggtaxsplitsmpl) == "Accession")] <- aggregatefeatures_label
                    aggtaxsplit <- rbind(aggtaxsplit, aggtaxsplitsmpl)
                }
                taxsplit <- aggtaxsplit
            }

            data(JAMStaxtable)
            data(Gram)
            tt <- JAMStaxtable[ , which(colnames(JAMStaxtable) %in% c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "LKT"))]
            tt <- tt[!(duplicated(tt$LKT)), ]
            LKTcolumns <- colnames(taxsplit)[!(colnames(taxsplit) %in% c("Sample", "Accession", colnames(curr_pt)))]
            tt <- subset(tt, LKT %in% LKTcolumns)
            tt <- tt[ , c(stratify_by_taxlevel, "Phylum")]
            Gram$Kingdom <- NULL
            tt <- left_join(tt, Gram, by = "Phylum")
            phycols <- setNames(as.character(tt$PhylumColour), as.character(tt$Phylum))[unique(tt$Phylum)]
            phycols <- c(phycols, Remainder = "#000000")
        }

        flog.info("Plotting results...")
        flog.info(paste("There are", nrow(countmat), "features to plot."))

        for (feat in rownames(countmat)){
            dat <- data.frame(Sample = names(countmat[feat, ]), PPM = as.numeric(countmat[feat, ]), Compareby = cl, stringsAsFactors = FALSE)
            rownames(dat) <- dat$Sample

            #if there is an explicit order to compareby then set it to that
            if (!is.null(compareby_order)){
                dat$Compareby <- factor(dat$Compareby, levels = compareby_order)
            }

            if (!is.null(shapeby)){
                dat$Shape <- curr_pt[rownames(dat) , shapeby]
            }

            if (!is.null(fillby)){
                dat$Fill <- curr_pt[rownames(dat) , fillby]
            }

            if (!is.null(colourby)){
                if (colourby == "GenomeCompleteness"){
                    pctgencompdf <- t(genomecompletenessdf) * 100
                    dat$Colour <- pctgencompdf[rownames(dat), feat]
                    #cap to 400%
                    dat$Colour[which(dat$Colour > 400)] <- 400
                } else if (colourby == "PctFromCtgs"){
                    pctctgsdf <- t(PctFromCtgsdf)
                    dat$Colour <- pctctgsdf[rownames(dat), feat]
                } else {
                    dat$Colour <- curr_pt[rownames(dat) , colourby]
                }
            }

            #Start building a plot
            p <- ggplot(dat, aes(x = Compareby, y = PPM))

            if (matstats$Method[1] %in% c("spearman", "pearson")){
                #Make a scatterplot
                p <- p + geom_point()
                p <- p + geom_smooth(method = lm, aes(group=1), se = FALSE)
                if (!(is.null(shapeby))){
                    p <- add_shape_to_plot_safely(p = p, shapevec = dat$Shape, shapeby = shapeby, cdict = cdict)
                }
                rotang <- 0

            } else {
                #Code for a boxplot
                if(length(discretenames) < nrow(curr_pt)){
                    jitfact <- -( 0.3 / nrow(colData(currobj))) * (length(discretenames)) + 0.25
                } else {
                    jitfact <- 0
                }

                if (!overlay_boxplot){
                    p <- p + geom_boxplot(outlier.shape = NA)
                }

                if (!(is.null(shapeby))){
                    p <- p + geom_jitter(position = position_jitter(width = jitfact, height = 0.0), aes(shape = Shape))
                    p <- add_shape_to_plot_safely(p = p, shapevec = dat$Shape, shapeby = shapeby, cdict = cdict)
                } else {
                    p <- p + geom_jitter(position = position_jitter(width = jitfact, height = 0.0))
                }

                if (overlay_boxplot){
                    p <- p + geom_boxplot(outlier.shape = NA)
                }

                if ((length(discretenames) > 1) && (length(discretenames) <= max_pairwise_cats)){
                    if (is.null(signiflabel)){
                        signiflabel <- "p.format"
                    }
                    #Add pval
                    my_comparisons <- combn(discretenames, m = 2, simplify = FALSE)
                    p <- p + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, paired = FALSE, label = signiflabel)
                } else {
                    flog.warn("There are too many combinations to plot significance.")
                }
                rotang <- 90
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

            if (!is.null(colourby)){
                p <- p + aes(colour = Colour)

                if (colourby == "GenomeCompleteness"){
                    p <- p + scale_fill_gradientn(aesthetics = "colour", colours = c("white", "forestgreen", "blue", "firebrick1", "black"),  values = scales::rescale(c(0, 100, 200, 300, 400), to = c(0, (400/max(dat$Colour)))))
                } else {
                    if (is.numeric(dat$Colour)){
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
            overallpmeth <- matstats[feat, "Method"]
            overallp <- paste0("pval=", round(matstats[feat, "pval"], 4))
            overalladjp <- paste0("padj_fdr=", round(matstats[feat, "padj_fdr"], 4))
            if("stat" %in% colnames(matstats)){
                overallstat <- paste0("stat=", round(matstats[feat, "stat"], 4))
            } else {
                overallstat <- NULL
            }
            stattit <- paste(overallpmeth, overallp, overalladjp, overallstat, ffeatmsg, sep = " | ")

            if ("correl" %in% colnames(matstats)){
                correlstat <- paste0("corr_coeff=", round(matstats[feat, "correl"], 3))
                stattit <- paste(stattit, correlstat, sep = "\n")
            }

            if ("l2fc" %in% colnames(matstats)){
                l2fcmsg <- paste0("Log2FC=", round(matstats[feat, "l2fc"], 3))
                l2fcmeaning <- paste("Positive l2fc means increased in", discretenames[1])
                l2fcmsg <- paste(l2fcmsg, l2fcmeaning, sep = " | ")
                stattit <- paste(stattit, l2fcmsg, sep = "\n")
            }

            #Add description to feature, if applicable
            if (analysis != "LKT"){
                featdesc <- rowData(currobj)[feat, "Description"]
                featname <- paste(feat, featdesc)
            } else {
                featname <- feat
            }

            msgs <- c(maintit, featname, stattit)
            plotit <- paste0(msgs, collapse = "\n")

            p <- p + ggtitle(plotit)

            if (!(is.null(colourby))){
                p <- p + labs(colour = colourby)
            }

            if (!(is.null(shapeby))){
                p <- p + labs(shape = shapeby)
            }

            if (uselog == TRUE){
                ytit <- "Relative Abundance in PPM"
                #p <- p + coord_trans(y = "log2", clip = "off")
                p <- p + scale_y_continuous(trans = scales::pseudo_log_trans(base = 2), breaks = scales::trans_breaks("log2", function(x) {((2 ^ x) - 1)}), labels = scales::trans_format("log2", function(x) {((2 ^ x) - 1)}))
            } else {
                ytit <- "Relative Abundance in PPM"
            }

            p <- p + labs(x = compareby, y = ytit)
            p <- p + theme(axis.text.x = element_text(angle = rotang, size = rel(1), colour = "black"))
            p <- p + theme(plot.title = element_text(size = 10))

            if (!return_plots){
                #print plot on the fly
                print(p)
            }

            gvec[[plotcount]] <- p
            names(gvec)[plotcount] <- paste(maintit, feat, sep = " | ")
            plotcount <- plotcount + 1

            if (all(c(dump_interpro_descriptions_to_plot, analysis == "Interpro"))){
                data(InterproDict)
                infotable <- as.data.frame(t(as.data.frame(InterproDict[feat, c("Abstract", "Citations")])))
                if (nchar(paste0(InterproDict[feat, c("Abstract", "Citations")], collapse = "")) > 2500){
                    fontsize <- 7
                } else {
                    fontsize <- 10
                }
                print_table(tb = infotable, tabletitle = paste0(InterproDict[feat, c("Accession", "Description")], collapse = " "), fontsize = fontsize, numrows = 20)
            }

            if (!is.null(stratify_by_taxlevel)){
                currtaxsplit <- subset(taxsplit, Accession == feat)

                #for (grp in unique(currtaxsplit$Compareby)){
                    p <- NULL
                    dat <- NULL
                    #currtaxsplitgrp <- subset(currtaxsplit, Compareby == grp)
                    currtaxsplitgrp <- currtaxsplit
                    LKTcolumns <- colnames(currtaxsplitgrp)[!(colnames(currtaxsplitgrp) %in% unique(c(colnames(curr_pt), c("Sample", "Accession", "Compareby"))))]

                    #Eliminate empties
                    LKTsToKeep <- names(which(colSums(currtaxsplitgrp[ , LKTcolumns]) > 0))
                    currtaxsplitgrp <- currtaxsplitgrp[ , c("Sample", LKTsToKeep)]
                    dat <- currtaxsplitgrp %>% gather(Taxon, PPM, 2:ncol(currtaxsplitgrp))

                    #Start building a plot
                    dat <- left_join(dat, as.data.frame(curr_pt), by = "Sample")
                    colnames(dat)[which(colnames(dat) == compareby)] <- "Compareby"
                    if (!(is.null(shapeby))){
                        colnames(dat)[which(colnames(dat) == shapeby)] <- "Shape"
                    }
                    if (!(is.null(colourby))){
                        colnames(dat)[which(colnames(dat) == colourby)] <- "Colour"
                    }

                    #Colour in by phylum
                    Taxon2phylum <- tt[which(tt[ , stratify_by_taxlevel] %in% dat$Taxon), ]
                    colnames(Taxon2phylum)[which(colnames(Taxon2phylum) == stratify_by_taxlevel)] <- "Taxon"
                    dat <- left_join(dat, Taxon2phylum, by = "Taxon")

                    #Order and aggregate if more than 30
                    tally <- aggregate(PPM ~ Taxon, data = dat, FUN = "sum")
                    tally <- tally[order(tally$PPM, decreasing = TRUE), ]

                    orddat <- NULL
                    for (Txn in tally$Taxon[1:min(maxnumtaxa, nrow(tally))]){
                        datsplit <- subset(dat, Taxon == Txn)
                        orddat <- rbind(orddat, datsplit)
                    }
                    #Aggregate if there are leftovers
                    if (nrow(tally) > maxnumtaxa){
                        TaxaToAgg <- tally[maxnumtaxa:nrow(tally), ]$Taxon
                        datremainder <- subset(dat, Taxon %in% TaxaToAgg)
                        remaindertally <- aggregate(PPM ~ Sample, data = datremainder, FUN = "sum")
                        remaindertally$Taxon <- "Other_Taxa"
                        datremainder$Taxon <- NULL
                        datremainder$PPM <- NULL
                        datremainder$Phylum <- NULL
                        datremainder$Gram <- NULL
                        datremainder <- datremainder[!(duplicated(datremainder)), ]
                        aggremainder <- left_join(remaindertally, datremainder, by = "Sample")
                        aggremainder$Phylum <- "Remainder"
                        aggremainder$Gram <- "Remainder"
                        orddat <- rbind(orddat, aggremainder)
                    }

                    dat <- orddat

                    dat$Taxon <- factor(dat$Taxon, levels = unique(dat$Taxon))

                    p <- ggplot(dat, aes(x = Taxon, y = PPM))

                    p <- p + geom_boxplot(aes(fill = Phylum), outlier.shape = NA)
                    #Rescale to exclude outliers
                    if (!is.null(rescale_axis_quantiles)){
                        p <- p + scale_y_continuous(limits = quantile(dat$PPM, rescale_axis_quantiles))
                    }
                    p <- p + scale_fill_manual(values = phycols[(names(phycols) %in% dat$Phylum)])

                    if ((length(unique(dat$Taxon))) < (nrow(dat))){
                        jitfact <- -( 0.3 / (nrow(dat))) * (length(unique(dat$Taxon))) + 0.25
                    } else {
                        jitfact <- 0
                    }


                    #p <- p + geom_violin()

                    if (plot_points_on_taxonomy == TRUE){
                        if (!(is.null(shapeby))){
                            p <- p + geom_jitter(position = position_jitter(width = jitfact, height = 0.0), aes(shape = Shape))
                            p <- add_shape_to_plot_safely(p = p, shapevec = dat$Shape, shapeby = shapeby, cdict = cdict)
                        } else {
                            p <- p + geom_jitter(position = position_jitter(width = jitfact, height = 0.0))
                        }
                    }


                    if (!is.null(colourby)){
                        p <- p + aes(colour = Colour)

                        if (colourby == "GenomeCompleteness"){
                            p <- p + scale_fill_gradientn(aesthetics = "colour", colours = c("white", "forestgreen", "blue", "firebrick1", "black"),  values = scales::rescale(c(0, 100, 200, 300, 400), to = c(0, (400/max(dat$Colour)))))
                        } else {
                            if (is.numeric(dat$Colour)){
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
                    }

                    if (horizontal == TRUE){
                        p <- p + coord_flip()
                    }

                    if (wrap_facet){
                        p <- p + facet_wrap( ~ Compareby)
                    } else {
                        p <- p + facet_grid( ~ Compareby)
                    }

                    p <- p + theme_minimal()
                    p <- p + theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size=1))
                    plotitstrat <- paste0(c(maintit, featname), collapse = "\n")

                    p <- p + ggtitle(plotitstrat)

                    if (!(is.null(colourby))){
                        p <- p + labs(colour = colourby)
                    }

                    if (!(is.null(shapeby))){
                        p <- p + labs(shape = shapeby)
                    }

                    if (uselog == TRUE){
                        ytit <- "Relative Abundance in PPM"
                        #p <- p + coord_trans(y = "log2", clip = "off")
                        p <- p + scale_y_continuous(trans = scales::pseudo_log_trans(base = 2), breaks = scales::trans_breaks("log2", function(x) {((2 ^ x) - 1)}), labels = scales::trans_format("log2", function(x) {((2 ^ x) - 1)}))

                    } else {
                        ytit <- "Relative Abundance in PPM"
                    }
                    p <- p + labs(x = "Contributing Taxon", y = ytit)

                    if (!horizontal){
                        p <- p + theme(axis.text.x = element_text(colour = "black", angle = rotang, size = rel(1)))
                    }
                    #df <- data.frame(x = factor(levels(dat$Taxon)), colour = factor(dat$Phcol[match(levels(dat$Taxon), dat$Taxon)]))
                    #p + geom_tile(data = df, aes(x = x, y = 2, fill = colour))
                    #p <- p + theme(axis.text.x = element_text(colour = phcol[dat$Phylum[!duplicated(dat$Phylum)]]))
                    p <- p + theme(plot.title = element_text(size = 10))

                    if (!return_plots){
                        #print plot on the fly
                        print(p)
                    }

                    gvec[[plotcount]] <- p
                    names(gvec)[plotcount] <- paste(maintit, feat, stratify_by_taxlevel, sep = " | ")
                    plotcount <- plotcount + 1
                #}#End loop for plotting stratify_by_taxlevel within each group
            }#End conditional of stratifying by taxonomy
        }#End loop for plotting each feature

    }#End loop for each subset

    if (return_plots){
        #Return plots, as nothing was printed
        return(gvec)
    }
}
