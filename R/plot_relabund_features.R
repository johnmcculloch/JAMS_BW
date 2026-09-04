#' plot_relabund_features(ExpObj = NULL, glomby = NULL, samplesToKeep = NULL, featuresToKeep = NULL, only_allow_CSBs = FALSE, aggregatefeatures = FALSE, aggregatefeatures_label = "Sum_of_wanted_features", subsetby = NULL, compareby = NULL, wilcox_paired_by = NULL, compareby_order = NULL, invertbinaryorder = FALSE, colourby = NULL, shapeby = NULL, fillby = NULL, connectby = NULL, facetby = NULL, wrap_facet = FALSE, overlay_boxplot = FALSE, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, ntop = NULL, minabscorrcoeff = NULL, adjustpval = TRUE, padjmeth = "fdr", showonlypbelow = NULL, showonlypadjusted = FALSE, maxl2fc = NULL, minl2fc = NULL, addtit = NULL, PPM_normalize_to_bases_sequenced = FALSE, log2tran_main_plot = FALSE, log2tran_strat_plot = FALSE, statsonlog = FALSE, y_axis_range = NULL, cdict = NULL, stratify_by_taxlevel = NULL, maxnumplots = NULL, signiflabel = "p.format", max_pairwise_cats = 4, dump_interpro_descriptions_to_plot = FALSE, numthreads = 1, nperm = 99, ignoreunclassified = TRUE, class_to_ignore = "N_A", maxnumtaxa = 20, horizontal = TRUE, plot_points_on_taxonomy = FALSE, use_cladogram_for_stratification = TRUE, return_taxon_stratification_df = FALSE, return_plots = FALSE, rescale_axis_quantiles = NULL, fun_for_l2fc = "geom_mean", ...)
#'
#' Generates relative abundance plots per feature annotated by the metadata using as input a SummarizedExperiment object
#'
#' When compareby = NULL, a single boxplot is drawn per feature showing the spread of relative abundance across ALL samples in the (sub)set, with no statistical comparison. If stratify_by_taxlevel is also set, the taxon-stratified boxplot(s) are likewise drawn across all samples as a single group. This is useful for inspecting which taxa contribute to a functional feature within a single group.
#' @export

plot_relabund_features <- function(ExpObj = NULL, glomby = NULL, samplesToKeep = NULL, featuresToKeep = NULL, only_allow_CSBs = FALSE, aggregatefeatures = FALSE, aggregatefeatures_label = "Sum_of_wanted_features", subsetby = NULL, compareby = NULL, wilcox_paired_by = NULL, compareby_order = NULL, invertbinaryorder = FALSE, colourby = NULL, shapeby = NULL, fillby = NULL, connectby = NULL, facetby = NULL, wrap_facet = FALSE, overlay_boxplot = FALSE, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, ntop = NULL, minabscorrcoeff = NULL, adjustpval = TRUE, padjmeth = "fdr", showonlypbelow = NULL, showonlypadjusted = FALSE, maxl2fc = NULL, minl2fc = NULL, addtit = NULL, PPM_normalize_to_bases_sequenced = FALSE, log2tran_main_plot = FALSE, log2tran_strat_plot = FALSE, statsonlog = FALSE, y_axis_range = NULL, cdict = NULL, stratify_by_taxlevel = NULL, maxnumplots = NULL, signiflabel = "p.format", max_pairwise_cats = 4, dump_interpro_descriptions_to_plot = FALSE, numthreads = 1, nperm = 99, ignoreunclassified = TRUE, class_to_ignore = "N_A", maxnumtaxa = 20, horizontal = TRUE, plot_points_on_taxonomy = FALSE, use_cladogram_for_stratification = TRUE, show_prevalence_in_cladogram = TRUE, return_taxon_stratification_df = FALSE, return_plots = FALSE, rescale_axis_quantiles = NULL, fun_for_l2fc = "geom_mean", ...){

    #Account for JAMS2 spaces
    taxonomic_spaces <- c("LKT", "Contig_LKT", "ConsolidatedGenomeBin", "MB2bin", "16S")

    #SINGLE_GROUP: Decide whether we are in single-group (no comparison) mode.
    single_group <- is.null(compareby)
    if (single_group){
        flog.info("compareby is NULL: plotting a single boxplot per feature across all samples with no statistical comparison.")
        #A dummy label used for the constant grouping on the x-axis.
        single_group_label <- "All_samples"
    }

    #SINGLE_GROUP: Do not pass a NULL compareby to variables_to_fix (it would be dropped anyway, but be explicit).
    variables_to_fix <- c(compareby, subsetby, colourby, shapeby)

    #Vet experiment object
    obj <- ExpObjVetting(ExpObj = ExpObj, samplesToKeep = samplesToKeep, featuresToKeep = NULL, only_allow_CSBs = only_allow_CSBs, glomby = glomby, variables_to_fix = variables_to_fix, class_to_ignore = class_to_ignore)

    analysis <- metadata(obj)$analysis
    if (!is.null(glomby)){
        analysisname <- glomby
    } else {
        analysisname <- analysis
    }

    presetlist <- declare_filtering_presets(analysis = analysis, applyfilters = applyfilters, featcutoff = featcutoff, GenomeCompletenessCutoff = GenomeCompletenessCutoff, maxl2fc = maxl2fc, minl2fc = minl2fc, minabscorrcoeff = minabscorrcoeff)

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
            #SINGLE_GROUP: Only interrogate compareby for stat type when there is a compareby.
            if (single_group){
                stattype <- "variance"
            } else if (can_be_made_numeric(curr_pt[ , compareby])){
                stattype <- "spearman"
            } else {
                stattype <- "auto"
            }

            currobj <- filter_experiment(SEobj = obj, featcutoff = presetlist$featcutoff, samplesToKeep = samplesToKeep, featuresToKeep = NULL, only_allow_CSBs = FALSE, normalization = "relabund", PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, GenomeCompletenessCutoff = presetlist$GenomeCompletenessCutoff)

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
            countmat <- as.matrix(assays(currobj)$PPM)

            if (ignoreunclassified == TRUE){
                dunno <- c(paste(analysis, "none", sep = "_"), "LKT__d__Unclassified", "LKT__Unclassified")
                rowsToKeep <- which(!(rownames(countmat) %in% dunno))
                countmat <- countmat[rowsToKeep, , drop = FALSE]
                if (nrow(countmat) < 1){
                    #abort, nothing is left over
                    flog.info("None of the wanted features were found in SummarizedExperiment object when using the current filtration parameters.")
                    return(NULL)
                }
            }

            #Eliminate non-relevant features if not adjusting p-value
            #SINGLE_GROUP: With no stats there is no p-value adjustment; cull to wantedfeatures directly.
            if (single_group || adjustpval == FALSE){
                keepnow <- wantedfeatures[wantedfeatures %in% rownames(countmat)]
                countmat <- countmat[keepnow, , drop = FALSE]
            }

            if ("GenomeCompleteness" %in% names(assays(currobj))){
                genomecompletenessdf <- as.matrix(assays(currobj)$GenomeCompleteness)
            } else {
                genomecompletenessdf <- NULL
            }

            #Aggregate if appropriate
            if (aggregatefeatures == TRUE){
                #If aggregating features, then cull to wantedfeatures now and aggregate
                wantedfeatures <- wantedfeatures[wantedfeatures %in% rownames(countmat)]
                countmat <- countmat[wantedfeatures, , drop = FALSE]
                originalwantedfeatures <- wantedfeatures
                #Count matrix should not be in log2 at this point
                aggcountmat <- colSums(countmat)
                aggcountmat <- t(as.matrix(aggcountmat))
                rownames(aggcountmat) <- aggregatefeatures_label
                countmat <- aggcountmat

                #Also aggregate Genome Completeness, if applicable
                if ("GenomeCompleteness" %in% names(assays(currobj))){
                    genomecompletenessdf <- genomecompletenessdf[wantedfeatures, , drop = FALSE]
                    agggenomecompletenessdf <- colSums(genomecompletenessdf)
                    agggenomecompletenessdf <- t(as.matrix(agggenomecompletenessdf))
                    rownames(agggenomecompletenessdf) <- aggregatefeatures_label
                    genomecompletenessdf <- agggenomecompletenessdf
                }

                #Also get mean percentage from contigs, if applicable
                if ("PctFromCtgs" %in% names(assays(currobj))){
                    PctFromCtgsdf <- PctFromCtgsdf[wantedfeatures, , drop = FALSE]
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

            #SINGLE_GROUP: In single-group mode there is nothing to compare. Build a minimal
            #matstats data frame by variance so downstream ordering / ntop / titles still work,
            #and construct a constant classesdf so the plotting scaffolding is unchanged.
            if (single_group){

                classesdf <- data.frame(Sample = colnames(countmat), cl = rep(single_group_label, ncol(countmat)), stringsAsFactors = FALSE)
                rownames(classesdf) <- classesdf$Sample
                discretenames <- single_group_label

                #Variance-based matstats purely for feature ordering; no p-values involved.
                matstats <- calculate_matrix_stats(countmatrix = countmat, uselog = log2tran_main_plot, statsonlog = FALSE, stattype = "variance", classesdf = NULL)

            } else {

                if (!is.null(wilcox_paired_by)){
                    flog.info(paste("Will attempt to pair samples by", wilcox_paired_by, "for Mann-Whitney-Wilcoxon test"))
                }

                classesdf <- make_classes_df(curr_pt = colData(currobj), compareby = compareby, wilcox_paired_by = wilcox_paired_by)

                discretenames <- sort(unique(classesdf$cl))
                if ("wilcox_pairs" %in% colnames(classesdf)){
                    flog.info(paste("Mann-Whitney-Wilcoxon test between", discretenames[1], "and", discretenames[2], "will be paired by", wilcox_paired_by))
                }

                matstats <- calculate_matrix_stats(countmatrix = countmat, uselog = log2tran_main_plot, statsonlog = FALSE, stattype = stattype, classesdf = classesdf, invertbinaryorder = invertbinaryorder, numthreads = numthreads, threshPA = threshPA, fun_for_l2fc = fun_for_l2fc)
            }

            ffeatmsg <- paste0("Number of features assessed = ", nrow(matstats))

            #Cull to only features wanted
            wantedfeatures <- wantedfeatures[wantedfeatures %in% rownames(matstats)]
            matstats <- matstats[wantedfeatures, , drop = FALSE]
            #Reorder matrix by p-value
            #SINGLE_GROUP: variance matstats has no pval column; order by SD instead.
            if (single_group){
                if ("SD" %in% colnames(matstats)){
                    matstats <- matstats[order(matstats$SD, decreasing = TRUE), , drop = FALSE]
                }
            } else {
                matstats <- matstats[order(matstats$pval), , drop = FALSE]
            }
            topcats <- nrow(matstats)
            if (!(is.null(ntop))) {
                topcats <- min(topcats, ntop)
            }

            #SINGLE_GROUP: p-value based row selection is meaningless without a comparison.
            if (!single_group && !is.null(showonlypbelow)){
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

            matstats <- matstats[rowcutoff, , drop = FALSE]

            #Filter by l2fc if applicable
            #SINGLE_GROUP: no l2fc exists without a two-class comparison; skip.
            if (!single_group && all(c((!is.null(presetlist$minl2fc)), ("absl2fc" %in% colnames(matstats))))){

                matstats <- subset(matstats, absl2fc >= presetlist$minl2fc)

                if (nrow(matstats) < 1){
                    #abort, nothing is left over
                    flog.warn("None of the wanted features were found in the SummarizedExperiment object when using the current log2 foldchange filtration parameters.")

                    return(NULL)
                }
            }

            #Filter by correlation coefficient, if applicable
            #SINGLE_GROUP: no correlation without a continuous compareby; skip.
            if (!single_group && all(c((!is.null(presetlist$minabscorrcoeff)), ("abscorrel" %in% colnames(matstats))))){
                matstats <- subset(matstats, abscorrel >= presetlist$minabscorrcoeff)

                if (nrow(matstats) < 1){
                    #abort, nothing is left over
                    flog.warn("None of the wanted features were found in the SummarizedExperiment object when using the current absolute correlation coefficient filtration parameters.")

                    return(NULL)
                }
            }

            #Redefine countmat to include only features matching filtering criteria
            countmat <- countmat[rownames(matstats), , drop = FALSE]

            if ("GenomeCompleteness" %in% names(assays(currobj))){
                #genomecompletenessdf <- as.matrix(assays(currobj)$GenomeCompleteness)
                genomecompletenessdf <- genomecompletenessdf[rownames(matstats), drop = FALSE]
                if (class(genomecompletenessdf)[1] != "matrix"){
                    genomecompletenessdf <- t(as.matrix(genomecompletenessdf))
                    rownames(genomecompletenessdf) <- rownames(matstats)
                }
            } else {
                genomecompletenessdf <- NULL
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
            if (all(c("allfeaturesbytaxa_matrix") %in% names(metadata(currobj)))){
                taxsplit_list <- retrieve_features_by_taxa(FuncExpObj = currobj, glomby = stratify_by_taxlevel, only_allow_CSBs = only_allow_CSBs, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, assay_for_matrix = "BaseCounts", wantedfeatures = featnamesforsubset, wantedsamples = colnames(countmat), asPPM = TRUE, append_metatada = TRUE, PPMthreshold = 0, include_samples_with_zero = TRUE, return_taxonomy_table = TRUE)
                taxsplit <- taxsplit_list$taxsplit

                #SINGLE_GROUP: assign a constant Compareby so the stratified plotting works unchanged.
                if (single_group){
                    taxsplit$Compareby <- single_group_label
                } else {
                    taxsplit$Compareby <- taxsplit[ , which(colnames(taxsplit) == compareby)]
                }

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
        }

        flog.info("Plotting results...")
        flog.info(paste("There are", nrow(countmat), "features to plot."))

        for (feat in rownames(countmat)){
            dat <- data.frame(Sample = classesdf$Sample, PPM = as.numeric(countmat[feat, classesdf$Sample]), Compareby = classesdf$cl, stringsAsFactors = FALSE)
            rownames(dat) <- dat$Sample

            #if there is an explicit order to compareby then set it to that
            if (!single_group && !is.null(compareby_order)){
                dat$Compareby <- factor(dat$Compareby, levels = compareby_order)
            }

            if (!is.null(shapeby)){
                dat$Shape <- curr_pt[rownames(dat) , shapeby]
            }

            if (!is.null(fillby)){
                dat$Fill <- curr_pt[rownames(dat) , fillby]
            }

            if (!is.null(connectby)){
                dat$Connect <- curr_pt[rownames(dat) , connectby]
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

            #SINGLE_GROUP: force a boxplot path; there is no spearman/pearson scatter without compareby.
            if (!single_group && matstats$Method[1] %in% c("spearman", "pearson")){
                #Make a scatterplot
                p <- p + geom_point()
                p <- p + geom_smooth(method = lm, aes(group=1), se = FALSE)
                if (!(is.null(shapeby))){
                    p <- p + aes(shape = Shape)
                    p <- add_shape_to_plot_safely(p = p, shapevec = dat$Shape, shapeby = shapeby, cdict = cdict)
                }
                rotang <- 0

            } else {
                #Code for a boxplot
                if (length(discretenames) < nrow(curr_pt)){
                    jitfact <- -( 0.3 / nrow(colData(currobj))) * (length(discretenames)) + 0.25
                } else {
                    jitfact <- 0
                }

                if (!overlay_boxplot){
                    p <- p + geom_boxplot(outlier.shape = NA)
                }

                #Jitter only if not connecting samples across boxplot
                if (!is.null(connectby)){
                    jitfact <- 0
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

                #SINGLE_GROUP: only add significance comparisons when there is a real comparison.
                if (!single_group && (length(discretenames) > 1) && (length(discretenames) <= max_pairwise_cats)){
                    if (is.null(signiflabel)){
                        signiflabel <- "p.format"
                    }
                    #Add pval
                    my_comparisons <- combn(discretenames, m = 2, simplify = FALSE)
                    p <- p + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = signiflabel)
                } else if (!single_group){
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
                        p <- p + scale_color_gradient(low = "blue", high = "red")
                    } else {
                        #if there is a colour dictionary, then use that
                        if (!(is.null(cdict))){
                            ct <- cdict[[colourby]]
                            groupcols <- setNames(as.character(ct$Colour), as.character(ct$Name))
                            p <- p + scale_color_manual(values = groupcols)
                        } else {
                            #Use colour table if available
                            if ("ctable" %in% names(metadata(currobj))){
                                discretenames <- sort(unique(dat$Colour))
                                colourshave <- discretenames[discretenames %in% rownames(metadata(currobj)$ctable)]
                                cores <- as.vector(rainbow(length(discretenames)))
                                names(cores) <- discretenames
                                cores[colourshave] <- metadata(currobj)$ctable[colourshave, "Hex"]
                                p <- p + scale_color_manual(values = cores)
                            }
                        }
                    }
                }
            }

            if (!is.null(connectby)){
                p <- p + geom_line(aes(group=Connect))
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
            #SINGLE_GROUP: no statistics to report; build a descriptive spread-only subtitle.
            if (single_group){
                nsampmsg <- paste0("n = ", ncol(countmat), " samples (single group, no comparison)")
                stattit <- paste(nsampmsg, ffeatmsg, sep = " | ")
            } else {
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
            }

            #Add description to feature, if applicable
            if (!(analysis %in% taxonomic_spaces)){
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

            if (log2tran_main_plot == TRUE){
                ytit <- "Relative Abundance in PPM"
                p <- p + scale_y_continuous(trans = scales::pseudo_log_trans(base = 2), breaks = scales::trans_breaks("log2", function(x) {((2 ^ x) - 1)}), labels = scales::trans_format("log2", function(x) {((2 ^ x) - 1)}))
            } else {
                ytit <- "Relative Abundance in PPM"
                if (!is.null(y_axis_range)){
                    p <- p + expand_limits(y=c(0, y_axis_range))
                }
            }

            #SINGLE_GROUP: label the x-axis meaningfully when there is no compareby.
            xtit <- if (single_group) single_group_label else compareby
            p <- p + labs(x = xtit, y = ytit)
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
                #Maintain only relevant columns
                annot_cols <- c("Compareby", "Shape", "Fill", "Connect", "Colour")[c("Compareby", "Shape", "Fill", "Connect", "Colour") %in% colnames(currtaxsplit)]
                LKTcolumns <- colnames(currtaxsplit)[!(colnames(currtaxsplit) %in% unique(c(colnames(curr_pt), "Accession", annot_cols)))]
                LKTsToKeep <- names(which(colSums(currtaxsplit[ , LKTcolumns]) > 0))
                currtaxsplit <- currtaxsplit[ , c("Sample", "Accession", annot_cols, LKTsToKeep)]

                #Determine how the stratified information is going to be plot
                if (use_cladogram_for_stratification != TRUE){
                    p <- NULL
                    dat <- NULL
                    currtaxsplitgrp <- currtaxsplit
                    LKTcolumns <- colnames(currtaxsplitgrp)[!(colnames(currtaxsplitgrp) %in% unique(c(colnames(curr_pt), c("Sample", "Accession", "Compareby", "Shape", "Fill", "Connect", "Colour"))))]

                    #Eliminate empties
                    LKTsToKeep <- names(which(colSums(currtaxsplitgrp[ , LKTcolumns]) > 0))
                    currtaxsplitgrp <- currtaxsplitgrp[ , c("Sample", LKTsToKeep)]
                    dat <- currtaxsplitgrp %>% gather(Taxon, PPM, 2:ncol(currtaxsplitgrp))

                    #Start building a plot
                    dat <- left_join(dat, as.data.frame(curr_pt), by = "Sample")

                    #SINGLE_GROUP: constant Compareby when there is no compareby variable.
                    if (single_group){
                        dat$Compareby <- single_group_label
                    } else {
                        dat$Compareby <- dat[ , which(colnames(dat) == compareby)]
                    }

                    #if there is an explicit order to compareby then set it to that
                    if (!single_group && !is.null(compareby_order)){
                        dat$Compareby <- factor(dat$Compareby, levels = compareby_order)
                    }

                    if (!(is.null(shapeby))){
                        dat$Shape <- dat[ , which(colnames(dat) == shapeby)]
                    }
                    if (!(is.null(colourby))){
                        dat$Colour <- dat[ , which(colnames(dat) == colourby)]
                    }

                    #Order and aggregate if more than 30
                    tally <- aggregate(PPM ~ Taxon, data = dat, FUN = "sum")
                    tally <- tally[order(tally$PPM, decreasing = TRUE), ]

                    if (return_taxon_stratification_df){
                        gvec[[plotcount]] <- as.data.frame(dat)
                        names(gvec)[plotcount] <- paste("Taxonomic_stratification_of", feat, sep = "_")
                        plotcount <- plotcount + 1
                    }

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
                        datremainder <- datremainder[!(duplicated(datremainder)), ]
                        aggremainder <- left_join(remaindertally, datremainder, by = "Sample")
                        orddat <- rbind(orddat, aggremainder[ , colnames(orddat)])
                    }

                    dat <- orddat

                    dat$Taxon <- factor(dat$Taxon, levels = unique(dat$Taxon))

                    if (!(is.null(colourby))){
                        dat$Colour <- dat[ , which(colnames(dat) == colourby)]
                    }

                    p <- ggplot(dat, aes(x = Taxon, y = PPM))

                    if (!overlay_boxplot){
                            p <- p + geom_boxplot(outlier.shape = NA)
                    }

                    #Rescale to exclude outliers
                    if (!is.null(rescale_axis_quantiles)){
                        p <- p + scale_y_continuous(limits = quantile(dat$PPM, rescale_axis_quantiles))
                    }

                    if ((length(unique(dat$Taxon))) < (nrow(dat))){
                        jitfact <- -( 0.3 / (nrow(dat))) * (length(unique(dat$Taxon))) + 0.25
                    } else {
                        jitfact <- 0
                    }

                    if (plot_points_on_taxonomy == TRUE){
                        if (!(is.null(shapeby))){
                            p <- p + geom_jitter(position = position_jitter(width = jitfact, height = 0.0), aes(shape = Shape))
                            p <- add_shape_to_plot_safely(p = p, shapevec = dat$Shape, shapeby = shapeby, cdict = cdict)
                        } else {
                            p <- p + geom_jitter(position = position_jitter(width = jitfact, height = 0.0))
                        }
                    }

                    if (overlay_boxplot){
                            p <- p + geom_boxplot(outlier.shape = NA)
                    }

                    if (!is.null(colourby)){
                        p <- p + aes(colour = Colour)

                        if (colourby == "GenomeCompleteness"){
                            p <- p + scale_fill_gradientn(aesthetics = "colour", colours = c("white", "forestgreen", "blue", "firebrick1", "black"),  values = scales::rescale(c(0, 100, 200, 300, 400), to = c(0, (400/max(dat$Colour)))))
                        } else {
                            if (is.numeric(dat$Colour)){
                                p <- p + scale_color_gradient(low = "blue", high = "red")
                            } else {
                                #if there is a colour dictionary, then use that
                                if (!(is.null(cdict))){
                                    ct <- cdict[[colourby]]
                                    groupcols <- setNames(as.character(ct$Colour), as.character(ct$Name))
                                    p <- p + scale_color_manual(values = groupcols)
                                } else {
                                    #Use colour table if available
                                    if ("ctable" %in% names(metadata(currobj))){
                                        discretenames <- sort(unique(dat$Colour))
                                        colourshave <- discretenames[discretenames %in% rownames(metadata(currobj)$ctable)]
                                        cores <- as.vector(rainbow(length(discretenames)))
                                        names(cores) <- discretenames
                                        cores[colourshave] <- metadata(currobj)$ctable[colourshave, "Hex"]
                                        p <- p + scale_color_manual(values = cores)
                                    }
                                }
                            }
                        }
                    }

                    if (horizontal == TRUE){
                        p <- p + coord_flip()
                    }

                    #SINGLE_GROUP: no need to facet by a constant single group.
                    if (!single_group){
                        if (wrap_facet){
                            p <- p + facet_wrap( ~ Compareby)
                        } else {
                            p <- p + facet_grid( ~ Compareby)
                        }
                    }

                    p <- p + theme_minimal()
                    p <- p + theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
                    plotitstrat <- paste0(c(maintit, featname), collapse = "\n")

                    p <- p + ggtitle(plotitstrat)

                    if (!(is.null(colourby))){
                        p <- p + labs(colour = colourby)
                    }

                    if (!(is.null(shapeby))){
                        p <- p + labs(shape = shapeby)
                    }

                    if (log2tran_strat_plot == TRUE){
                        ytit <- "Relative Abundance in PPM"
                        p <- p + scale_y_continuous(trans = scales::pseudo_log_trans(base = 2))

                    } else {
                        ytit <- "Relative Abundance in PPM"
                    }
                    p <- p + labs(x = "Contributing Taxon", y = ytit)

                    p <- p + theme(axis.text.x = element_text(colour = "black", angle = rotang, size = rel(0.85)))
                    p <- p + theme(plot.title = element_text(size = 10))

                } else {
                    #use cladogram stratification plotting function
                    plotitstrat <- paste0(c(maintit, featname), collapse = "\n")
                    data(PhyCols)
                    p <- plot_taxsplit_tree(taxsplit_df = currtaxsplit, strat_tt = taxsplit_list$taxtable, layout = "rectangular", taxon_hilight_palette = PhyCols, low_col = "#fbfbfc", high_col = "#1b0178", show_prevalence = show_prevalence_in_cladogram, plotitstrat = plotitstrat)
                }#End if statement for plotting as boxplot or cladogram

                if (!return_plots){
                    #print plot on the fly
                    print(p)
                }
                #Bank plot to gvec
                gvec[[plotcount]] <- p
                names(gvec)[plotcount] <- paste(maintit, feat, stratify_by_taxlevel, sep = " | ")
                plotcount <- plotcount + 1
            }#End conditional of stratifying by taxonomy
        }#End loop for plotting each feature
    }#End loop for each subset

    if (return_plots){
        #Return plots, as nothing was printed
        return(gvec)
    }
}