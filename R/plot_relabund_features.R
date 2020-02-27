#' plot_relabund_features(ExpObj = NULL, glomby = NULL, samplesToKeep = NULL, featuresToKeep = NULL, aggregatefeatures = FALSE, aggregatefeatures_label = "Sum_of_wanted_features", subsetby = NULL, compareby = NULL, colourby = NULL, shapeby = NULL, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL, ntop = NULL, adjustpval = TRUE, padjmeth = "fdr", showonlypbelow = NULL, showonlypadjusted = FALSE, maxl2fc = NULL, minl2fc = NULL, addtit = NULL, PPM_normalize_to_bases_sequenced = FALSE, uselog = FALSE, statsonlog = FALSE, cdict = NULL, maxnumplots = NULL, signiflabel = "p.format", max_pairwise_cats = 4, numthreads = 1, nperm = 99, ignoreunclassified = TRUE, class_to_ignore = "N_A", ...)
#'
#' Generates relative abundance plots per feature annotated by the metadata using as input a SummarizedExperiment object
#' @export

plot_relabund_features <- function(ExpObj = NULL, glomby = NULL, samplesToKeep = NULL, featuresToKeep = NULL, aggregatefeatures = FALSE, aggregatefeatures_label = "Sum_of_wanted_features", subsetby = NULL, compareby = NULL, colourby = NULL, shapeby = NULL, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL, ntop = NULL, adjustpval = TRUE, padjmeth = "fdr", showonlypbelow = NULL, showonlypadjusted = FALSE, maxl2fc = NULL, minl2fc = NULL, addtit = NULL, PPM_normalize_to_bases_sequenced = FALSE, uselog = FALSE, statsonlog = FALSE, cdict = NULL, maxnumplots = NULL, signiflabel = "p.format", max_pairwise_cats = 4, numthreads = 1, nperm = 99, ignoreunclassified = TRUE, class_to_ignore = "N_A", ...){

    variables_to_fix <- c(compareby, subsetby, colourby, shapeby)

    #Vet experiment object
    obj <- ExpObjVetting(ExpObj = ExpObj, samplesToKeep = samplesToKeep, featuresToKeep = NULL, glomby = glomby, variables_to_fix = variables_to_fix, class_to_ignore = class_to_ignore)

    analysis <- metadata(obj)$analysis
    analysisname <- analysis

    presetlist <- declare_filtering_presets(analysis = analysis, applyfilters = applyfilters, featcutoff = featcutoff, GenomeCompletenessCutoff = GenomeCompletenessCutoff, PctFromCtgscutoff = PctFromCtgscutoff, maxl2fc = maxl2fc, minl2fc = minl2fc)

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
                    flog.info("None of the wanted features were not found in SummarizedExperiment object when using the current filtration parameters.")
                    return(NULL)
                }
            }

            #Aggregate if appropriate
            if (aggregatefeatures == TRUE){
                #Count matrix should not be in log2 at this point
                aggcountmat <- colSums(countmat)
                aggcountmat <- t(as.matrix(aggcountmat))
                rownames(aggcountmat) <- aggregatefeatures_label
                countmat <- aggcountmat
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

            if (length(rowcutoff) == 0){
                #abort, nothing is left over
                flog.warn("None of the wanted features were not found in SummarizedExperiment object when using the current p-value filtration parameters.")

                return(NULL)
            }

            matstats <- matstats[rowcutoff, ]

            #Filter by l2fc if applicable
            if (all(c((!is.null(minl2fc)), ("absl2fc" %in% colnames(matstats))))){
                matstats <- subset(matstats, absl2fc >= minl2fc)
                if (nrow(matstats) < 1){
                    #abort, nothing is left over
                    flog.warn("None of the wanted features were not found in SummarizedExperiment object when using the current log2 foldchange filtration parameters.")

                    return(NULL)
                }
            }

            #Redefine countmat to include only features matching filtering criteria
            countmat <- countmat[rownames(matstats), ]
            if (class(countmat) != "matrix"){
                countmat <- t(as.matrix(countmat))
                rownames(countmat) <- wantedfeatures
            }

            if ("GenomeCompleteness" %in% names(assays(currobj))){
                genomecompletenessdf <- as.matrix(assays(currobj)$GenomeCompleteness)
                genomecompletenessdf <- genomecompletenessdf[wantedfeatures, ]
                if (class(genomecompletenessdf) != "matrix"){
                    genomecompletenessdf <- t(as.matrix(genomecompletenessdf))
                    rownames(genomecompletenessdf) <- wantedfeatures
                }
            } else {
                genomecompletenessdf <- NULL
            }

            if ("PctFromCtgs" %in% names(assays(currobj))){
                PctFromCtgsdf <- as.matrix(assays(currobj)$PctFromCtgs)
                PctFromCtgsdf <- PctFromCtgsdf[wantedfeatures, ]
                if (class(PctFromCtgsdf) != "matrix"){
                    PctFromCtgsdf <- t(as.matrix(PctFromCtgsdf))
                    rownames(PctFromCtgsdf) <- wantedfeatures
                }
            } else {
                PctFromCtgsdf <- NULL
            }

        } else {

            #abort, nothing is left over
            flog.warn("None of the wanted features were not found in SummarizedExperiment object when using the current filtration parameters.")

            return(NULL)
        }

        for (feat in rownames(countmat)){
            dat <- data.frame(Sample = names(countmat[feat, ]), PPM = as.numeric(countmat[feat, ]), Compareby = cl, stringsAsFactors = FALSE)
            rownames(dat) <- dat$Sample

            if (!is.null(shapeby)){
                dat$Shape <- curr_pt[rownames(dat) , shapeby]
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
                    p <- p + aes(shape = Shape)
                    numshapes <- length(unique(dat$Shape))
                    p <- p + scale_shape_manual(values = 15:(numshapes + 15))
                }
                rotang <- 0

            } else {
                #Code for a boxplot
                if(length(discretenames) < nrow(curr_pt)){
                    jitfact <- -( 0.3 / nrow(colData(currobj))) * (length(discretenames)) + 0.25
                } else {
                    jitfact <- 0
                }

                p <- p + geom_boxplot(outlier.shape = NA)

                if (!(is.null(shapeby))){
                    p <- p + geom_jitter(position = position_jitter(width = jitfact, height = 0.0), aes(shape = Shape))
                } else {
                    p <- p + geom_jitter(position = position_jitter(width = jitfact, height = 0.0))
                }

                if ((length(discretenames) > 1) && (length(discretenames) <= max_pairwise_cats)){
                    if (missing(signiflabel)){
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
            p <- p + theme_minimal()
            #Build plot title
            overallpmeth <- matstats[feat, "Method"]
            overallp <- paste0("pval=", round(matstats[feat, "pval"], 4))
            overalladjp <- paste0("padj_fdr=", round(matstats[feat, "padj_fdr"], 4))
            stattit <- paste(overallpmeth, overallp, overalladjp, ffeatmsg, sep = " | ")

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
                p <- p + scale_y_continuous(trans = scales::log2_trans(), breaks = scales::trans_breaks("log2", function(x) {2^x}), labels = scales::trans_format("log2", function(x) {2^x}))
            } else {
                ytit <- "Relative Abundance in PPM"
            }

            p <- p + labs(x = compareby, y = ytit)
            p <- p + theme(axis.text.x = element_text(angle = rotang, size = rel(1), colour = "black"))
            p <- p + theme(plot.title = element_text(size = 10))

            gvec[[plotcount]] <- p
            names(gvec)[plotcount] <- paste(maintit, feat, sep = " | ")
            plotcount <- plotcount + 1

        }#End loop for plotting each feature

    }#End loop for each subset

    return(gvec)

}
