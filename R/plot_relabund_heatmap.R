#' plot_relabund_heatmap(ExpObj = NULL, glomby = NULL, hmtype = NULL, samplesToKeep = NULL, featuresToKeep = NULL, subsetby = NULL, compareby = NULL, invertbinaryorder = FALSE, hmasPA = FALSE, threshPA = 0, ntop = NULL, splitcolsby = NULL, ordercolsby = NULL, cluster_samples_per_heatmap = FALSE, cluster_features_per_heatmap = FALSE, colcategories = NULL, textby = NULL, cluster_rows = TRUE, max_rows_in_heatmap = 50, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL, discard_SDoverMean_below = NULL, maxl2fc = NULL, minl2fc = NULL, adjustpval = FALSE, showonlypbelow = NULL, showpval = TRUE, showl2fc = TRUE, showGram = TRUE, secondaryheatmap = "GenomeCompleteness", addtit = NULL, PPM_normalize_to_bases_sequenced = FALSE, scaled = FALSE, cdict = NULL, maxnumheatmaps = NULL, numthreads = 1, nperm = 99, statsonlog = FALSE, ignoreunclassified = TRUE, returnstats = FALSE, class_to_ignore = "N_A", ...)
#'
#' Plots relative abundance heatmaps annotated by the metadata using as input a SummarizedExperiment object
#' @export

plot_relabund_heatmap <- function(ExpObj = NULL, glomby = NULL, hmtype = NULL, samplesToKeep = NULL, featuresToKeep = NULL, subsetby = NULL, compareby = NULL, invertbinaryorder = FALSE, hmasPA = FALSE, threshPA = 0, ntop = NULL, splitcolsby = NULL, ordercolsby = NULL, cluster_samples_per_heatmap = FALSE, cluster_features_per_heatmap = FALSE, colcategories = NULL, textby = NULL, cluster_rows = TRUE, max_rows_in_heatmap = 50, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL, discard_SDoverMean_below = NULL, maxl2fc = NULL, minl2fc = NULL, adjustpval = FALSE, showonlypbelow = NULL, showpval = TRUE, showl2fc = TRUE, showGram = TRUE, secondaryheatmap = "GenomeCompleteness", addtit = NULL, PPM_normalize_to_bases_sequenced = FALSE, scaled = FALSE, cdict = NULL, maxnumheatmaps = NULL, numthreads = 1, nperm = 99, statsonlog = FALSE, ignoreunclassified = TRUE, returnstats = FALSE, class_to_ignore = "N_A", ...){

    #Test for silly stuff
    if ((hmtype %in% c("comparative", "PA")) && (is.null(compareby))){
        stop("If rows are to be selected by highest significant difference between classes in a discrete category, you must determine the category using the argument *classesvector*")
    }

    #Fix metadata classes and remove classes to ignore, if present
    if (hmtype %in% c("comparative", "PA")){
        variables_to_fix <- c(compareby, subsetby)
    } else {
        variables_to_fix <- subsetby
    }

    #Vet experiment object
    obj <- ExpObjVetting(ExpObj = ExpObj, samplesToKeep = samplesToKeep, featuresToKeep = featuresToKeep, glomby = glomby, variables_to_fix = variables_to_fix, class_to_ignore = class_to_ignore)

    analysis <- metadata(obj)$analysis
    if (!is.null(glomby)){
        analysisname <- glomby
    } else {
        analysisname <- analysis
    }

    presetlist <- declare_filtering_presets(analysis = analysis, applyfilters = applyfilters, featcutoff = featcutoff, GenomeCompletenessCutoff = GenomeCompletenessCutoff, PctFromCtgscutoff = PctFromCtgscutoff, maxl2fc = maxl2fc, minl2fc = minl2fc)

    if (!(is.null(subsetby))){
        subset_points <- sort(unique(colData(obj)[, which(colnames(colData(obj)) == subsetby)]))
    } else {
        subset_points <- "none"
    }

    if (is.null(colcategories)){
        if (!is.null(cdict)){
            colcategories <- names(cdict)
        } else {
            #Include anything with between 2 and 10 classes in it
            colcategories <- colnames(colData(obj))[which((sapply(1:ncol(colData(obj)), function (x) { length(levels(as.factor(colData(obj)[, x]))) })) < 10 & (sapply(1:ncol(colData(obj)), function (x) { length(levels(as.factor(colData(obj)[, x]))) })) > 1)]
        }
    }

    #Initialize Stats Vector list
    svec <- NULL
    svec <- vector("list", length = 1000)
    s <- 1
    n <- 1

    #subset by metadata column
    for (sp in 1:length(subset_points)){
        if (!(is.null(subsetby))){
            samplesToKeep <- rownames(colData(obj))[which(colData(obj)[ , subsetby] == subset_points[sp])]
            flog.info(paste("Plotting within", subset_points[sp]))
            colcategories <- colcategories[!(colcategories %in% subsetby)]
            subsetname <- subset_points[sp]
        } else {
            samplesToKeep <- rownames(colData(obj))
            subsetname <- "no_sub"
        }

        #See if there are enough samples and features to go ahead
        proceed <- TRUE
        curr_pt <- colData(obj)[samplesToKeep, ]

        if ((dim(curr_pt)[1] * dim(curr_pt)[2]) < 4){
            #There are less than 4 cells, a heatmap is meaningless.
            proceed <- FALSE
        }

        if (all(c(!is.null(compareby), (hmtype != "exploratory")))) {
            #investigate if there are at least two of each class of things to compare to go ahead
            if (length(unique(curr_pt[ , compareby])) == 2){
                #Comparison is binary. Is there at least 2 of each class to get a p-value?
                if ((min(table(curr_pt[ , compareby]))) < 2){
                    proceed <- FALSE
                }
            }
        }

        if (proceed){

            hmtypemsg <- "Relative Abundance Heatmap"
            if (hmtype == "exploratory"){
                stattype <- "variance"
                if (is.null(ntop)){
                    ntop <- max_rows_in_heatmap
                }
            } else if (hmtype == "PA"){
                stattype <- "PA"
            } else {
                stattype <- "auto"
            }

            currobj <- filter_experiment(ExpObj = obj, featcutoff = presetlist$featcutoff, samplesToKeep = samplesToKeep, featuresToKeep = featuresToKeep, asPPM = TRUE, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, GenomeCompletenessCutoff = presetlist$GenomeCompletenessCutoff, PctFromCtgscutoff = presetlist$PctFromCtgscutoff, discard_SDoverMean_below = discard_SDoverMean_below)

        } else {
            flog.info("Unable to make heatmap with the current metadata for this comparison.")
            return(NULL)
        }

        #There must be at least two samples for a heatmap and at least two features
        if ((nrow(currobj) * ncol(currobj)) >= 4){

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
            }

            #Rename rows to include description if not taxonomic data or MetaCyc which has enormous descriptions
            #if (!(analysis %in% c("LKT", "MetaCyc"))){
            if (analysis  != "LKT"){
                feattable <- rowData(currobj)
                feattable$Feature <- paste(feattable$Accession, feattable$Description, sep = "-")
                rownames(countmat) <- feattable$Feature[match(rownames(countmat), feattable$Accession)]
            }
            matrixSamples <- colnames(countmat)
            matrixRows <- rownames(countmat)

            topcats <- nrow(countmat)
            if (!(is.null(ntop))) {
                topcats <- min(topcats, ntop)
            }

            #Calculate matrix stats and get new matrix.
            if (stattype == "variance"){
                matstats <- calculate_matrix_stats(countmatrix = countmat, stattype = "variance")
                #Bank raw stats
                matstatsbank <- as.data.frame(matstats)
                if ("GenomeCompleteness" %in% names(assays(currobj))){
                    matstatsbank$Taxa <- rownames(matstatsbank)
                    genomecompletenessmat <- as.matrix(assays(currobj)$GenomeCompleteness)
                    genomecompletenessdf <- as.data.frame(genomecompletenessmat)
                    genomecompletenessdf$MedianGenomeComp <- rowMedians(genomecompletenessmat)
                    genomecompletenessdf$SDGenomeComp <- rowSds(genomecompletenessmat)
                    genomecompletenessdf$Taxa <- rownames(genomecompletenessdf)
                    genomecompletenessdfmedians <- genomecompletenessdf[, c("Taxa", "MedianGenomeComp", "SDGenomeComp")]
                    matstatsbank <- left_join(matstatsbank, genomecompletenessdfmedians, by = "Taxa")
                    rownames(matstatsbank) <- matstatsbank$Taxa
                    matstatsbank$Taxa <- NULL
                }
                svec[[s]] <- matstatsbank
                ffeatmsg <- paste0("Number of features assessed = ", nrow(matstats))

                matstats$Colour <- rep("black", nrow(matstats))
                countmat2 <- countmat[rownames(matstats), ]

                #Transform to log2 space
                countmat2 <- convert_matrix_log2(mat = countmat2, transformation = "to_log2")

                if (all(c(cluster_rows, (any(!(c(cluster_samples_per_heatmap, cluster_features_per_heatmap))))))){
                    flog.info("Clustering samples and features using entire matrix to obtain sample and feature order for all heatmaps.")
                    #create a heatmap from the entire count matrix for getting column order.
                    htfull <- Heatmap(countmat2, name = "FullHM", column_title = "FullHM", column_names_gp = gpar(fontsize = 1), cluster_rows = TRUE, show_row_dend = FALSE, row_names_side = "left", row_names_gp = gpar(fontsize = 1, col = "black"), use_raster = TRUE)

                    fullheatmap_column_order <- column_order(htfull)
                    fullheatmap_column_dend <- column_dend(htfull)
                    fullheatmap_row_order <- row_order(htfull)
                    fullheatmap_row_dend <- row_dend(htfull)
                }

                if (!cluster_features_per_heatmap){
                    countmat2 <- countmat2[fullheatmap_row_order, ]
                    matstats <- matstats[rownames(countmat2), ]
                }

                #Create a list of matrices each of maximum 50 rows
                rowlist <- split(1:topcats, ceiling(seq_along(1:topcats) / max_rows_in_heatmap))
                matlist <- lapply(1:length(rowlist), function(x){countmat2[rowlist[[x]], ]})
                statslist <- lapply(1:length(rowlist), function(x){ matstats[rowlist[[x]], ] })
                rowlblcol_list <- lapply(1:length(rowlist), function(x){ rep("black", length(rowlist[[x]])) })
                stattit <- paste("Top", topcats, "most variant features across samples")
                statmsg <- stattype

            } else {
                cl <- colData(currobj)[ , which(colnames(colData(currobj)) == compareby)]
                discretenames <- sort(unique(cl))

                matstats <- calculate_matrix_stats(countmatrix = countmat, uselog = FALSE, statsonlog = statsonlog, stattype = stattype, classesvector = cl, invertbinaryorder = invertbinaryorder, numthreads = numthreads, nperm = nperm, threshPA = threshPA)

                ffeatmsg <- paste0("Number of features assessed = ", nrow(matstats))

                #Bank raw stats
                matstatsbank <- as.data.frame(matstats)
                if ("GenomeCompleteness" %in% names(assays(currobj))){
                    matstatsbank$Taxa <- rownames(matstatsbank)
                    genomecompletenessmat <- as.matrix(assays(currobj)$GenomeCompleteness)
                    genomecompletenessdf <- as.data.frame(genomecompletenessmat)
                    genomecompletenessdf$MedianGenomeComp <- rowMedians(genomecompletenessmat)
                    genomecompletenessdf$SDGenomeComp <- rowSds(genomecompletenessmat)
                    genomecompletenessdf$Taxa <- rownames(genomecompletenessdf)
                    genomecompletenessdfmedians <- genomecompletenessdf[, c("Taxa", "MedianGenomeComp", "SDGenomeComp")]
                    matstatsbank <- left_join(matstatsbank, genomecompletenessdfmedians, by = "Taxa")
                    rownames(matstatsbank) <- matstatsbank$Taxa
                    matstatsbank$Taxa <- NULL
                }
                svec[[s]] <- matstatsbank

                #Account for the fact that if there is only a single class within compareby variable, stats may have been coerced to variance.
                if (matstats$Method[1] == "variance") {
                    matstats$Colour <- rep("black", nrow(matstats))
                    countmat2 <- countmat[rownames(matstats), ]
                    #Transform to log2 space
                    countmat2 <- convert_matrix_log2(mat = countmat2, transformation = "to_log2")

                    if (all(c(cluster_rows, (any(!(c(cluster_samples_per_heatmap, cluster_features_per_heatmap))))))){
                        flog.info("Clustering samples and features using entire matrix to obtain sample and feature order for all heatmaps.")
                        #create a heatmap from the entire count matrix for getting column order.
                        htfull <- Heatmap(countmat2, name = "FullHM", column_title = "FullHM", column_names_gp = gpar(fontsize = 1), cluster_rows = TRUE, show_row_dend = FALSE, row_names_side = "left", row_names_gp = gpar(fontsize = 1, col = "black"), use_raster = TRUE)

                        fullheatmap_column_order <- column_order(htfull)
                        fullheatmap_column_dend <- column_dend(htfull)
                        fullheatmap_row_order <- row_order(htfull)
                        fullheatmap_row_dend <- row_dend(htfull)
                    }

                    if (all(c(cluster_rows, (!cluster_features_per_heatmap)))){
                        countmat2 <- countmat2[fullheatmap_row_order, ]
                        matstats <- matstats[rownames(countmat2), ]
                    }

                    #Create a list of matrices each of maximum 50 rows
                    topcats <- min(nrow(countmat2), max_rows_in_heatmap)
                    rowlist <- split(1:topcats, ceiling(seq_along(1:topcats) / max_rows_in_heatmap))
                    matlist <- lapply(1:length(rowlist), function(x){ countmat2[rowlist[[x]], ] })
                    statslist <- lapply(1:length(rowlist), function(x){ matstats[rowlist[[x]], ] })
                    rowlblcol_list <- lapply(1:length(rowlist), function(x){rep("black", length(rowlist[[x]]))})
                    stattit <- paste("Top", topcats, "most variant features across samples")
                    statmsg <- paste("Var", compareby, sep = "_")

                } else {
                    #Stats are not variance
                    if (all(c((!(is.null(presetlist$minl2fc))), ("l2fc" %in% colnames(matstats)), (hmtype != "PA")))) {
                        #Check if there is enough leftover after filter.
                        if (length(which(matstats$absl2fc > presetlist$minl2fc)) < 2){
                            flog.info(paste("There are less than 2 features which have >", presetlist$minl2fc, "l2fc."))
                            #Redefine min l2fc to whatever is the lowest available keeping between 2 and 40 features to plot.
                            presetlist$minl2fc <- matstats[order(matstats$absl2fc, decreasing = TRUE),][min(max(length(matstats$absl2fc), 2), 40), ]$absl2fc
                            flog.info(paste("Resetting minimum l2fc to", round(presetlist$minl2fc, 2)))
                        }
                        minl2fcmsg <- paste("log2foldchange >", round(presetlist$minl2fc, 2))
                        matstats <- subset(matstats, absl2fc > presetlist$minl2fc)
                    } else {
                        minl2fcmsg <- "log2foldchange > 0"
                    }

                    if (all(c((!(is.null(presetlist$maxl2fc))), ("l2fc" %in% colnames(matstats)), (hmtype != "PA")))) {
                        #Check if there is enough leftover after filter.
                        if (length(which(matstats$absl2fc < presetlist$maxl2fc)) < 2){
                            flog.info(paste("There are less than 2 features which have <", presetlist$maxl2fc, "l2fc."))
                            #Redefine max l2fc to whatever is the lowest available keeping between 2 and 40 features to plot.
                            presetlist$maxl2fc <- matstats[order(matstats$absl2fc, decreasing = FALSE),][min(max(length(matstats$absl2fc), 2), 40), ]$absl2fc
                            flog.info(paste("Resetting maximum l2fc to", round(presetlist$maxl2fc, 2)))
                        }
                        maxl2fcmsg <- paste("log2foldchange <", round(presetlist$maxl2fc, 2))
                        matstats <- subset(matstats, absl2fc < presetlist$maxl2fc)
                    } else {
                        maxl2fcmsg <- "log2foldchange > Inf"
                    }

                    if ("l2fc" %in% colnames(matstats)){
                        matstats$Colour <- ifelse(matstats$l2fc < 0, "#900000", "#000000")
                        statmsg <- paste("MWW", compareby, sep="_")
                    } else if ("OddsRatio" %in% colnames(matstats)){
                        matstats$Colour <- ifelse(matstats$OddsRatio < 1, "#900000", "#000000")
                        statmsg <- paste("Fisher", compareby, sep="_")
                    } else {
                        matstats$Colour <- rep("#000000", nrow(matstats))

                        if (matstats$Method[1] == "permanova"){
                            statmsg <- paste("PERMANOVA", compareby, sep="_")
                        } else if (matstats$Method[1] == "anova"){
                            statmsg <- paste("ANOVA", compareby, sep="_")
                        } else if (matstats$Method[1] == "variance"){
                            statmsg <- paste("Var", compareby, sep="_")
                        }
                    }

                    if (matstats$Method[1] == "fisher"){
                        binary_directionality <- paste("Odds Ratio > 1 means enriched in", (discretenames[1:2][c(!invertbinaryorder, invertbinaryorder)]))
                    } else {
                        binary_directionality <- paste("Positive l2fc means increased in", (discretenames[1:2][c(!invertbinaryorder, invertbinaryorder)]))
                    }

                    #Obtain a matrix that represents cells in the heatmap
                    countmat2 <- as.matrix(countmat[rownames(matstats), ])
                    if (all(c((hmtype == "PA"), hmasPA))){
                        #Transform to Presence/Absence according to threshold space
                        countmat2 <- convert_matrix_PA(mat = countmat2, threshPA = threshPA)
                    } else {
                        #Transform to log2 space
                        countmat2 <- convert_matrix_log2(mat = countmat2, transformation = "to_log2")
                    }

                    #Decide row and sample ordering
                    if (all(c(cluster_rows, (any(!(c(cluster_samples_per_heatmap, cluster_features_per_heatmap))))))){
                        flog.info("Clustering samples and features using entire matrix to obtain sample and feature order for all heatmaps.")
                        #create a heatmap from the entire count matrix for getting column order.
                        htfull <- Heatmap(countmat2, name = "FullHM", column_title = "FullHM", column_names_gp = gpar(fontsize = 1), cluster_rows = TRUE, show_row_dend = FALSE, row_names_side = "left", row_names_gp = gpar(fontsize = 1, col = "black"), use_raster = TRUE)

                        fullheatmap_column_order <- column_order(htfull)
                        fullheatmap_column_dend <- column_dend(htfull)
                        fullheatmap_row_order <- row_order(htfull)
                        fullheatmap_row_dend <- row_dend(htfull)
                    }

                    if (all(c(cluster_rows, (!cluster_features_per_heatmap)))) {
                        countmat2 <- countmat2[fullheatmap_row_order, ]
                        matstats <- matstats[rownames(countmat2), ]
                    }

                    #If adjustpval is set to auto, then find out which is best and re-set it to either TRUE or FALSE
                    if (("pval" %in% colnames(matstats)) && adjustpval == "auto"){
                        propsigadj <- length(which(matstats$padj_fdr < 0.05)) / length(matstats$padj_fdr)
                        propsignonadj <- length(which(matstats$pval < 0.05)) / length(matstats$pval)
                        fracsigadj <- propsigadj / propsignonadj
                        if ((!is.na(fracsigadj)) && (fracsigadj > 0.2)){
                            adjustpval <- TRUE
                        } else {
                            adjustpval <- FALSE
                        }
                        flog.info(paste("P-value adjustment set to auto. Upon adjustment using FDR, the proportion still below 0.05 after adjustment is", round((fracsigadj * 100), 1), "% of unadjusted p-values, so p-value adjustment is set to", as.character(adjustpval)))
                    }

                    if (("pval" %in% colnames(matstats)) && adjustpval != TRUE){
                        sigmeas <- "pval"
                    } else if (("pval" %in% colnames(matstats)) && adjustpval == TRUE){
                        sigmeas <- "padj_fdr"
                    }

                    if (is.null(showonlypbelow)){
                        if (nrow(countmat2) < topcats){
                            topcats <- nrow(countmat2)
                            flog.info("There are not enough features matching the criteria imposed.")
                            flog.info(paste("Showing top", topcats, "features."))
                        }
                        countmat2 <- countmat2[1:topcats, ]
                        #Create a list of matrices each of maximum ~50 rows
                        if((topcats %% max_rows_in_heatmap) == 1){
                            chunksize <- (max_rows_in_heatmap - 2)
                        } else {
                            chunksize <- max_rows_in_heatmap
                        }

                        rowlist <- split(1:topcats, ceiling(seq_along(1:topcats) / chunksize))
                        matlist <- lapply(1:length(rowlist), function(x){ countmat2[rowlist[[x]], ] })
                        statslist <- lapply(1:length(rowlist), function(x){ matstats[rowlist[[x]], ] })
                        rowlblcol_list <- lapply(1:length(rowlist), function(x) { matstats$Colour[rowlist[[x]]] })

                        if (matstats$Method[1] == "MannWhitneyWilcoxon") {
                            stattit <- paste("Top", topcats, "different between", compareby, "using MannWhitneyWilcoxon")
                        } else if (matstats$Method[1] == "permanova"){
                            stattit <- paste("Top", topcats, "different between", compareby, "using PERMANOVA")
                        } else if (matstats$Method[1] == "anova"){
                            stattit <- paste("Top", topcats, "different between", compareby, "using ANOVA")
                        } else if (matstats$Method[1] == "fisher"){
                            stattit <- paste("Top", topcats, "Present/Absent between", compareby, "using Fishers test")
                        }

                    } else {
                        if (adjustpval != TRUE){
                            rowcutoff <- which(matstats$pval < showonlypbelow)
                        } else {
                            rowcutoff <- which(matstats$padj_fdr < showonlypbelow)
                        }

                        #Limit number of features to requested number or number available
                        rowcutoff <- rowcutoff[1:(min(topcats, length(rowcutoff)))]

                        #Must have at least two rows in a matrix to plot a heatmap
                        if (length(rowcutoff) > 1){
                            #Account for the fact that if the remainder of a chunk of 50 is 1 then a heatmap with a single feature cannot be drawn. Oh my, I have seen everything havent I...
                            if((length(rowcutoff) %% max_rows_in_heatmap) == 1){
                                chunksize <- (max_rows_in_heatmap - 2)
                            } else {
                                chunksize <- max_rows_in_heatmap
                            }

                            rowlist <- split(rowcutoff, ceiling(seq_along(rowcutoff) / chunksize))
                            matlist <- lapply(1:length(rowlist), function(x){ countmat2[rowlist[[x]], ] })
                            statslist <- lapply(1:length(rowlist), function(x){ matstats[rowlist[[x]], ] })
                            if (matstats$Method[1] == "MannWhitneyWilcoxon") {
                                stattit <- paste(sigmeas, "<", showonlypbelow, "different between", compareby, "using MannWhitneyWilcoxon")
                            } else if (matstats$Method[1] == "permanova"){
                                stattit <- paste(sigmeas, "<", showonlypbelow, "different between", compareby, "using PERMANOVA")
                            } else if (matstats$Method[1] == "anova"){
                                stattit <- paste(sigmeas, "<", showonlypbelow, "different between", compareby, "using ANOVA")
                            } else if (matstats$Method[1] == "fisher"){
                                stattit <- paste(sigmeas, "<", showonlypbelow, "Present/Absent between", compareby, "using Fishers test")
                            }
                        } else {
                            #allIhave<-min(nrow(countmat2), 30)
                            #stattit<-paste("No values", "p <", showonlypbelow, "between", compareby,"\nShowing top", allIhave, "with any p")
                            #matlist<-list(countmat2[1:allIhave, ])
                            #rowlblcol_list<-list(matstats$Colour[1:allIhave])
                            #Report nothing was found
                            matlist <- NULL
                            if (!(is.na(matstats$Method[1]))) {
                                if (matstats$Method[1] == "MannWhitneyWilcoxon") {
                                    stattit <- c(paste(sigmeas, "<", showonlypbelow, "different between"), paste(compareby, "using MannWhitneyWilcoxon"))
                                } else if (matstats$Method[1] == "permanova"){
                                    stattit <- c(paste(sigmeas, "<", showonlypbelow, "different between"), paste(compareby, "using PERMANOVA"))
                                } else if (matstats$Method[1] == "anova"){
                                    stattit <- c(paste(sigmeas, "<", showonlypbelow, "different between"), paste(compareby, "using ANOVA"))
                                } else if (matstats$Method[1] == "fisher"){
                                    stattit <- c(paste(sigmeas, "<", showonlypbelow, "Present/Absent between"), paste(compareby, "using Fishers test"))
                                    if (threshPA != 0){
                                        stattit <- c(stattit, paste("Presence is considered PPM >=", threshPA))
                                    }
                                } #End conditional for getting stats title
                            } else {
                                #Stats failed for some reason, like there are no stat features fulfilling the filtering criteria
                                stattit <- NULL
                            }
                        } #End conditional that there are two rows or more to plot
                    } #End conditional for only showing feats with certain pval
                } #End conditional that stats method is not variance
            } #End conditional for choosing stats method

            #Name stats in svec
            stattitle <- paste(analysisname, statmsg, subsetname, sep = "_")
            names(svec)[s] <- stattitle
            s <- s + 1

            #If the matrix list is not empty, plot the heatmap. If there is nothing to plot, say so.
            if (length(matlist) > 0){
                #Cycle through list of matrices to transform into heatmaps.
                if (!(is.null(maxnumheatmaps))){
                    #Prune graphics list to largest number allowed
                    mhm <- min(maxnumheatmaps, length(matlist))
                } else {
                    mhm <- length(matlist)
                }

                for (hm in 1:mhm){
                    #Regenerate the current matrix being plot from matrix list
                    mathm <- matlist[[hm]]
                    rownames(mathm) <- strtrim(rownames(mathm), 60)
                    stathm <- statslist[[hm]]
                    rowlblcol <- stathm$Colour
                    ht1fs <- 10

                    #Plot the heatmap
                    fontsizey <- min(6, round((((-1 / 150) * (nrow(mathm))) + 1) * 5, 1))
                    fontsizex <- as.numeric(unlist(round((((-1 / 150) * (ncol(mathm))) + 1) * 5, 0)))

                    hmdf <- as.data.frame(matrix(data = 0, nrow = nrow(colData(currobj)), ncol = length(colcategories)))
                    cores <- vector("list", length = length(colcategories))
                    for (g in 1:length(colcategories)){
                        hmdf[ , g] <- colData(currobj)[ , which(colnames(colData(currobj)) == colcategories[g])]
                        colnames(hmdf)[g] <- colcategories[g]

                        #Test if variable can be coerced to numeric
                        if (!(can_be_made_numeric(hmdf[ , g], cats_to_ignore = class_to_ignore))){
                            if (is.null(cdict)){
                                cores[[g]] <- as.vector(rainbow(length(unique(hmdf[ ,g]))))
                                names(cores[[g]]) <- sort(unique(hmdf[ ,g]))
                            } else {
                                ct <- cdict[[colcategories[g]]]
                                ct <- subset(ct, Name %in% hmdf[ , g])
                                cores[[g]] <- as.vector(ct$Hex)
                                names(cores[[g]]) <- as.vector(ct$Name)
                            }
                        } else {
                            #Variable is (or can be made) numeric, but check for variance in the numbers
                            numvals <- as.numeric(hmdf[, g][which(!(hmdf[, g] %in% class_to_ignore))])
                            hmdf[, g] <- as.numeric(hmdf[, g])
                            if ((max(numvals) - min(numvals)) > 0 ){
                                #If values contain class to ignore, make them black else, only span whicte and dark blue
                                cores[[g]] <- colorRamp2(c(min(numvals), max(numvals)), c("#3a8aa7", "#780078"), space = "HSV")
                            } else {
                                #Not enough variance, so make them discrete
                                cores[[g]] <- as.vector(rainbow(length(unique(hmdf[, g]))))
                                names(cores[[g]]) <- sort(unique(hmdf[, g]))
                            }
                        }
                        names(cores)[g] <- colcategories[g]
                    }

                    #Make colour scale for relabund heatmap
                    if (hmasPA == FALSE) {
                        if (scaled == TRUE) {
                            #Scale the matrix
                            sampordernames <- colnames(mathm)
                            mathm <- t(apply(mathm, 1, scale))
                            colnames(mathm) <- sampordernames

                            #This is the colour spectrum we are aiming to span
                            PctHmColours <- c("#1307FC", "#FFFFFF", "#F70C00")
                            #Let us see what the distribution looks like to fit it to the colour spectrum
                            #Transform to PPM
                            #RelabundBreakPts <- c(min(mathm), ((min(mathm) + max(mathm)) / 2), max(mathm))
                            RelabundBreakPts <- c(min(mathm), 0, max(mathm))
                            relabundscalename <- "scaling"
                            RelabundBreakPtsLbls <- round(RelabundBreakPts, 2)
                            HMrelabundBreaks <- RelabundBreakPtsLbls
                        } else {
                            #This is the colour spectrum we are aiming to span
                            PctHmColours <- c("blue4", "blue", "slategray1", "khaki", "orange", "tomato", "red", "magenta2", "magenta4")
                            if (analysis == "LKT"){
                                PctBreakPts <- c(0.0001, 0.001, 0.1, 1, 2.5, 5, 10, 50, 100)
                                RelabundBreakPts <- signif(PctBreakPts, digits = 3)
                                relabundscalename <- "Relative Abundance (%)"
                                RelabundBreakPtsLbls <- as.character(paste0(RelabundBreakPts, "%"))
                                HMrelabundBreaks <- Pct2log2PPM(PctBreakPts)
                            } else {
                                #Let us see what the distribution looks like to fit it to the colour spectrum
                                #Transform to PPM
                                quantprobs <- round(((2:(length(PctHmColours) - 1)) * (1 / length(PctHmColours))), 1)
                                #quantprobs <- c(0.10, 0.20, 0.35, 0.50, 0.85, 0.95, 0.99)
                                nonlogmat <- convert_matrix_log2(mat = mathm, transformation = "from_log2")
                                #Transform to get quantiles
                                countmatdistrib <- apply(nonlogmat, function (x) { quantile(x, probs = quantprobs) }, MARGIN = 2)
                                #take median values of these
                                medpoints <- rowMedians(countmatdistrib)
                                PPMrelabundBreaks <- c(0, medpoints, max(nonlogmat))
                                HMrelabundBreaks <- log2(PPMrelabundBreaks + 1)
                                RelabundBreakPts <- PPMrelabundBreaks
                                relabundscalename <- "Parts Per Million"
                                RelabundBreakPtsLbls <- as.character(paste(RelabundBreakPts, "PPM"))
                            }
                        }
                        relabundheatmapCols <- colorRamp2(HMrelabundBreaks, PctHmColours)
                    } else {
                        #Plot cells within heatmap as present/absent
                        relabundheatmapCols <- colorRamp2(c(0, 1), c("blue4", "red"))
                        RelabundBreakPtsLbls <- c("Absent", paste("Present >=", threshPA, "PPM") )
                        HMrelabundBreaks <- c(0, 1)
                        relabundscalename <- "Pres/Abs"
                    }

                    if (!is.null(textby)){
                        hmdf_txt <- as.character(colData(currobj)[ , which(colnames(colData(currobj)) == textby[1])])
                        ha_column <- HeatmapAnnotation(which = "column", df = hmdf, col = cores, textby = anno_text(hmdf_txt, gp = gpar(fontsize = fontsizex, col = "black")), annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 7, col = "black"))
                    } else {
                        ha_column <- HeatmapAnnotation(which = "column", df = hmdf, col = cores, annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 7, col = "black"))
                    }

                    #Build plot title
                    nsampmsg <- paste0("Number of samples in heatmap = ", ncol(mathm))
                    msgs <- c(maintit, stattit, presetlist$filtermsg, paste(nsampmsg, ffeatmsg, sep = " | "))
                    plotit <- paste0(msgs, collapse = "\n")

                    if (all(c((any(!(c(is.null(presetlist$minl2fc), is.null(presetlist$maxl2fc))))), ("l2fc" %in% colnames(stathm))))){
                        l2fcmsg <- paste(minl2fcmsg, maxl2fcmsg, sep = " | ")
                        plotit <- paste(plotit, l2fcmsg, sep = "\n")
                    }

                    if (any(c("l2fc", "OddsRatio") %in% colnames(stathm))){
                        plotit <- paste(plotit, binary_directionality, sep = "\n")
                    }

                    #Add plot number if there is more than one heatmap matrix.
                    if (length(matlist) > 1){
                        hmcounter <- paste(hm, length(matlist), sep = "/")
                        plotit <- paste(plotit, paste("Heatmap", hmcounter), sep = "\n")
                    }

                    #Switch plot title if doing a dual heatmap with either genome completeness or percentage from contigs
                    #Get genome completeness hmdf if in taxonomic space
                    ht1fs <- 8
                    hm1tit <- plotit

                    if (is.null(secondaryheatmap)) {
                        #Cheap trick
                        secondaryheatmap <- FALSE
                    }

                    if (secondaryheatmap == "GenomeCompleteness") {
                        if ("GenomeCompleteness" %in% names(assays(currobj))){
                            gchmdf <- genomecompletenessdf[rownames(mathm), colnames(mathm)]
                            gchmdf <- as.matrix(gchmdf)
                            gchmdf <- gchmdf * 100
                            #Cap genome completeness to 400% max, to preserve scale across heatmaps
                            gchmdf[which(gchmdf[] > 400)] <- 400
                            if (ncol(mathm) < 51){
                                ht1fs <- 7
                                dualHMtit <- plotit
                                hm1tit <- "Relative Abundance"
                            }
                        }
                    } else if (secondaryheatmap == "PctFromCtgs"){
                        if ("PctFromCtgs" %in% names(assays(currobj))){
                            pctctgsdf <- as.data.frame(assays(currobj)$PctFromCtgs)
                            gchmdf <- pctctgsdf[rownames(mathm), colnames(mathm)]
                            gchmdf <- as.matrix(gchmdf)
                            if (ncol(mathm) < 51){
                                ht1fs <- 7
                                dualHMtit <- plotit
                                hm1tit <- "Relative Abundance"
                            }
                        }
                    }

                    #Include row annotations for pvalues and l2fc if required.
                    if (("l2fc" %in% colnames(stathm))) {
                        #Get l2fc amplitude for the current heatmap
                        l2fcamplitude <- round(stathm$l2fc, 1)
                        maximuml2fctoshow = 15
                        #Take into account that l2fc may contain only infinites
                        l2fcamplitudeculled <- l2fcamplitude
                        l2fcamplitudeculled[which(l2fcamplitudeculled == Inf)] <- maximuml2fctoshow
                        l2fcamplitudeculled[which(l2fcamplitudeculled == -Inf)] <- (maximuml2fctoshow * -1)
                        l2fcamplitudeshifted <- l2fcamplitudeculled + 15 #Offset amplitude to contain only positive values, maintaining their relativity
                    }

                    #Include annotations for Odds Ratio if PA
                    if ("OddsRatio" %in% colnames(stathm)){
                        ORamplitude <- paste0(round(stathm$OddsRatio, 2), paste0("(", paste(round(stathm$OR95lwr, 2), round(stathm$OR95upr, 2), sep = "-"), ")"))
                    }

                    if (showl2fc == TRUE){
                        showl2fc = "text"
                    }

                    #This is not the most eficient way of doing it, but will keep it like this for now as it is working.
                    if (all(c(("pval" %in% colnames(stathm)), showpval))){
                        statannotnonadj <- as.character(signif(stathm[, "padj_none"], 2))
                        statannotadj <- as.character(signif(stathm[, "padj_fdr"], 2))
                        statannot <- paste(statannotnonadj, statannotadj, sep = " | ")

                        #If showing Log2FC then do so
                        if (all(c(("l2fc" %in% colnames(stathm)), (showl2fc %in% c("text", "plot"))))) {
                            if (showl2fc == "text"){
                                #Show pvalue AND l2fc
                                row_ha <- HeatmapAnnotation(which = "row", Pval = anno_text(statannot, gp = gpar(fontsize = fontsizey)), Log2FC = anno_text(l2fcamplitude, gp = gpar(fontsize = fontsizey)), gap = unit(3, "mm"))
                            } else {
                                row_ha <- HeatmapAnnotation(which = "row", Pval = anno_text(statannot, gp = gpar(fontsize = fontsizey)), Log2FC = anno_points(l2fcamplitudeshifted, ylim = c(0, 30), width = unit(0.8, "cm"), axis_param = list(side = "bottom", at = c(0, 15, 30), labels = c("<-15", "0", ">15"), labels_rot = 90)), annotation_name_gp = gpar(fontsize = 6, col = "black"))
                            }
                        } else if ("OddsRatio" %in% colnames(stathm)) {
                            #Show Pval and Odds Ratio
                            row_ha <- rowAnnotation(Pval = anno_text(statannot, gp = gpar(fontsize = fontsizey)), OR = anno_text(ORamplitude, gp = gpar(fontsize = fontsizey)))
                        } else {
                            #Show only pval
                            row_ha <- rowAnnotation(Pval = anno_text(statannot, gp = gpar(fontsize = fontsizey)))
                        }

                    } else {
                        #OK, not showing p-values
                        if (all(c(("l2fc" %in% colnames(stathm)), (showl2fc %in% c("text", "plot"))))) {
                            #If showing Log2FC then do so
                            if (showl2fc == "text"){
                                row_ha <- rowAnnotation(Log2FC = anno_text(l2fcamplitude, gp = gpar(fontsize = fontsizey)))
                            } else {
                                row_ha <- rowAnnotation(Log2FC = anno_points(l2fcamplitudeshifted, ylim = c(0, 30), width = unit(0.8, "cm"), axis_param = list(side = "bottom", at = c(0, 15, 30), labels = c("<-15", "0", ">15"), labels_rot = 90)), annotation_name_gp = gpar(fontsize = 6, col = "black"))
                            }
                        } else if ("OddsRatio" %in% colnames(stathm)) {
                            #Show only Odds Ratio
                            row_ha <- rowAnnotation(OR = anno_text(ORamplitude, gp = gpar(fontsize = fontsizey)))
                        } else if ("Rank" %in% colnames(stathm)){
                            #Include variance rank if Present
                            VarRank <- as.character(stathm$Rank)
                            row_ha <- rowAnnotation(Rank = anno_text(VarRank, gp = gpar(fontsize = fontsizey)))
                        } else {
                            row_ha <- NULL
                        }
                    }

                    #Plot Gram and phyla, if applicable
                    if (all(c(showGram, (analysisname %in% c("LKT", "Species", "Genus", "Family", "Order", "Class"))))){
                        data(Gram)
                        tt <- as.data.frame(rowData(currobj))
                        tt <- tt[rownames(mathm), c(analysisname, "Phylum")]
                        tt <- left_join(tt, Gram, by = "Phylum")
                        phycols <- setNames(as.character(Gram$PhylumColour), as.character(Gram$Phylum))[unique(tt$Phylum)]
                        gramcols <- setNames(as.character(Gram$GramColour), as.character(Gram$Gram))[unique(tt$Gram)]
                        hatax <- rowAnnotation(Phylum = tt$Phylum, Gram = tt$Gram, col = list(Phylum = phycols, Gram = gramcols),  annotation_name_gp = gpar(fontsize = 6, col = "black"), show_legend = TRUE)
                    } else {
                        hatax <- NULL
                    }

                    #Determine column order explicitly if required and draw heatmap
                    if (!is.null(ordercolsby)) {
                        column_order <- order(colData(currobj)[ , which(colnames(colData(currobj)) == ordercolsby)])
                        cluster_column_slices <- FALSE
                    } else {
                        column_order <- NULL
                        cluster_column_slices <- TRUE
                    }

                    #Define if groups will be split
                    if (!(is.null(splitcolsby))){

                        splitcol <- as.character(colData(currobj)[ , which(colnames(colData(currobj)) == splitcolsby)])
                        colgroups <- sort(unique(splitcol))
                        column_split <- factor(splitcol, levels = colgroups)

                        ht1 <- Heatmap(mathm, name = relabundscalename, cluster_columns = cluster_column_slices, column_split = column_split, column_order = column_order, cluster_column_slices = cluster_column_slices, show_column_dend = TRUE, column_dend_side = "top", column_gap = unit(3, "mm"), column_title = hm1tit, column_title_gp = gpar(fontsize = ht1fs), top_annotation = ha_column, col = relabundheatmapCols, column_names_gp = gpar(fontsize = fontsizex), right_annotation = row_ha, left_annotation = hatax, cluster_rows = cluster_rows, show_row_dend = FALSE, row_names_side = "left", row_names_gp = gpar(fontsize = fontsizey, col = rowlblcol),  row_names_max_width = unit(6, "cm"), heatmap_legend_param = list(direction = "horizontal", title = relabundscalename, labels = RelabundBreakPtsLbls, at = HMrelabundBreaks, title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 6)))

                    } else {

                        column_split <- NULL
                        if (any(cluster_samples_per_heatmap, !cluster_rows)){

                            ht1 <- Heatmap(mathm, name = relabundscalename, column_title = hm1tit, column_title_gp = gpar(fontsize = ht1fs), top_annotation = ha_column, col = relabundheatmapCols, column_names_gp = gpar(fontsize = fontsizex), column_dend_height = unit(5, "mm"), right_annotation = row_ha, left_annotation = hatax, cluster_rows = cluster_rows, show_row_dend = FALSE, row_names_side = "left", row_names_gp = gpar(fontsize = fontsizey, col = rowlblcol), heatmap_legend_param = list(direction = "horizontal", title = relabundscalename, labels = RelabundBreakPtsLbls, at = HMrelabundBreaks, title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 4)), row_names_max_width = unit(6, "cm"))

                        } else {

                            #Coerce column order to the order obtained using the full countmatrix
                            ht1 <- Heatmap(mathm, name = relabundscalename, column_title = hm1tit, column_title_gp = gpar(fontsize = ht1fs), top_annotation = ha_column, col = relabundheatmapCols, column_names_gp = gpar(fontsize = fontsizex), cluster_columns = fullheatmap_column_dend, column_dend_height = unit(5, "mm"), right_annotation = row_ha, left_annotation = hatax, cluster_rows = cluster_rows, show_row_dend = FALSE, row_names_side = "left", row_names_gp = gpar(fontsize = fontsizey, col = rowlblcol), heatmap_legend_param = list(direction = "horizontal", title = relabundscalename, labels = RelabundBreakPtsLbls, at = HMrelabundBreaks, title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 4)), row_names_max_width = unit(6, "cm"))

                        }
                    }

                    #Make a genome completeness heatmap if in taxonomic space
                    if (all(c((secondaryheatmap %in% c("GenomeCompleteness", "PctFromCtgs")), (any(c(("GenomeCompleteness" %in% names(assays(currobj))), ("PctFromCtgs" %in% names(assays(currobj))))))))) {

                        if (secondaryheatmap == "GenomeCompleteness"){
                            #Draw heatmap with completeness
                            GCheatmapCols <- colorRamp2(c(0, 100, 200, 300, 400), c("white", "forestgreen", "blue", "firebrick1", "black"))
                            GCha_column <- HeatmapAnnotation(df = hmdf, col = cores, show_annotation_name = FALSE)
                            RelabundRowOrder <- row_order(ht1)
                            ht2 <- Heatmap(gchmdf, name = "GenComp", column_split = column_split, column_title = "% Genome completeness", column_title_gp = gpar(fontsize = ht1fs), top_annotation = GCha_column, col = GCheatmapCols, column_names_gp = gpar(fontsize = fontsizex), right_annotation = NULL, left_annotation = hatax, cluster_rows = FALSE, column_order = unlist(column_order(ht1)), row_order = RelabundRowOrder, show_row_dend = FALSE, row_names_side = "left", row_names_gp = gpar(fontsize = fontsizey, col = rowlblcol), heatmap_legend_param = list(direction = "horizontal", title = "% GenComp", labels = c("0%", "100%", "200%", "300%", "> 400%"), title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 4)), row_names_max_width = unit(6, "cm"))
                        } else if (secondaryheatmap == "PctFromCtgs"){
                            #Draw heatmap with percentage from contigs
                            GCheatmapCols <- colorRamp2(c(0, 100), c("white", "midnightblue"))
                            GCha_column <- HeatmapAnnotation(df = hmdf, col = cores, show_annotation_name = FALSE)
                            RelabundRowOrder <- row_order(ht1)
                            ht2 <- Heatmap(gchmdf, name = "PctFromCtgs", column_split = column_split, column_title = "% Taxonomic info from Contigs", column_title_gp = gpar(fontsize = ht1fs), top_annotation = GCha_column, col = GCheatmapCols, column_names_gp = gpar(fontsize = fontsizex), right_annotation = NULL, left_annotation = hatax, cluster_rows = FALSE, column_order = unlist(column_order(ht1)), row_order = RelabundRowOrder, show_row_dend = FALSE, row_names_side = "left", row_names_gp = gpar(fontsize = fontsizey, col = rowlblcol), heatmap_legend_param = list(direction = "horizontal", title = "PctFromCtgs", title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 4)), row_names_max_width = unit(6, "cm"))
                        }

                        #Plot heatmaps side by side if there are 50 samples or less. Else plot one on each page.
                        if (ncol(mathm) < 51){
                            ht_list = ht1 + ht2
                            draw(ht_list, heatmap_legend_side = "bottom", annotation_legend_side = "right", ht_gap = unit(0.2, "cm"),padding = unit(c(2, 2, 2, 2), "mm"), column_title = dualHMtit, column_title_gp = gpar(fontsize = ht1fs))
                            drawseparateGChm <- FALSE
                        } else {
                            par(oma = c(10, 7, 3, 10) + 0.1, xpd = TRUE)
                            draw(ht1, heatmap_legend_side = "bottom", annotation_legend_side = "right", padding = unit(c(2, 2, 2, 2), "mm"))
                            drawseparateGChm <- TRUE
                        }

                    } else {
                        #Draw the heatmap
                        par(oma = c(10, 7, 3, 10) + 0.1, xpd = TRUE)
                        draw(ht1, heatmap_legend_side = "bottom", annotation_legend_side = "right", padding = unit(c(2, 2, 2, 2), "mm"))
                        drawseparateGChm <- FALSE
                    }

                    #Print what the column annotations are
                    #for(an in colnames(hmdf)) {
                    #    decorate_annotation(an, {
                    #        grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right", gp = gpar(fontsize = 10, col = "black"))
                    #    } )
                    #}
                    if (all(c(("pval" %in% colnames(stathm)), showpval))) {
                        decorate_annotation("Pval", {
                            grid.text("p | FDR p", y = unit(1, "npc") + unit(5, "mm"), just = "bottom", gp = gpar(fontsize = 5, col = "black"))
                        })
                    }
                    if (all(c(("l2fc" %in% colnames(stathm)), (showl2fc %in% c("text", "plot", TRUE))))) {
                        decorate_annotation("Log2FC", {
                            grid.text("Log2FC", y = unit(1, "npc") + unit(5, "mm"), just = "bottom", gp = gpar(fontsize = 5, col = "black"))
                        })
                    }
                    if ("OddsRatio" %in% colnames(stathm)) {
                        decorate_annotation("OR", {
                            grid.text("Odds Ratio", y = unit(1, "npc") + unit(5, "mm"), just = "bottom", gp = gpar(fontsize = 5, col = "black"))
                        })
                    }
                    if ("Rank" %in% colnames(stathm)) {
                        decorate_annotation("Rank", {
                            grid.text("Rank", y = unit(1, "npc") + unit(5, "mm"), just = "bottom", gp = gpar(fontsize = 5, col = "black"))
                        })
                    }

                    #gvec[[n]]<-recordPlot()
                    n <- n + 1

                    if (drawseparateGChm == TRUE){
                        draw(ht2, heatmap_legend_side = "bottom", annotation_legend_side = "right", padding = unit(c(2, 2, 2, 2), "mm"))
                        n <- n + 1
                    }

                    if(n > 100){
                        stop("There are too many combinations. I think you have had enough plots. I am stopping here.")
                    }
                }
            } else {
                #There is no valid matrix, so print out the conditions which led the matrix to be empty.
                plotit <- c(stattit, paste(presetlist$filtermsg, sep=", "))
                plot.new()
                grid.table(c(maintit, "No features fulfilling", plotit), rows = NULL, cols = NULL, theme = ttheme_default(base_size = 10))
                #gvec[[n]]<-recordPlot()
                n <- n + 1
            } #End conditional if there are heatmaps to plot from stats matrix

        } else {
            if (ncol(colData(currobj)) < 2){
                #There is only one sample in subset, so nothing was done.
                flog.info(paste(as.character(colnames(colData(currobj))[samplesToKeep]), "was the only sample within", subsetname, "subset. Must have at least two samples for a heatmap."))
            }
            if (nrow(colData(currobj)) < 2){
                #There is only one sample in subset, so nothing was done.
                flog.info(paste(as.character(rownames(assays(currobj)$BaseCounts)), "was the only feature within", subsetname, "subset. Must have at least two features for a heatmap."))
            }
        } #End conditional if there is more than a two samples and or two features in subset

    } #End subsetby loop

    #Redefine stats list as ones only containing data
    svec <- svec[sapply(svec, function(x){ !(is.null(x)) } )]

    if (returnstats == TRUE){
        return(svec)
    } else {
        flog.info("Heatmaps generated.")
    }
}
