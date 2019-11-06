#' plot_correlation_heatmap(mgseqobj = NULL, glomby = NULL, heatpalette = "smart", stattype = "spearman", cluster_rows = FALSE, subsetby = NULL, maxnumfeatallowed = 500, featmaxatleastPPM = 0, featcutoff = c(0, 0), applyfilters = NULL, featuresToKeep = NULL, samplesToKeep = NULL, genomecompleteness = NULL, list.data = NULL, showGram = FALSE, showphylum = FALSE, addtit = NULL, mgSeqnorm = FALSE, cdict = NULL, ignoreunclassified = TRUE, returnstats = FALSE, ...)
#'
#' Plots correlation heatmaps annotated by the metadata or a correlelogram of features
#' @export

plot_correlation_heatmap <- function(mgseqobj = NULL, glomby = NULL, stattype = "spearman", subsetby = NULL, maxnumfeatallowed = 10000, minabscorrcoeff = NULL, ntopvar = NULL, featmaxatleastPPM = 0, featcutoff = c(0, 0), applyfilters = NULL, featuresToKeep = NULL, samplesToKeep = NULL, genomecompleteness = NULL, list.data = NULL, showGram = TRUE, showphylum = TRUE, addtit = NULL, mgSeqnorm = FALSE, cdict = NULL, ignoreunclassified = TRUE, class_to_ignore = NULL) {

    #Define other functions
    filter_correlations <- function(corrmat = NULL, mincorrelcoeff = NULL){

        if(nrow(corrmat) != ncol(corrmat)){
            stop("Correlation matrix must have equal numbers of rows and columns.")
        }

        featsIwant <- NULL

        for (rw in 1:nrow(corrmat)){
            featint <- rownames(corrmat)[rw]
            #print(paste("Checking:", featint))
            correlations <- corrmat[which(rownames(corrmat) != featint), featint]

            if(max(abs(correlations)) >= mincorrelcoeff){
                feat <- featint
            } else {
                feat <- NULL
            }

            featsIwant <- append(featsIwant, feat)

        }

        corrmat <- corrmat[featsIwant, featsIwant]

        return(corrmat)
    }

    #Get appropriate object to work with
    obj <- mgseqobj

    #Exclude samples and features if specified
    if (!(is.null(samplesToKeep))){
        obj <- obj[, samplesToKeep]
    }

    if (!(is.null(featuresToKeep))){
        obj <- obj[featuresToKeep, ]
    }

    analysis <- attr(obj, "analysis")
    analysisname <- analysis

    if (!(is.null(glomby))){
        obj <- agglomerate_features(mgseqobj = obj, glomby = glomby)
        if (analysis != "LKT"){
            analysisname <- attr(obj, "analysis")
        } else {
            analysisname <- glomby
        }
    }

    if (!is.null(applyfilters)){
        if (applyfilters == "stringent"){
            if (analysis == "LKT"){
                featcutoff <- c(2000, 15)
                genomecompleteness <- 0.3
                #minabscorrcoeff <- 0.8
            } else {
                featcutoff <- c(50, 15)
                genomecompleteness <- NULL
                #minabscorrcoeff <- 0.8
            }
        } else if (applyfilters == "moderate"){
            if (analysis == "LKT"){
                featcutoff <- c(500, 10)
                genomecompleteness <- 0.1
                #minabscorrcoeff <- 0.5
            } else {
                featcutoff <- c(10, 5)
                genomecompleteness <- NULL
                #minabscorrcoeff <- 0.5
            }
        }
    }

    if ((analysis != "LKT") && (!(is.null(genomecompleteness)))){
        warning("Genome completeness only makes sense for taxa. Please choose a taxonomic (non functional) analysis.")
        genomecompleteness <- NULL
    }

    if ((!(is.null(genomecompleteness))) && (is.null(list.data))){
        stop("Genome completeness can only be obtained by supplying a list.data object.")
    }

    if (!(is.null(subsetby))){
        subset_points <- sort(unique((pData(obj)[, which(colnames(pData(obj)) == subsetby)])))
    } else {
        subset_points <- "none"
    }

    #Initialize Stats and Graph Vector lists
    #svec <- NULL
    #svec <- vector("list", length = 1000)
    #s <- 1
    #n <- 1

    numfeats <- nrow(MRcounts(obj))

    #subset by metadata column
    for (sp in 1:length(subset_points)) {
        if (!(is.null(subsetby))){
            samplesToKeep <- which((pData(obj)[, which(colnames(pData(obj)) == subsetby)]) == subset_points[sp])
            print(paste("Plotting within", subset_points[sp]))
            subsetname <- subset_points[sp]
        } else {
            samplesToKeep = rownames(pData(obj))
            subsetname <- "no_sub"
        }

        #There must be at least two samples for a heatmap and at least two features
        if ((length(samplesToKeep) > 1) && (numfeats > 1)){
            #Discard features which do not match certain criteria
            if (!(is.null(featcutoff))){
                thresholdPPM <- featcutoff[1]
                sampcutoffpct <- min(featcutoff[2], 100)
                cutoffmsg <- paste("Feature must be >", thresholdPPM, "PPM in at least ", sampcutoffpct, "% of samples", sep = "")
            } else {
                cutoffmsg <- "Feature must be > 0 PPM in at least 0% of samples"
                featcutoff <- c(0, 0)
            }

            if (all(c((!(is.null(featmaxatleastPPM))), (featmaxatleastPPM > 0)))) {
                minPPMmsg <- paste("Highest feature must be >", featmaxatleastPPM, "PPM", sep = " ")
            } else {
                minPPMmsg <- NULL
                featmaxatleastPPM <- 0
            }

            currobj <- filter_experiment(mgseqobj = obj, featmaxatleastPPM = featmaxatleastPPM, featcutoff = featcutoff, samplesToKeep = samplesToKeep, asPA = FALSE, asPPM = TRUE, mgSeqnorm = mgSeqnorm)

            #Compose an appropriate title for the plot
            if (length(unique(subset_points)) > 1){
                maintit <- paste("Feature Correlation Heatmap", analysisname, paste("within", subset_points[sp]), sep = " | ")
            } else {
                maintit <- paste("Feature Correlation Heatmap", analysisname, sep = " | ")
            }
            if (!is.null(addtit)) {
                maintit <- paste(addtit, maintit, sep = "\n")
            }

            #Get counts matrix
            countmat <- MRcounts(currobj, norm = FALSE, log = FALSE)

            #Protect against rows with empty data
            rowsToKeep <- which(rowSums(countmat) > 0 & rownames(countmat) != "")
            countmat <- countmat[rowsToKeep, ]

            if (ignoreunclassified == TRUE){
                dunno <- c(paste(analysis, "none", sep = "_"), "LKT__d__Unclassified", "LKT__Unclassified")
                rowsToKeep <- which(!(rownames(countmat) %in% dunno))
                countmat <- countmat[rowsToKeep, ]
            }

            if (!is.null(ntopvar)){
                ntop <- min(ntopvar, nrow(countmat))
                featsds <- rowSds(countmat)
                featIndices <- names(featsds[order(featsds, decreasing = TRUE)[1:ntop]])
                countmat <- countmat[featIndices, ]
                ntopvarmsg <- paste("Top", ntop, "most variant features across samples")
            } else {
                ntopvarmsg <- NULL
            }

            #Rename rows to include description if not taxonomic data
            if (analysis != "LKT"){
                feattable <- fData(currobj)
                feattable$Feature <- paste(feattable$Accession, feattable$Description, sep = "-")
                rownames(countmat) <- feattable$Feature[match(rownames(countmat), feattable$Accession)]
            } else {
                #get genome completeness for taxonomic objects
                genomecompletenessdf <- get_genome_completeness(pheno = pData(currobj), list.data = list.data)
            }
            matrixSamples <- colnames(countmat)
            matrixRows <- rownames(countmat)

            #Discard taxa below required level of completeness
            if (!(is.null(genomecompleteness))){
                print(paste("Genome completeness must be", genomecompleteness, "in at least one sample"))
                featuresToKeep2 <- rownames(genomecompletenessdf)[which(rowMax(as.matrix(genomecompletenessdf)) >= genomecompleteness)]
                countmat <- countmat[(rownames(countmat)[(rownames(countmat) %in% featuresToKeep2)]), ]
                completenessmsg <- paste("Genome completeness >", genomecompleteness)
            }

            docorrelations <- TRUE
            if (!(is.null(maxnumfeatallowed))) {
                if (nrow(countmat) > maxnumfeatallowed){
                    print(paste("There are", nrow(countmat), "features pairwise correlate, which is more than",  maxnumfeatallowed, "allowed. This would entail", (nrow(countmat) ^ 2) , "comparisons. If you are sure you want that many, set maxnumfeatallowed to a higher value."))
                    docorrelations <- FALSE
                }
            }

            if (docorrelations == TRUE){
                #Calculate matrix stats and get new matrix with correlations.
                matstats <- calculate_matrix_stats(countmatrix = countmat, uselog = FALSE, statsonlog = FALSE, stattype = stattype, classesvector = NULL)

                if (!is.null(minabscorrcoeff)){
                    print(paste("Eliminating features which do not correlate with other features with a coefficient of at least", minabscorrcoeff))
                    matstats <- filter_correlations(corrmat= matstats, mincorrelcoeff = minabscorrcoeff)
                    minabscorrcoeffmsg <- paste("Largest correlation coefficient at least", minabscorrcoeff)
                } else {
                    minabscorrcoeffmsg <- NULL
                }

                #Plot heatmap
                #Set scale
                #This is the colour spectrum we are aiming to span
                CorrHmColours <- c("blue4", "lightgoldenrodyellow", "red1")
                heatmapCols <- colorRamp2(c(-1, 0, 1), CorrHmColours)

                fontsizey <- min(6, round((((-1 / 300) * (nrow(matstats))) + 1) * 4, 0))

                #Add genome completeness info if LKT
                if (analysis == "LKT"){
                    genomecompletenessstats <- as.matrix(genomecompletenessdf[rownames(matstats), ])
                    gcl <- lapply(1:nrow(genomecompletenessstats), function (x){ (as.numeric(genomecompletenessstats[x, ][which(genomecompletenessstats[x, ] != 0)])) * 100 })

                    data(Gram)
                    #Get Phyla
                    if (analysisname %in% c("LKT", "Species", "Genus", "Family", "Order", "Class")){
                        tt <- fData(currobj)
                        tt <- tt[rownames(matstats), c(analysisname, "Phylum")]
                        Gram$Kingdom <- NULL
                        tt <- left_join(tt, Gram)
                        tt$Gram[which(!(tt$Gram %in% c("positive", "negative")))] <- "not_sure"
                        phcol <- colorRampPalette((brewer.pal(9, "Set1")))(length(unique(tt$Phylum)))
                        names(phcol) <- unique(tt$Phylum)
                        phcol[which(names(phcol) == "p__Unclassified")] <- "#000000"
                        phcol <- phcol[!duplicated(phcol)]
                    }

                    ha1 <- rowAnnotation(Pct_Genome_Compl = anno_boxplot(gcl, width = unit(4, "cm"), pch = 20, size = unit(1, "mm"),  axis_param = list(labels_rot = 90)), Gram = tt$Gram, Phylum = tt$Phylum, col = list(Gram = c("positive" = "#7D00C4", "negative" = "#FC0345", "not_sure" = "#B8B8B8"), Phylum = phcol),  annotation_name_gp = gpar(fontsize = 6, col = "black"))

                    ha2 <- HeatmapAnnotation(Phylum = tt$Phylum, Gram = tt$Gram, col = list(Phylum = phcol, Gram = c("positive" = "#7D00C4", "negative" = "#FC0345", "not_sure" = "#B8B8B8")),  annotation_name_gp = gpar(fontsize = 6, col = "black"), show_legend = FALSE)

                } else {
                    ha1 <- NULL
                    ha2 <- NULL
                }

                #Build plot title
                stattit <- paste("Correlation measure =", stattype)
                plotit <- paste(maintit, stattit, cutoffmsg, minPPMmsg, ntopvarmsg, minabscorrcoeffmsg, sep = "\n")

                ht1 <- Heatmap(matstats, name = paste(stattype, "correlation coefficient"), column_title = plotit, column_title_gp = gpar(fontsize = 10), col = heatmapCols, column_dend_height = unit(5, "mm"), cluster_rows = TRUE, show_row_dend = FALSE, row_names_gp = gpar(fontsize = fontsizey), column_names_gp = gpar(fontsize = fontsizey), heatmap_legend_param = list(direction = "horizontal", legend_width = unit(6, "cm"), title = paste(stattype, "correlation coefficient"), labels = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1), at = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1), title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 6)), left_annotation = ha1, bottom_annotation = ha2)

                draw(ht1, heatmap_legend_side = "bottom", annotation_legend_side = "left", padding = unit(c(2, 20, 2, 2), "mm"))
            } #End conditional of going ahead and doing correlations if there arent too many features
        } #End conditional if there are any features left over after filtering
    } #End for loop for plotting within each subset point

    #Redefine stats list as ones only containing data
    #svec <- svec[sapply(svec, function(x){ !(is.null(x)) } )]

    #if (returnstats == TRUE){
    #    return(svec)
    #} else {
        return(print("Heatmaps generated."))
    #}
}
