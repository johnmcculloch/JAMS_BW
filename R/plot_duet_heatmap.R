#' plot_duet_heatmap(mgseqobj = NULL, glomby = NULL, colcategories = NULL, compareby = NULL, adjustpval = NULL, showl2fc = TRUE, showpval = TRUE, stattype = "spearman", subsetby = NULL, maxnumfeatallowed = 10000, minabscorrcoeff = NULL, ntopvar = NULL, featmaxatleastPPM = 0, featcutoff = c(0, 0), applyfilters = NULL, featuresToKeep = NULL, samplesToKeep = NULL, genomecompleteness = NULL, list.data = NULL, showGram = TRUE, showphylum = TRUE, addtit = NULL, mgSeqnorm = FALSE, cdict = NULL, ignoreunclassified = TRUE, class_to_ignore = NULL)
#'
#' Plots correlation heatmaps annotated by the metadata or a correlelogram of features
#' DISCLAIMER: This function is experimental
#' @export

plot_duet_heatmap <- function(mgseqobj = NULL, glomby = NULL, colcategories = NULL, compareby = NULL, adjustpval = NULL, showl2fc = TRUE, showpval = TRUE, stattype = "spearman", subsetby = NULL, maxnumfeatallowed = 10000, minabscorrcoeff = NULL, ntopvar = NULL, featmaxatleastPPM = 0, featcutoff = c(0, 0), applyfilters = NULL, featuresToKeep = NULL, samplesToKeep = NULL, genomecompleteness = NULL, list.data = NULL, showGram = TRUE, showphylum = TRUE, addtit = NULL, mgSeqnorm = FALSE, cdict = NULL, ignoreunclassified = TRUE, class_to_ignore = NULL) {


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
                maintit <- paste("Feature Duet Heatmap", analysisname, paste("within", subset_points[sp]), sep = " | ")
            } else {
                maintit <- paste("Feature Duet Heatmap", analysisname, sep = " | ")
            }
            if (!is.null(addtit)) {
                maintit <- paste(addtit, maintit, sep = "\n")
            }

            #Get counts matrix
            countmat <- MRcounts(currobj, norm = FALSE, log = TRUE)

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
                corrmatstats <- calculate_matrix_stats(countmatrix = countmat, uselog = FALSE, statsonlog = FALSE, stattype = stattype, classesvector = NULL)

                if (!is.null(minabscorrcoeff)){
                    print(paste("Eliminating features which do not correlate with other features with a coefficient of at least", minabscorrcoeff))
                    corrmatstats <- filter_correlations(corrmat= corrmatstats, mincorrelcoeff = minabscorrcoeff)
                    minabscorrcoeffmsg <- paste("Largest correlation coefficient at least", minabscorrcoeff)
                } else {
                    minabscorrcoeffmsg <- NULL
                }

                #Plot heatmap
                #Set scale
                #This is the colour spectrum we are aiming to span
                CorrHmColours <- c("blue4", "lightgoldenrodyellow", "red1")
                heatmapCols <- colorRamp2(c(-1, 0, 1), CorrHmColours)

                fontsizey <- max(1, (min(5, round((((-1 / 300) * (nrow(corrmatstats))) + 1) * 3, 0))))

                #Add genome completeness info if LKT
                if (analysis == "LKT"){
                    genomecompletenessstats <- as.matrix(genomecompletenessdf[rownames(corrmatstats), ])
                    gcl <- lapply(1:nrow(genomecompletenessstats), function (x){ (as.numeric(genomecompletenessstats[x, ][which(genomecompletenessstats[x, ] != 0)])) * 100 })

                    data(Gram)
                    #Get Phyla
                    if (analysisname %in% c("LKT", "Species", "Genus", "Family", "Order", "Class")){
                        tt <- fData(currobj)
                        tt <- tt[rownames(corrmatstats), c(analysisname, "Phylum")]
                        Gram$Kingdom <- NULL
                        tt <- left_join(tt, Gram)
                        tt$Gram[which(!(tt$Gram %in% c("positive", "negative")))] <- "not_sure"
                        phcol <- colorRampPalette((brewer.pal(9, "Set1")))(length(unique(tt$Phylum)))
                        names(phcol) <- unique(tt$Phylum)
                        phcol[which(names(phcol) == "p__Unclassified")] <- "#000000"
                        phcol <- phcol[!duplicated(phcol)]
                    }

                    #ha1 <- rowAnnotation(Pct_Genome_Compl = anno_boxplot(gcl, width = unit(4, "cm"), pch = 20, size = unit(1, "mm"),  axis_param = list(labels_rot = 90)), Gram = tt$Gram, Phylum = tt$Phylum, col = list(Gram = c("positive" = "#7D00C4", "negative" = "#FC0345", "not_sure" = "#B8B8B8"), Phylum = phcol),  annotation_name_gp = gpar(fontsize = 6, col = "black"))
                    ha1 <- NULL

                    ha2 <- HeatmapAnnotation(Phylum = tt$Phylum, Gram = tt$Gram, col = list(Phylum = phcol, Gram = c("positive" = "#7D00C4", "negative" = "#FC0345", "not_sure" = "#B8B8B8")), annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 6, col = "black"), show_legend = TRUE, annotation_legend_param = list(nrow = 3, direction = "horizontal"))

                } else {
                    ha1 <- NULL
                    ha2 <- NULL
                }

                #Build plot title
                stattit <- paste("Correlation measure =", stattype)
                plotit <- paste(maintit, stattit, cutoffmsg, minPPMmsg, ntopvarmsg, minabscorrcoeffmsg, sep = "\n")

                ht1 <- Heatmap(corrmatstats, name = paste(stattype, "correlation coefficient"), column_title = plotit, column_title_gp = gpar(fontsize = 8), col = heatmapCols, column_dend_height = unit(5, "mm"), cluster_rows = TRUE, show_row_dend = FALSE, row_names_gp = gpar(fontsize = fontsizey), column_names_gp = gpar(fontsize = fontsizey), heatmap_legend_param = list(direction = "horizontal", legend_width = unit(6, "cm"), title = paste(stattype, "correlation coefficient"), labels = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1), at = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1), title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 6)), left_annotation = ha1, bottom_annotation = ha2)

                CorrelRowOrder <- row_order(ht1)
                CorrelRowNames <- rownames(corrmatstats)[CorrelRowOrder]

                #Get relabund matrix with Correlation matrix order
                #mathm <- countmat[CorrelRowNames, ]
                mathm <- countmat[rownames(corrmatstats), ]

                fontsizex <- as.numeric(unlist(round((((-1 / 150) * (ncol(mathm))) + 1) * 5, 0)))

                cl <- pData(currobj)[ , which(colnames(pData(currobj)) == compareby)]
                discretenames <- sort(unique(cl))

                stathm <- calculate_matrix_stats(countmatrix = mathm, uselog = TRUE, statsonlog = TRUE, stattype = "auto", classesvector = cl, invertbinaryorder = FALSE)
                stathm <- stathm[rownames(mathm), ]
                #rownames(mathm) <- strtrim(rownames(mathm), 60)


                if ("l2fc" %in% colnames(stathm)){
                    stathm$Colour <- ifelse(stathm$l2fc < 0, "#900000", "#000000")
                    statmsg <- paste("MWW", compareby, sep="_")
                } else if ("oddsRatio" %in% colnames(stathm)){
                    stathm$Colour <- ifelse(stathm$oddsRatio < 1, "#900000", "#000000")
                    statmsg <- paste("Fisher", compareby, sep="_")
                } else {
                    stathm$Colour <- rep("#000000", nrow(stathm))

                    if (stathm$Method[1] == "permanova"){
                        statmsg <- paste("PERMANOVA", compareby, sep="_")
                    } else if (stathm$Method[1] == "anova"){
                        statmsg <- paste("ANOVA", compareby, sep="_")
                    } else if (stathm$Method[1] == "variance"){
                        statmsg <- paste("Var", compareby, sep="_")
                    }
                }

                rowlblcol <- stathm$Colour


                l2fcmeaning <- paste("Positive l2fc means increased in", discretenames[1])

                if (("pval" %in% colnames(stathm)) && adjustpval == "auto"){
                    propsigadj <- length(which(stathm$padj_fdr < 0.05)) / length(stathm$padj_fdr)
                    propsignonadj <- length(which(stathm$pval < 0.05)) / length(stathm$pval)
                    fracsigadj <- propsigadj / propsignonadj
                    if ((!is.na(fracsigadj)) && (fracsigadj > 0.2)){
                        adjustpval <- TRUE
                    } else {
                        adjustpval <- FALSE
                    }
                }

                if (("pval" %in% colnames(stathm)) && adjustpval != TRUE){
                    sigmeas <- "pval"
                } else if (("pval" %in% colnames(stathm)) && adjustpval == TRUE){
                    sigmeas <- "padj_fdr"
                }


                if (stathm$Method[1] == "MannWhitneyWilcoxon") {
                    stattithm <- paste(sigmeas, "of difference between", compareby, "using MannWhitneyWilcoxon")
                } else if (stathm$Method[1] == "permanova"){
                    stattithm <- paste(sigmeas, "of difference between", compareby, "using PERMANOVA")
                } else if (stathm$Method[1] == "anova"){
                    stattithm <- paste(sigmeas, "of difference between", compareby, "using ANOVA")
                } else if (stathm$Method[1] == "fisher"){
                    stattithm <- paste(sigmeas, "of Present/Absent between", compareby, "using Fishers test")
                }

                #Name stats in svec
                stattitlerelabund <- paste(analysisname, statmsg, subsetname, sep = "_")

                hmdf <- as.data.frame(matrix(data = 0, nrow = nrow(pData(currobj)), ncol = length(colcategories)))
                cores <- vector("list", length = length(colcategories))
                for (g in 1:length(colcategories)){
                    hmdf[ , g] <- pData(currobj)[ , which(colnames(pData(currobj)) == colcategories[g])]
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

                hmasPA <- FALSE
                #Make colour scale for relabund heatmap
                if ( hmasPA == FALSE ) {
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
                        nonlogmat <- round(((2 ^ (mathm)) - 1), 0)
                        #Transform tget quantiles
                        countmatdistrib <- apply(nonlogmat, function (x) {quantile(x, probs = quantprobs) }, MARGIN = 2)
                        #take median values of these
                        medpoints <- rowMedians(countmatdistrib)
                        PPMrelabundBreaks <- c(0, medpoints, max(nonlogmat))
                        HMrelabundBreaks <- log2(PPMrelabundBreaks + 1)
                        RelabundBreakPts <- PPMrelabundBreaks
                        relabundscalename <- "Parts Per Million"
                        RelabundBreakPtsLbls <- as.character(paste(RelabundBreakPts, "PPM"))
                    }
                    relabundheatmapCols <- colorRamp2(HMrelabundBreaks, PctHmColours)
                } else {
                    relabundheatmapCols <- colorRamp2(c(0, 1), c("blue4", "red"))
                    relabundscalename <- "Pres/Abs"
                }


                ha_column_relabund <- HeatmapAnnotation(df = hmdf, col = cores, annotation_name_side = "right", annotation_name_gp = gpar(fontsize = 7, col = "black"))


                #Build plot title
                nsampmsg <- paste0("Number of samples in heatmap = ", ncol(mathm))
                msgs <- c("Relative Abundance", stattithm, nsampmsg)
                plotitrelabund <- paste0(msgs, collapse = "\n")
                hm2tit <- plotitrelabund

                #This is not the most eficient way of doing it, but will keep it like this for now as it is working.
                if (all(c(("pval" %in% colnames(stathm)), showpval))){
                    statannot <- as.character(signif(stathm[, sigmeas], 2))
                    #If showing Log2FC then do so
                    if (all(c(("l2fc" %in% colnames(stathm)), (showl2fc %in% c("text", "plot"))))) {
                        if (showl2fc == "text"){
                            #Show pvalue AND l2fc
                            row_ha <- HeatmapAnnotation(which = "row", Pval = anno_text(statannot, gp = gpar(fontsize = fontsizey)), Log2FC = anno_text(l2fcamplitude, gp = gpar(fontsize = fontsizey)), gap = unit(3, "mm"))
                        } else {
                            row_ha <- HeatmapAnnotation(which = "row", Pval = anno_text(statannot, gp = gpar(fontsize = fontsizey)), Log2FC = anno_points(l2fcamplitudeshifted, ylim = c(0, 30), width = unit(0.8, "cm"), axis_param = list(side = "bottom", at = c(0, 15, 30), labels = c("<-15", "0", ">15"), labels_rot = 90)), annotation_name_gp = gpar(fontsize = 6, col = "black"))
                        }
                    } else if ("oddsRatio" %in% colnames(stathm)) {
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
                    } else if ("oddsRatio" %in% colnames(stathm)) {
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

                ht2fs <- 6

                ht2 <- Heatmap(mathm, row_order = 1:nrow(mathm), name = relabundscalename, column_title = hm2tit, column_title_gp = gpar(fontsize = ht2fs), top_annotation = ha_column_relabund, col = relabundheatmapCols, column_names_gp = gpar(fontsize = fontsizex), column_dend_height = unit(5, "mm"), right_annotation = row_ha, cluster_rows = FALSE, show_row_dend = FALSE, row_names_side = "right", row_names_gp = gpar(fontsize = fontsizey, col = rowlblcol), heatmap_legend_param = list(direction = "horizontal", title = relabundscalename, labels = RelabundBreakPtsLbls, at = HMrelabundBreaks, title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize = 6)), row_names_max_width = unit(6, "cm"))

                draw(ht1 + ht2, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = FALSE, padding = unit(c(2, 20, 2, 2), "mm"))

                #Now, make a relabund heatmap but with the row orders for

            } #End conditional of going ahead and doing correlations if there arent too many features
        } #End conditional if there are any features left over after filtering
    } #End for loop for plotting within each subset point

    #Redefine stats list as ones only containing data
    #svec <- svec[sapply(svec, function(x){ !(is.null(x)) } )]

    #if (returnstats == TRUE){
    #    return(svec)
    #} else {
    message("Heatmaps generated.")
    #}
}

