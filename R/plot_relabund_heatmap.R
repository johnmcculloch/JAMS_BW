#' plot_relabund_heatmap(mgseqobj = NULL, glomby = NULL, heatpalette = "diverging", hmtype = NULL, hmasPA = FALSE, compareby = NULL, invertbinaryorder = FALSE, ntop = NULL, ordercolsby = NULL, colcategories = NULL, cluster_rows = FALSE, subsetby = NULL, applyfilters = NULL, featmaxatleastPPM = 0, featcutoff = c(0, 0), samplesToKeep = NULL, featuresToKeep = NULL, adjustpval = FALSE, showonlypbelow = NULL, showpval = TRUE, showl2fc = TRUE, maxl2fc = NULL, minl2fc = NULL, genomecompleteness = NULL, list.data = NULL, addtit = NULL,  scaled = FALSE, mgSeqnorm = FALSE, cdict = NULL, maxnumheatmaps = NULL, numthreads = 4, nperm = 99, statsonlog = TRUE, ignoreunclassified = TRUE, returnstats = FALSE, ...)
#'
#' Plots relative abundance heatmaps annotated by the metadata
#' @export

plot_relabund_heatmap <- function(mgseqobj = NULL, glomby = NULL, heatpalette = "smart", hmtype = NULL, hmasPA = FALSE, compareby = NULL, invertbinaryorder = FALSE, ntop = NULL, ordercolsby = NULL, colcategories = NULL, cluster_rows = FALSE, subsetby = NULL, applyfilters = NULL, featmaxatleastPPM = 0, featcutoff = c(0, 0), samplesToKeep = NULL, featuresToKeep = NULL, adjustpval = FALSE, showonlypbelow = NULL, showpval = TRUE, showl2fc = TRUE, maxl2fc = NULL, minl2fc = NULL, genomecompleteness = NULL, list.data = NULL, addtit = NULL,  scaled = FALSE, mgSeqnorm = FALSE, cdict = NULL, maxnumheatmaps = NULL, numthreads = 4, nperm = 99, statsonlog = TRUE, ignoreunclassified = TRUE, returnstats = FALSE, ...){

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

    #Test for silly stuff
    if ((hmtype %in% c("comparative", "PA")) && (is.null(compareby))){
        stop("If rows are to be selected by highest significant difference between classes in a discrete category, you must determine the category using the argument *classesvector*")
    }

    if (!is.null(applyfilters)){
        if (applyfilters == "stringent"){
            if (analysis == "LKT"){
                featcutoff <- c(2000, 15)
                genomecompleteness <- 0.3
                minl2fc <- 2
            } else {
                featcutoff <- c(50, 15)
                genomecompleteness <- NULL
                minl2fc <- 2.5
            }
        } else if (applyfilters == "moderate"){
            if (analysis == "LKT"){
                featcutoff <- c(500, 10)
                genomecompleteness <- 0.1
                minl2fc <- 1.5
            } else {
                featcutoff <- c(10, 5)
                genomecompleteness <- NULL
                minl2fc <- 2
            }
        } else {
          featcutoff <- c(0, 0)
          genomecompleteness <- NULL
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

    #Initialize Stats Vector list
    svec <- NULL
    svec <- vector("list", length = 1000)
    s <- 1
    n <- 1

    numfeats <- nrow(MRcounts(obj))

    #subset by metadata column
    for (sp in 1:length(subset_points)){
        if (!(is.null(subsetby))){
            samplesToKeep <- which((pData(obj)[,which(colnames(pData(obj)) == subsetby)]) == subset_points[sp])
            print(paste("Plotting within", subset_points[sp]))
            colcategories <- colcategories[!(colcategories %in% subsetby)]
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
                cutoffmsg <- paste("Feature must be >", thresholdPPM, "PPM in at least ", sampcutoffpct, "% of samples", sep="")
            } else {
                cutoffmsg <- "Feature must be > 0 PPM in at least 0% of samples"
                featcutoff <- c(0,0)
            }

            if (all(c((!(is.null(featmaxatleastPPM))), (featmaxatleastPPM > 0)))) {
                minPPMmsg <- paste("Highest feature must be >", featmaxatleastPPM, "PPM", sep = " ")
            } else {
                minPPMmsg <- "Highest feature must be > 0 PPM"
                featmaxatleastPPM <- 0
            }

            if (hmtype == "PA"){
                hmtypemsg <- "Presence/Absence Heatmap"
                asPA <- hmasPA
                stattype <- "PA"
            } else {
                hmtypemsg <- "Relative Abundance Heatmap"
                asPA <- FALSE
                if (hmtype != "comparative"){
                    stattype <- "variance"
                    if(is.null(ntop)){
                        ntop <- 50
                    }
                } else {
                    stattype <- "auto"
                }
            }

            currobj <- filter_experiment(mgseqobj = obj, featmaxatleastPPM = featmaxatleastPPM, featcutoff = featcutoff, samplesToKeep = samplesToKeep, asPA = asPA, asPPM = TRUE, mgSeqnorm = mgSeqnorm)

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
            countmat <- MRcounts(currobj, norm = FALSE, log = !(asPA))

            #Protect against rows with empty data
            rowsToKeep <- which(rowSums(countmat) > 0 & rownames(countmat) != "")
            countmat <- countmat[rowsToKeep, ]

            if (ignoreunclassified == TRUE){
               dunno <- c(paste(analysis, "none", sep = "_"), "LKT__d__Unclassified", "LKT__Unclassified")
               rowsToKeep <- which(!(rownames(countmat) %in% dunno))
               countmat <- countmat[rowsToKeep, ]
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

            topcats <- nrow(countmat)
            if (!(is.null(ntop))) {
                topcats <- min(topcats, ntop)
            }

            #Calculate matrix stats and get new matrix.
            if (stattype == "variance"){
                matstats <- calculate_matrix_stats(countmatrix = countmat, stattype = "variance")
                #Bank raw stats
                matstatsbank <- as.data.frame(matstats)
                if (!(is.null(genomecompleteness))){
                    matstatsbank$Taxa <- rownames(matstatsbank)
                    genomecompletenessdf$MedianGenomeComp <- rowMedians(as.matrix(genomecompletenessdf))
                    genomecompletenessdf$SDGenomeComp <- rowSds(as.matrix(genomecompletenessdf))
                    genomecompletenessdf$Taxa <- rownames(genomecompletenessdf)
                    genomecompletenessdfmedians <- genomecompletenessdf[, c("Taxa", "MedianGenomeComp", "SDGenomeComp")]
                    matstatsbank <- left_join(matstatsbank, genomecompletenessdfmedians)
                    rownames(matstatsbank) <- matstatsbank$Taxa
                    matstatsbank$Taxa <- NULL
                }
                svec[[s]] <- matstatsbank

                matstats$Colour <- rep("black", nrow(matstats))
                countmat2 <- countmat[rownames(matstats), ]
                #Create a list of matrices each of maximum 50 rows
                rowlist <- split(1:topcats, ceiling(seq_along(1:topcats) / 50))
                matlist <- lapply(1:length(rowlist), function(x){countmat2[rowlist[[x]], ]})
                statslist <- lapply(1:length(rowlist), function(x){ matstats[rowlist[[x]], ] })
                rowlblcol_list <- lapply(1:length(rowlist), function(x){ rep("black", length(rowlist[[x]])) })
                stattit <- paste("Top", topcats, "most variant features across samples")
                statmsg <- stattype

            } else {
                cl = pData(currobj)[ , which(colnames(pData(currobj)) == compareby)]
                discretenames <- sort(unique(cl))

                matstats <- calculate_matrix_stats(countmatrix = countmat, uselog = TRUE, statsonlog = statsonlog, stattype = stattype, classesvector = cl, invertbinaryorder = invertbinaryorder, numthreads = numthreads, nperm = nperm)

                #Bank raw stats
                matstatsbank <- as.data.frame(matstats)
                if (!(is.null(genomecompleteness))){
                    matstatsbank$Taxa <- rownames(matstatsbank)
                    genomecompletenessdf$MedianGenomeComp <- rowMedians(as.matrix(genomecompletenessdf))
                    genomecompletenessdf$SDGenomeComp <- rowSds(as.matrix(genomecompletenessdf))
                    genomecompletenessdf$Taxa <- rownames(genomecompletenessdf)
                    genomecompletenessdfmedians <- genomecompletenessdf[, c("Taxa", "MedianGenomeComp", "SDGenomeComp")]
                    matstatsbank <- left_join(matstatsbank, genomecompletenessdfmedians)
                    rownames(matstatsbank) <- matstatsbank$Taxa
                    matstatsbank$Taxa <- NULL
                }
                svec[[s]] <- matstatsbank

                #Account for the fact that if there is only a single class within compareby variable, stats may have been coerced to variance.
                if (matstats$Method[1] == "variance") {
                    matstats$Colour <- rep("black", nrow(matstats))
                    countmat2 <- countmat[rownames(matstats), ]
                    #Create a list of matrices each of maximum 50 rows
                    topcats <- min(nrow(countmat2), 50)
                    rowlist <- split(1:topcats, ceiling(seq_along(1:topcats) / 50))
                    matlist <- lapply(1:length(rowlist), function(x){ countmat2[rowlist[[x]], ] })
                    statslist <- lapply(1:length(rowlist), function(x){ matstats[rowlist[[x]], ] })
                    rowlblcol_list <- lapply(1:length(rowlist), function(x){rep("black", length(rowlist[[x]]))})
                    stattit <- paste("Top", topcats, "most variant features across samples")
                    statmsg <- paste("Var", compareby, sep = "_")

                } else {

                    if (!(is.null(minl2fc)) & ("l2fc" %in% colnames(matstats))){
                        #Check if there is enough leftover after filter.
                        if (length(which(matstats$absl2fc > minl2fc)) < 2){
                            print(paste("There are less than 2 features which have >", minl2fc, "l2fc."))
                            minl2fc <- round((max((matstats$absl2fc[order(matstats$absl2fc, decreasing = TRUE)][min(length(matstats$absl2fc), 40)]), 0.5)), 1)
                        }
                        minl2fcmsg <- paste("log2foldchange >", minl2fc)
                        matstats <- subset(matstats, absl2fc > minl2fc)
                    } else {
                        minl2fcmsg <- "log2foldchange > 0"
                    }

                    if (!(is.null(maxl2fc)) & ("l2fc" %in% colnames(matstats))){
                        matstats <- subset(matstats, absl2fc < maxl2fc)
                        maxl2fcmsg <- paste("log2foldchange <", maxl2fc)
                    } else {
                        maxl2fcmsg <- "log2foldchange < Inf"
                    }

                    if ("l2fc" %in% colnames(matstats)){
                        matstats$Colour <- ifelse(matstats$l2fc < 0, "#900000", "#000000")
                        statmsg <- paste("MWW", compareby, sep="_")
                    } else if ("oddsRatio" %in% colnames(matstats)){
                        matstats$Colour <- ifelse(matstats$oddsRatio < 1, "#900000", "#000000")
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

                    countmat2 <- as.matrix(countmat[rownames(matstats), ])

                    #If adjustpval is set to auto, then find out which is best and re-set it to either TRUE or FALSE
                    if (("pval" %in% colnames(matstats)) && adjustpval == "auto"){
                        propsigadj <- length(which(matstats$padj_fdr < 0.05)) / length(matstats$padj_fdr)
                        propsignonadj <- length(which(matstats$pval < 0.05)) / length(matstats$pval)
                        fracsigadj <- propsigadj / propsignonadj
                        if ((!is.na(fracsigadj)) && (fracsigadj > 0.2)){
                            adjustpval = TRUE
                        } else {
                            adjustpval = FALSE
                        }
                    }

                    if (("pval" %in% colnames(matstats)) && adjustpval != TRUE){
                        sigmeas <- "pval"
                    } else if (("pval" %in% colnames(matstats)) && adjustpval == TRUE){
                        sigmeas <- "padj_fdr"
                    }

                    if (is.null(showonlypbelow)){
                        if (nrow(countmat2) < topcats){
                            topcats <- nrow(countmat2)
                            print("There are not enough features matching the criteria imposed.")
                            print(paste("Showing top", topcats, "features."))
                        }
                        countmat2 <- countmat2[1:topcats, ]
                        #Create a list of matrices each of maximum ~50 rows
                        if((topcats %% 50) == 1){
                            chunksize <- 48
                        } else {
                            chunksize <- 50
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
                            if((length(rowcutoff) %% 50) == 1){
                                chunksize <- 48
                            } else {
                                chunksize <- 50
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
                            matlist<-NULL
                            if (matstats$Method[1] == "MannWhitneyWilcoxon") {
                                stattit <- c(paste(sigmeas, "<", showonlypbelow, "different between"), paste(compareby, "using MannWhitneyWilcoxon"))
                            } else if (matstats$Method[1] == "permanova"){
                                stattit <- c(paste(sigmeas, "<", showonlypbelow, "different between"), paste(compareby, "using PERMANOVA"))
                            } else if (matstats$Method[1] == "anova"){
                                stattit <- c(paste(sigmeas, "<", showonlypbelow, "different between"), paste(compareby, "using ANOVA"))
                            } else if (matstats$Method[1] == "fisher"){
                                stattit <- c(paste(sigmeas, "<", showonlypbelow, "Present/Absent between"), paste(compareby, "using Fishers test"))
                            } #End conditional to define statistics title
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
                    fontsizey <- min(6, round((((-1 / 150) * (nrow(mathm))) + 1) * 6, 0))
                    fontsizex <- as.numeric(unlist(round((((-1 / 150) * (ncol(mathm))) + 1) * 5, 0)))

                    hmdf <- as.data.frame(matrix(data = 0, nrow = nrow(pData(currobj)), ncol = length(colcategories)))
                    cores <- vector("list", length = length(colcategories))
                    for (g in 1:length(colcategories)){
                        hmdf[ , g] <- pData(currobj)[ , which(colnames(pData(currobj)) == colcategories[g])]
                        colnames(hmdf)[g] <- colcategories[g]
                        if (!(is.numeric(hmdf[ ,g]))){
                            if (is.null(cdict)){
                                cores[[g]] <- as.vector(rainbow(length(unique(hmdf[ ,g]))))
                                names(cores[[g]]) <- sort(unique(hmdf[ ,g]))
                            } else {
                                ct <- cdict[[colcategories[g]]]
                                cores[[g]] <- as.vector(ct$Hex)
                                names(cores[[g]]) <- as.vector(ct$Name)
                            }
                        } else {
                            #Variable is numeric, but check for variance in the numbers
                            if ((max(hmdf[, g]) - min(hmdf[, g])) > 0 ){
                                cores[[g]] <- colorRamp2(c(0, max(hmdf[, g])), c("white", "midnightblue"))
                            } else {
                                cores[[g]] <- as.vector(rainbow(length(unique(hmdf[, g]))))
                                names(cores[[g]]) <- sort(unique(hmdf[, g]))
                            }
                        }
                        names(cores)[g] <- colcategories[g]
                    }

                    if (heatpalette == "sequential"){
                        heatmapCols <- colorRampPalette((brewer.pal(9, "YlOrRd")))(50)
                    } else if (heatpalette == "diverging"){
                        heatmapCols <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(50)
                    } else {
                        PctHmColours <- c("blue4", "blue", "slategray1", "khaki", "orange", "tomato", "red", "magenta2", "magenta4")
                        if (analysis == "LKT"){
                            RelabundBreakPts <- c(0.0001, 0.001, 0.1, 1, 2.5, 5, 10, 50, 100)
                            relabundscalename <- "Relative Abundance (%)"
                            RelabundBreakPtsLbls <- as.character(paste0(RelabundBreakPts, "%"))
                            HMrelabundBreaks <- Pct2log2PPM(RelabundBreakPts)
                        } else {
                            countmatdistrib <- sapply(1:nrow(mathm), function(x) { quantile(mathm[x, ], probs=c(0.20, 0.5, 0.95)) })
                            median20 <- median(countmatdistrib["20%", ])
                            median50 <- median(countmatdistrib["50%", ])
                            median95 <- median(countmatdistrib["95%", ])
                            HMrelabundBreaks <- c(0, (median20 / 2), median20, median50, ((median50 + median95) * (1/3)), ((median50 + median95) * (2/3)), median95,  (median95 * 2), max(mathm))

                            RelabundBreakPts <- round(((2 ^ (HMrelabundBreaks)) - 1), 0)
                            relabundscalename <- "Parts Per Million"
                            RelabundBreakPtsLbls <- as.character(paste(RelabundBreakPts, "PPM"))
                        }
                        heatmapCols <- colorRamp2(HMrelabundBreaks, PctHmColours)
                    }
                    ha_column <- HeatmapAnnotation(df = hmdf, col = cores, annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 7, col = "black"))

                    #Build plot title
                    plotit <- paste(maintit, stattit, cutoffmsg, minPPMmsg, sep = "\n")

                    if (all(c((any(!(c(is.null(minl2fc), is.null(maxl2fc))))), ("l2fc" %in% colnames(stathm))))){
                        l2fcmsg <- paste(minl2fcmsg, maxl2fcmsg, sep = " | ")
                        plotit <- paste(plotit, l2fcmsg, sep = "\n")
                    }

                    if (!(is.null(genomecompleteness))){
                        plotit <- paste(plotit, completenessmsg, sep = "\n")
                    }

                    #Add plot number if there is more than one heatmap matrix.
                    if (length(matlist) > 1){
                        hmcounter <- paste(hm, length(matlist), sep = "/")
                        plotit <- paste(plotit, paste("Heatmap", hmcounter), sep = "\n")
                    }

                    #Switch plot title if doing a dual heatmap with genome completeness
                    #Get genome completeness hmdf if in taxonomic space
                    ht1fs <- 10
                    hm1tit <- plotit
                    if (analysis == "LKT"){
                        gchmdf <- genomecompletenessdf[rownames(mathm), colnames(mathm)]
                        gchmdf <- as.matrix(gchmdf)
                        gchmdf <- gchmdf * 100
                        #Cull genome completeness to 400% max, to preserve scale across heatmaps
                        gchmdf[which(gchmdf[] > 400)] <- 400
                        if (ncol(mathm) < 51){
                            ht1fs <- 7
                            dualHMtit <- plotit
                            hm1tit <- "Relative Abundance"
                        }
                    }

                    if (scaled == TRUE) {
                        countmat2 <- t(apply(countmat2, 1, scale))
                        lname <- "scaling"
                    } else {
                        if(asPA == TRUE){
                            lname = "P/A"
                        } else {
                            lname = "log2(PPM)"
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
                    if ("oddsRatio" %in% colnames(stathm)){
                        ORamplitude <- paste0(round(stathm$oddsRatio, 2), paste0("(", paste(round(stathm$lower, 2), round(stathm$upper, 2), sep = "-"), ")"))
                    }

                    if (showl2fc == TRUE){
                        showl2fc = "text"
                    }

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

                    #Determine column order explicitly if required and draw heatmap
                    if (!(is.null(ordercolsby))){
                        co <- NULL
                        colgroups <- sort(unique(pData(currobj)[, which(colnames(pData(currobj)) == ordercolsby)]))
                        for(w in 1:length(unique(colgroups))){
                            wantedsamp <- pData(currobj)[]$Sample[which(pData(currobj)[, ordercolsby] == colgroups[w])]
                            #get temporary matrix with only one group to calculate distance
                            gmat <- mathm[, wantedsamp]
                            hcmat <- hclust(dist(t(gmat)))
                            gord <- hcmat$order
                            co <- append(co, colnames(gmat)[gord], after = length(co))
                        }
                        ht1 <- Heatmap(mathm, name = lname, cluster_columns = FALSE, column_order = co, column_title = hm1tit, column_title_gp = gpar(fontsize = ht1fs), top_annotation = ha_column, col = heatmapCols, column_names_gp = gpar(fontsize = fontsizex), right_annotation = row_ha, cluster_rows = cluster_rows, show_row_dend = FALSE, row_names_side="left", row_names_gp = gpar(fontsize = fontsizey, col = rowlblcol), heatmap_legend_param = list(direction = "horizontal", title = relabundscalename, labels = RelabundBreakPtsLbls, at = HMrelabundBreaks, title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 6)), row_names_max_width = unit(6, "cm"))
                    } else {
                        ht1 <- Heatmap(mathm, name = "Relative_Abundance", column_title = hm1tit, column_title_gp = gpar(fontsize = ht1fs), top_annotation = ha_column, col = heatmapCols, column_names_gp = gpar(fontsize = fontsizex), column_dend_height = unit(5, "mm"), right_annotation = row_ha, cluster_rows = cluster_rows, show_row_dend = FALSE, row_names_side = "left", row_names_gp = gpar(fontsize = fontsizey, col = rowlblcol), heatmap_legend_param = list(direction = "horizontal", title = relabundscalename, labels = RelabundBreakPtsLbls, at = HMrelabundBreaks, title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 6)), row_names_max_width = unit(6, "cm"))
                    }
                    #Make a genome completeness heatmap if in taxonomic space
                    if (analysis == "LKT"){
                        #Draw heatmap with completeness if taxonomic analysis
                        GCheatmapCols <- colorRamp2(c(0, 100, 200, 300, 400), c("white", "forestgreen", "blue", "firebrick1", "black"))
                        GCha_column <- HeatmapAnnotation(df = hmdf, col = cores, show_annotation_name = FALSE)
                        ht2 <- Heatmap(gchmdf, name = "GenComp", column_title = "% Genome completeness", column_title_gp = gpar(fontsize = ht1fs), top_annotation = GCha_column, col = GCheatmapCols, column_names_gp = gpar(fontsize = fontsizex), right_annotation = NULL, cluster_rows = FALSE, column_order = column_order(ht1), show_row_dend = FALSE, row_names_side = "left", row_names_gp = gpar(fontsize = fontsizey, col = rowlblcol), heatmap_legend_param = list(direction = "horizontal", title = "% GenComp", labels = c("0%", "100%", "200%", "300%", "> 400%"), title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 6)), row_names_max_width = unit(6, "cm"))

                        #Plot heatmaps side by side if there are 50 samples or less. Else plot one on each page.
                        if (ncol(mathm) < 51){
                            ht_list = ht1 + ht2
                            draw(ht_list, heatmap_legend_side = "bottom", annotation_legend_side = "right", ht_gap = unit(0.2, "cm"), padding = unit(c(2, 2, 2, 2), "mm"), column_title = dualHMtit, column_title_gp = gpar(fontsize = 10))
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
                            grid.text(sigmeas, y = unit(1, "npc") + unit(5, "mm"), just = "bottom", gp = gpar(fontsize = 5, col = "black"))
                        })
                    }
                    if (all(c(("l2fc" %in% colnames(stathm)), (showl2fc %in% c("text", "plot", TRUE))))) {
                        decorate_annotation("Log2FC", {
                            grid.text("Log2FC", y = unit(1, "npc") + unit(5, "mm"), just = "bottom", gp = gpar(fontsize = 5, col = "black"))
                        })
                    }
                    if ("oddsRatio" %in% colnames(stathm)) {
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
                plotit <- c(stattit, paste(cutoffmsg, minPPMmsg, sep=", "))
                plot.new()
                grid.table(c(maintit, "No features fulfilling", plotit), rows = NULL, cols = NULL, theme = ttheme_default(base_size = 10))
                #gvec[[n]]<-recordPlot()
                n <- n + 1
            } #End conditional if there are heatmaps to plot from stats matrix

        } else {

            if(length(samplesToKeep) < 2){
                #There is only one sample in subset, so nothing was done.
                print(paste(as.character(colnames(pData(obj))[samplesToKeep]), "was the only sample within", subsetname, "subset. Must have at least two samples for a heatmap."))
            }
            if(numfeats < 2){
                #There is only one sample in subset, so nothing was done.
                print(paste(as.character(rownames(MRcounts(obj))), "was the only feature within", subsetname, "subset. Must have at least two features for a heatmap."))
            }

        } #End conditional if there is more than a single sample and or single feature in subset

    } #End subset by loop

    #Redefine stats list as ones only containing data
    svec <- svec[sapply(svec, function(x){ !(is.null(x)) } )]

    if(returnstats == TRUE){
        return(svec)
    } else {
        return(print("Heatmaps generated."))
    }
}
