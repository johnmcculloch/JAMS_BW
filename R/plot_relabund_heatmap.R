#' plot_relabund_heatmap(mgseqobj=NULL, glomby=NULL, heatpalette="diverging", hmtype=NULL, hmasPA=FALSE, compareby=NULL, invertbinaryorder=FALSE, ntop=NULL, ordercolsby=NULL, colcategories=NULL, cluster_rows=FALSE, subsetby=NULL, applyfilters=NULL, featmaxatleastPPM=0, featcutoff=c(0, 0), samplesToKeep=NULL, featuresToKeep=NULL, adjustpval=FALSE, showonlypbelow=NULL, showpval=TRUE, showl2fc=TRUE, maxl2fc=NULL, minl2fc=NULL, genomecompleteness=NULL, list.data=NULL, addtit=NULL, scaled=FALSE, mgSeqnorm=FALSE, cdict=NULL, maxnumheatmaps=NULL, numthreads=4, nperm=99, statsonlog=TRUE, ignoreunclassified=TRUE, returnstats=TRUE, ...)
#'
#' Plots relative abundance heatmaps annotated by the metadata
#' @export

plot_relabund_heatmap <- function(mgseqobj=NULL, glomby=NULL, heatpalette="diverging", hmtype=NULL, hmasPA=FALSE, compareby=NULL, invertbinaryorder=FALSE, ntop=NULL, ordercolsby=NULL, colcategories=NULL, cluster_rows=FALSE, subsetby=NULL, applyfilters=NULL, featmaxatleastPPM=0, featcutoff=c(0, 0), samplesToKeep=NULL, featuresToKeep=NULL, adjustpval=FALSE, showonlypbelow=NULL, showpval=TRUE, showl2fc=TRUE, maxl2fc=NULL, minl2fc=NULL, genomecompleteness=NULL, list.data=NULL, addtit=NULL, scaled=FALSE, mgSeqnorm=FALSE, cdict=NULL, maxnumheatmaps=NULL, numthreads=4, nperm=99, statsonlog=TRUE, ignoreunclassified=TRUE, returnstats=TRUE, ...){

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
        obj <- agglomerate_features(mgseqobj=obj, glomby=glomby)
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
        subset_points <- sort(unique((pData(obj)[, which(colnames(pData(obj))==subsetby)])))
    } else {
        subset_points <- "none"
    }

    #Initialize Stats Vector list
    svec <- NULL
    svec <- vector("list", length = 1000)
    s <- 1
    n <- 1

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

        numfeats <- nrow(MRcounts(obj))

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

            if(!(is.null(featmaxatleastPPM))){
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
                maintit <- paste(hmtypemsg, analysisname, paste("within", subset_points[sp]), sep=" | ")
            } else {
                maintit <- paste(hmtypemsg, analysisname, sep=" | ")
            }
            if (!is.null(addtit)) {
                maintit <- paste(addtit, maintit, sep="\n")
            }

            #Get counts matrix
            countmat <- MRcounts(currobj, norm = FALSE, log = !(asPA))

            #Protect against rows with empty data
            rowsToKeep <- which(rowSums(countmat) > 0 & rownames(countmat) != "")
            countmat <- countmat[rowsToKeep, ]

            if (ignoreunclassified == TRUE){
               dunno <- c(paste(analysis, "none", sep="_"), "LKT__d__Unclassified")
               rowsToKeep <- which(!(rownames(countmat) %in% dunno))
               countmat <- countmat[rowsToKeep, ]
            }

            #Rename rows to include description if not taxonomic data
            if (analysis != "LKT"){
                feattable <- fData(currobj)
                feattable$Feature <- paste(feattable$Accession, feattable$Description, sep = "-")
                rownames(countmat) <- feattable$Feature[match(rownames(countmat), feattable$Accession)]
                #Allow up to 40 characters for readability. Full names can be seen in the spreadsheet.
                rownames(countmat) <- strtrim(rownames(countmat2), 40)
            }
            matrixSamples <- colnames(countmat)
            matrixRows <- rownames(countmat)

            #Discard taxa below required level of completeness
            if (!(is.null(genomecompleteness))){
                print(paste("Genome completeness must be", genomecompleteness, "in at least one sample"))
                genomecompletenessdf <- get_genome_completeness(pheno = pData(currobj), list.data = list.data)
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
                    genomecompletenessdf <- genomecompletenessdf[, c("Taxa", "MedianGenomeComp", "SDGenomeComp")]
                    matstatsbank <- left_join(matstatsbank, genomecompletenessdf)
                    rownames(matstatsbank) <- matstatsbank$Taxa
                    matstatsbank$Taxa <- NULL
                }
                svec[[s]] <- matstatsbank

                matstats$Colour <- rep("black", nrow(matstats))
                countmat2 <- countmat[rownames(matstats), ]
                #Add rank variance to matrix names
                for (r in 1:nrow(countmat2)){
                    feature <- rownames(countmat2)[r]
                    rank <- getOrdinalNumber1(as.numeric(which(rownames(countmat2) == feature)))
                    rownames(countmat2)[r] <- paste(feature, rank, sep = "-")
                }
                #Create a list of matrices each of maximum 50 rows
                rowlist <- split(1:topcats, ceiling(seq_along(1:topcats)/50))
                matlist <- lapply(1:length(rowlist), function(x){countmat2[rowlist[[x]], ]})
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
                    genomecompletenessdf <- genomecompletenessdf[, c("Taxa", "MedianGenomeComp", "SDGenomeComp")]
                    matstatsbank <- left_join(matstatsbank, genomecompletenessdf)
                    rownames(matstatsbank) <- matstatsbank$Taxa
                    matstatsbank$Taxa <- NULL
                }
                svec[[s]] <- matstatsbank

                #Account for the fact that if there is only a single class within compareby variable, stats may have been coerced to variance.
                if (matstats$Method[1] == "variance") {
                    matstats$Colour <- rep("black", nrow(matstats))
                    countmat2 = countmat[rownames(matstats), ]
                    #Add rank variance to matrix names
                    for (r in 1:nrow(countmat2)){
                        feature <- rownames(countmat2)[r]
                        rank <- getOrdinalNumber1(as.numeric(which(rownames(countmat2) == feature)))
                        rownames(countmat2)[r] <- paste(feature, rank, sep = "-")
                    }
                    #Create a list of matrices each of maximum 50 rows
                    topcats <- min(nrow(countmat2), 50)
                    rowlist <- split(1:topcats, ceiling(seq_along(1:topcats) / 50))
                    matlist <- lapply(1:length(rowlist), function(x){ countmat2[rowlist[[x]], ] })
                    rowlblcol_list <- lapply(1:length(rowlist), function(x){rep("black", length(rowlist[[x]]))})
                    stattit <- paste("Top", topcats, "most variant features across samples")
                    statmsg <- paste("Var", compareby, sep="_")
                } else {

                    if(!(is.null(minl2fc)) & ("l2fc" %in% colnames(matstats))){
                        #Check if there is enough leftover after filter.
                        if (length(which(matstats$absl2fc > minl2fc)) < 2){
                            print(paste("There are less than 2 features which have >", minl2fc, "l2fc."))
                            minl2fc <- round((max((matstats$absl2fc[order(matstats$absl2fc, decreasing=TRUE)][min(length(matstats$absl2fc),40)]), 0.5)), 1)
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

                    if ((showl2fc==TRUE) & ("l2fc" %in% colnames(matstats))){
                        rownames(countmat2) <- paste(rownames(countmat2), " | L2FC=  ", round(matstats$l2fc, 1),sep = "")
                    }

                    #Add Odds Ratio if PA
                    if ("oddsRatio" %in% colnames(matstats)){
                        rownames(countmat2) <- paste(rownames(countmat2), " | OR=  ", paste(round(matstats$oddsRatio, 1), "(", round(matstats$lower, 1), "-", round(matstats$upper, 1), ")", sep= ""), sep = "")
                    }

                    #If adjustpval is set to auto, then find out which is best and re-set it to either TRUE or FALSE
                    if (("pval" %in% colnames(matstats)) && adjustpval == "auto"){
                        propsigadj <- length(which(matstats$padj_fdr < 0.05)) / length(matstats$padj_fdr)
                        propsignonadj <- length(which(matstats$pval < 0.05)) / length(matstats$pval)
                        fracsigadj <- propsigadj / propsignonadj
                        if((!is.na(fracsigadj)) && (fracsigadj > 0.2)){
                            adjustpval = TRUE
                        } else {
                            adjustpval = FALSE
                        }
                    }

                    if (("pval" %in% colnames(matstats)) && adjustpval != TRUE){
                        sigmeas <- "pval"
                        if (showpval == TRUE){
                            rownames(countmat2) <- paste((paste(paste0(sigmeas, "="), round(matstats[,sigmeas], 3), "|", sep = "")), rownames(countmat2))
                        }
                    } else if (("pval" %in% colnames(matstats)) && adjustpval == TRUE){
                        sigmeas <- "padj_fdr"
                        if (showpval == TRUE){
                            rownames(countmat2) <- paste((paste(paste0(sigmeas, "="), round(matstats[,sigmeas], 3), "|", sep = "")), rownames(countmat2))
                        }
                    }

                    if (is.null(showonlypbelow)){
                        if (nrow(countmat2) < topcats){
                            topcats <- nrow(countmat2)
                            print("There are not enough features matching the criteria imposed.")
                            print(paste("Showing top", topcats, "features."))
                        }
                        countmat2 <- countmat2[1:topcats, ]
                        #Create a list of matrices each of maximum 50 rows
                        rowlist <- split(1:topcats, ceiling(seq_along(1:topcats)/50))
                        matlist <- lapply(1:length(rowlist), function(x){ countmat2[rowlist[[x]], ] })
                        rowlblcol_list <- lapply(1:length(rowlist), function(x){matstats$Colour[rowlist[[x]]]})

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
                            rowlist <- split(rowcutoff, ceiling(seq_along(rowcutoff)/50))
                            matlist <- lapply(1:length(rowlist), function(x){ countmat2[rowlist[[x]], ] })
                            rowlblcol_list <- lapply(1:length(rowlist), function(x){matstats$Colour[rowlist[[x]]]})
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
                    rowlblcol <- rowlblcol_list[[hm]]
                    #Plot the heatmap
                    fontsizey <- max(4, round((((-1/150)*(nrow(mathm)))+1)*5, 0))
                    fontsizex <- max(3, as.numeric(unlist(round((((-1/150)*(ncol(mathm)))+1)*5, 0))))

                    hmdf <- as.data.frame(matrix(data=0, nrow=nrow(pData(currobj)), ncol=length(colcategories)))
                    cores <- vector("list",length=length(colcategories))
                    for (g in 1:length(colcategories)){
                        hmdf[ , g]<-pData(currobj)[ , which(colnames(pData(currobj))==colcategories[g])]
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

                    if(is.null(heatpalette) || heatpalette=="sequential"){
                        heatmapCols = colorRampPalette((brewer.pal(9, "YlOrRd")))(50)
                    } else {
                        heatmapCols = colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(50)
                    }
                    ha_column = HeatmapAnnotation(df = hmdf, col = cores)

                    #Build plot title
                    plotit<-paste(maintit, stattit, cutoffmsg, minPPMmsg, sep="\n")

                    if(stattype=="binary"){
                        l2fcmsg<-paste(minl2fcmsg, maxl2fcmsg, sep=" | ")
                        plotit<-paste(plotit, l2fcmsg, sep="\n")
                    }

                    if(!(is.null(genomecompleteness))){
                        plotit<-paste(plotit, completenessmsg, sep="\n")
                    }

                    #Add plot number if there is more than one heatmap matrix.
                    if(length(matlist)>1){
                        hmcounter<-paste(hm, length(matlist), sep="/")
                        plotit<-paste(plotit, paste("Heatmap", hmcounter), sep="\n")
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

                    #Determine column order explicitly if required and draw heatmap
                    if (!(is.null(ordercolsby))){
                        co<-NULL
                        colgroups<-sort(unique(pData(currobj)[,which(colnames(pData(currobj))==ordercolsby)]))
                        for(w in 1:length(unique(colgroups))){
                            wantedsamp<-pData(currobj)[]$Sample[which(pData(currobj)[,ordercolsby]==colgroups[w])]
                            #get temporary matrix with only one group to calculate distance
                            gmat<-mathm[ ,wantedsamp]
                            hcmat<-hclust(dist(t(gmat)))
                            gord<-hcmat$order
                            co<-append(co, colnames(gmat)[gord], after=length(co))
                        }
                        ht1 = Heatmap(mathm, name = lname, cluster_columns = FALSE, column_order = co, column_title = plotit, column_title_gp = gpar(fontsize = 10), top_annotation = ha_column, col = heatmapCols, column_names_gp = gpar(fontsize = fontsizex), cluster_rows = cluster_rows, row_names_side="left", row_names_gp = gpar(fontsize = fontsizey, col = rowlblcol), heatmap_legend_param = list(title = lname, title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 10)))
                    } else {
                        ht1 = Heatmap(mathm, name = lname, column_title = plotit, column_title_gp = gpar(fontsize = 10), top_annotation = ha_column, col = heatmapCols, column_names_gp = gpar(fontsize = fontsizex), cluster_rows = cluster_rows, row_names_side="left", row_names_gp = gpar(fontsize = fontsizey, col = rowlblcol), heatmap_legend_param = list(title = lname, title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 10)))
                    }

                    #Draw the heatmap
                    par(oma = c(10,7,3,10)+0.1, xpd = TRUE)
                    draw(ht1, heatmap_legend_side = "right", annotation_legend_side = "right", padding = unit(c(4, 2, 2, 2), "mm"))

                    #Print what the column annotations are
                    for(an in colnames(hmdf)) {
                        decorate_annotation(an, {
                            grid.text(an, unit(0, "npc") - unit(2, "mm"), 0.5, default.units = "npc", just = "right", gp=gpar(fontsize=10, col="black"))
                        })
                    }
                    #gvec[[n]]<-recordPlot()
                    n <- n + 1

                    if(n > 100){
                        stop("There are too many combinations. I think you have had enough plots. I am stopping here.")
                    }
                }
            } else {
                #There is no valid matrix, so print out the conditions which led the matrix to be empty.

                plotit<-c(stattit, paste(cutoffmsg, minPPMmsg, sep=", "))
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
    svec<-svec[sapply(svec, function(x){!(is.null(x))})]

    if(returnstats == TRUE){
        return(svec)
    } else {
        return(print("Heatmaps generated."))
    }
}
