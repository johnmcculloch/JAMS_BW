#' plot_feature_relabund_devel()
#'
#' Plots correlation heatmaps annotated by the metadata or a correlelogram of features
#' @export

plot_feature_relabund_devel <- function(mgseqobj = NULL, glomby = NULL, mgSeqnorm = FALSE, features = NULL, compareby = NULL, colourby = NULL, shapeby = NULL, subsetby = NULL, uselog = TRUE, statmeth = NULL, samplesToKeep = NULL, signiflabel = "p.signif", plottitle = NULL, asPA = FALSE, subsetbytaxlevel = NULL, taxtable = NULL, list.data = NULL, cdict = NULL, max_categories = 3, ...) {

    #Get appropriate object to work with
    obj <- mgseqobj

    #Exclude samples and features if specified
    if (!(is.null(samplesToKeep))){
        obj <- obj[ , samplesToKeep]
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
                minabscorrcoeff <- 0.8
            } else {
                featcutoff <- c(50, 15)
                genomecompleteness <- NULL
                minabscorrcoeff <- 0.8
            }
        } else if (applyfilters == "moderate"){
            if (analysis == "LKT"){
                featcutoff <- c(500, 10)
                genomecompleteness <- 0.1
                minabscorrcoeff <- 0.5
            } else {
                featcutoff <- c(10, 5)
                genomecompleteness <- NULL
                minabscorrcoeff <- 0.5
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

    if ((analysis == "LKT") && (!(is.null(subsetbytaxlevel)))){
        stop("Plot taxa is only for functional (not taxonomic) analyses. It will additionally plot which taxa bear which of the funcitons of interest. If you are trying to plot only a certain taxon itself, then use the feature argument with a taxonomical metagenomeSeq experiment.")
    }

    if (!(is.null(subsetby))){
        subset_points <- sort(unique((pData(obj)[, which(colnames(pData(obj)) == subsetby)])))
    } else {
        subset_points <- "none"
    }

    #Initialize Stats and Graph Vector lists
    svec <- NULL
    svec <- vector("list", length = 1000)
    gvec <- NULL
    gvec <- vector("list", length = 1000)
    s <- 1
    n <- 1
    plotcount <- 1

    numfeats <- nrow(MRcounts(obj))

    #subset by metadata column
    for (sp in 1:length(subset_points)) {
        if (!(is.null(subsetby))){
            samplesToKeep <- which((pData(obj)[, which(colnames(pData(obj)) == subsetby)]) == subset_points[sp])
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
                cutoffmsg <- paste("Feature must be >", thresholdPPM, "PPM in at least ", sampcutoffpct, "% of samples", sep = "")
            } else {
                cutoffmsg <- "Feature must be > 0 PPM in at least 0% of samples"
                featcutoff <- c(0, 0)
            }

            if (all(c((!(is.null(featmaxatleastPPM))), (featmaxatleastPPM > 0)))) {
                minPPMmsg <- paste("Highest feature must be >", featmaxatleastPPM, "PPM", sep = " ")
            } else {
                minPPMmsg <- "Highest feature must be > 0 PPM"
                featmaxatleastPPM <- 0
            }

            currobj <- filter_experiment(mgseqobj = obj, featmaxatleastPPM = featmaxatleastPPM, featcutoff = featcutoff, samplesToKeep = samplesToKeep, asPA = FALSE, asPPM = TRUE, mgSeqnorm = mgSeqnorm)

            #Cull features to ones desired
            if (!is.null(features)) {
                currobj <- currobj[features, ]
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

            #Find out what kind of a thing we're comparing by.
            if (!is.null(compareby)){
                lookslike <- class(pData(currobj)[ , compareby])
                if (lookslike == "numeric"){
                    stattype <- "spearman"
                } else {
                    stattype <- "auto"
                }
            }

            #Coerce stattype to what the user wants, if not NULL. User has to know what (s)he is doing.
            if (!is.null(stattype)){

            }

                #Get correlation to a continuous variable in the metadata.
                #Check if variable is really numeric
                continuousvec <- pData(currobj)[ , which(colnames(pData(currobj)) == compareby)]
                if (!is.numeric(continuousvec)){
                    stop(paste("Variable", compareby, "of the metadata is NOT continuous. Use function plot_relabund_heatmap for discrete variables. Aborting now."))
                }
                print(paste("Calculating correlation coefficients of features to", compareby))
                matstats <- calculate_matrix_stats(countmatrix = countmat, uselog = FALSE, statsonlog = FALSE, stattype = stattype, classesvector = continuousvec)

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

                #Plot information obtained in matstats
                if (!is.null(minabscorrcoeff)){
                    matstats <- subset(matstats, abscorrel >= minabscorrcoeff)
                }

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

                if (!is.null(showonlypbelow)){
                    matstats <- matstats[(which(matstats[, sigmeas] <= showonlypbelow)), ]
                }

                #only plot if there is more than one thing to plot.
                if (nrow(matstats) > 0) {
                    maxnumfeats <- min(topcats, nrow(matstats))
                    #plot a scatterplot line by line
                    for (f in 1:maxnumfeats){
                        plotdf <- NULL
                        scatp <- NULL
                        featint <- rownames(matstats)[f]
                        currcts <- countmat[featint, ]
                        plotdf <- data.frame(Samples = names(currcts), Relabund = currcts, stringsAsFactors = FALSE)
                        plotdf$Continuous <- as.numeric(pData(currobj)[(match(plotdf$Samples, rownames(pData(currobj)))), compareby])
                        scatp <- ggplot(plotdf, aes(x = Continuous, y = Relabund)) + geom_point(pch = 19)
                        scatp <- scatp + geom_smooth(method = lm, color = "red", se = FALSE)
                        scatp <- scatp + theme_minimal()
                        correlmsg <- paste(corstat, round(matstats$correl[f], 2), sep = "=")
                        pvalmsg <- paste(sigmeas, round(matstats[, sigmeas][f], 3), sep = "=")
                        statmsg <- paste(correlmsg, pvalmsg)
                        maintit <- paste(hmtypemsg, (paste(featint, "by", compareby)), paste("Subset:", subset_points[sp]), statmsg, sep = "\n")
                        scatp <- scatp + ggtitle(maintit)
                        scatp <- scatp + labs(x = compareby, y = "Relative Abundance in PPM")
                        gvec[[plotcount]] <- scatp
                        plotcount <- plotcount + 1
                    }
                }

            } else {
                matstats <- calculate_matrix_stats(countmatrix = countmat, uselog = FALSE, statsonlog = FALSE, stattype = stattype, classesvector = NULL)
            }
        } #End conditional if there are any features left over after filtering
    } #End for loop for plotting within each subset point

    #Redefine stats list as ones only containing data
    svec <- svec[sapply(svec, function(x){ !(is.null(x)) } )]
    gvec <- gvec[sapply(gvec, function(x){ !(is.null(x)) } )]
    if (!is.null(gvec)){
        print(gvec)
    }

    if (returnstats == TRUE){
        return(svec)
    } else {
        return(print("Heatmaps generated."))
    }
}
