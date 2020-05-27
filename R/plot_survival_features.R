#' plot_relabund_features(ExpObj = NULL, glomby = NULL, samplesToKeep = NULL, featuresToKeep = NULL, aggregatefeatures = FALSE, aggregatefeatures_label = "Sum_of_wanted_features", subsetby = NULL, bin_names = NULL, survivaltime = NULL, Samples_to_censor = NULL, conf.int = FALSE, facetby = NULL, wrap_facet = FALSE, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL, ntop = NULL, minabscorrcoeff = NULL, adjustpval = TRUE, padjmeth = "fdr", showonlypbelow = NULL, showonlypadjusted = FALSE, addtit = NULL, PPM_normalize_to_bases_sequenced = FALSE, uselog = FALSE, dump_interpro_descriptions_to_plot = FALSE, ignoreunclassified = TRUE, class_to_ignore = "N_A", return_plots = TRUE, ...)
#'
#' Generates relative abundance plots per feature annotated by the metadata using as input a SummarizedExperiment object
#' @export

plot_Kaplan_Meier_features <- function(ExpObj = NULL, glomby = NULL, samplesToKeep = NULL, featuresToKeep = NULL, aggregatefeatures = FALSE, aggregatefeatures_label = "Sum_of_wanted_features", subsetby = NULL, bin_names = NULL, survivaltime = NULL, Samples_to_censor = NULL, conf.int = FALSE, facetby = NULL, wrap_facet = FALSE, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL, ntop = NULL, minabscorrcoeff = NULL, adjustpval = TRUE, padjmeth = "fdr", showonlypbelow = NULL, showonlypadjusted = FALSE, addtit = NULL, PPM_normalize_to_bases_sequenced = FALSE, uselog = FALSE, dump_interpro_descriptions_to_plot = FALSE, ignoreunclassified = TRUE, class_to_ignore = "N_A", return_plots = TRUE, ...){

    flog.warn("This function is experimental in JAMS. Use at your own risk.")
    require(survival)
    require(survminer)

    variables_to_fix <- survivaltime

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

            hmtypemsg <- "Kaplan-Meier Plot"

            currobj <- filter_experiment(ExpObj = obj, featcutoff = presetlist$featcutoff, samplesToKeep = samplesToKeep, featuresToKeep = NULL, asPPM = TRUE, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, GenomeCompletenessCutoff = presetlist$GenomeCompletenessCutoff, PctFromCtgscutoff = presetlist$PctFromCtgscutoff)

        } else {

            flog.info("Unable to make plots with the current metadata for this comparison.")
            return(NULL)

        }

        curr_pt <- as.data.frame(colData(currobj))
        pheno_survival <- curr_pt[ , unique(c("Sample", survivaltime))]

        #Add censoring data
        pheno_survival$Censoring <- 1
        if (!is.null(Samples_to_censor)){
            pheno_survival$Censoring[which(pheno_survival$Sample %in% Samples_to_censor)] <- 0
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

                wantedfeatures <- aggregatefeatures_label
            }

            matrixSamples <- colnames(countmat)
            matrixRows <- rownames(countmat)

            #Transform continuous into categorical
            continuous2discrete <- function(wantedfeatures_PPM_per_Sample = NULL, feat = NULL, bin_names = NULL){
                wantedtaxon_PPM_per_Sample_distrib <- quantile(wantedfeatures_PPM_per_Sample[ , feat], probs = round(((1:length(bin_names)) * (1 / length(bin_names))), 2))
                if (is.redundant(c(0, as.numeric(wantedtaxon_PPM_per_Sample_distrib)))){
                    flog.info(paste("Feaure", feat, "doesn't have enough variance for being binned into", length(bin_names), "categories and will be discarded."))
                    TaxonPPMcats <- NULL
                    Breakinfo <- NULL
                } else {
                    TaxonPPMcats <- cut(wantedfeatures_PPM_per_Sample[ , feat], breaks = c(0, as.numeric(wantedtaxon_PPM_per_Sample_distrib)), include.lowest = TRUE, labels = bin_names)
                    TaxonPPMcats <- as.factor(TaxonPPMcats)
                    PPMbreaks <- c(0, round(as.numeric(wantedtaxon_PPM_per_Sample_distrib), 0))
                    Breakinfo <- NULL
                    for (bn in 1:length(bin_names)){
                        Breakinfo[bn] <- paste(bin_names[bn], paste(paste(PPMbreaks[bn], "PPM"), paste(PPMbreaks[(bn + 1)], "PPM"), sep = " to "), sep = " = ")
                    }
                }
                discretelist <- list()
                discretelist$TaxonPPMcats <- TaxonPPMcats
                discretelist$Breakinfo <- Breakinfo

                return(discretelist)
            }

            #Make a list of Break infos
            wantedfeatures_Breakinfos <- as.data.frame(t(countmat[wantedfeatures, ]))
            wantedfeatures_Breakinfos_list <- lapply(colnames(wantedfeatures_Breakinfos), function (x) { continuous2discrete(wantedfeatures_PPM_per_Sample = wantedfeatures_Breakinfos, feat = x, bin_names = bin_names)[]$Breakinfo } )
            names(wantedfeatures_Breakinfos_list) <- colnames(wantedfeatures_Breakinfos)

            #Bin PPM into discrete
            wantedfeatures_PPM_per_Sample <- as.data.frame(t(countmat[wantedfeatures, ]))
            for (colm in colnames(wantedfeatures_PPM_per_Sample)){
                wantedfeatures_PPM_per_Sample[ , colm] <- continuous2discrete(wantedfeatures_PPM_per_Sample = wantedfeatures_PPM_per_Sample, feat = colm, bin_names = bin_names)[]$TaxonPPMcats
            }

            wantedfeatures_PPM_per_Sample$Sample <- rownames(wantedfeatures_PPM_per_Sample)
            #Left join PPM
            pheno_survival <- left_join(pheno_survival, wantedfeatures_PPM_per_Sample, by = "Sample")
            rownames(pheno_survival) <- pheno_survival$Sample

            surv_object <- Surv(time = pheno_survival[ , survivaltime], event = pheno_survival$Censoring)

            #Not the most elegant way of doing it, but being safe.
            fitlist <- list()
            pvals <- NULL
            validfeats <- colnames(pheno_survival)[!(colnames(pheno_survival) %in% c("Sample", survivaltime, "Censoring"))]
            for (nfeat in 1:length(validfeats)){
                feat <- validfeats[nfeat]
                #Use the old trick of regenerating the dataframe and changing variable to common name
                curr_pheno_survival <- pheno_survival[ ,c(survivaltime, "Censoring", feat)]
                colnames(curr_pheno_survival)[which(colnames(curr_pheno_survival) == feat)] <- "FeatStratum"
                featfit <- NULL
                featfit <- survfit(surv_object ~ FeatStratum, data = curr_pheno_survival)
                fitlist[[nfeat]] <- featfit
                names(fitlist)[nfeat] <- feat
                pv <- surv_pvalue(featfit, data = curr_pheno_survival, method = "survdiff", test.for.trend = FALSE)[]$pval
                pvals[nfeat] <- pv
                names(pvals)[nfeat] <- feat
            }
            matstats <- data.frame(Accession = names(pvals), pval = as.numeric(pvals), stringsAsFactors = FALSE)
            rownames(matstats) <- matstats$Accession

            #Adjust p-values
            adjlist <- lapply(p.adjust.methods, function(x){ p.adjust(matstats$pval, method = x, n = length(matstats$pval)) })
            names(adjlist) <- paste("padj", p.adjust.methods, sep = "_")
            adjdf <- as.data.frame(adjlist)
            matstats <- cbind(matstats, adjdf)
            matstats <- matstats[order(matstats$pval, decreasing = FALSE), ]

            matstats$Method <- rep("survdiff", nrow(matstats))

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

        } else {

            #abort, nothing is left over
            flog.warn("None of the wanted features were found in SummarizedExperiment object when using the current filtration parameters.")

            return(NULL)
        }

        flog.info("Plotting results...")
        feat <- NULL
        for (feat in rownames(matstats)){
            flog.info(paste("Plotting", feat))
            curr_pheno_survival <- pheno_survival[ ,c(survivaltime, "Censoring", feat)]
            colnames(curr_pheno_survival)[which(colnames(curr_pheno_survival) == feat)] <- "FeatStratum"

            fit1 <- fitlist[[feat]]
            breakinfo <- wantedfeatures_Breakinfos_list[[feat]]
            #Build plot title
            overallpmeth <- matstats[feat, "Method"]
            overallp <- paste0("pval=", signif(matstats[feat, "pval"], digits = 3))
            overalladjp <- paste0("padj_fdr=", signif(matstats[feat, "padj_fdr"], digits = 3))
            stattit <- paste(overallpmeth, overallp, overalladjp, ffeatmsg, sep = " | ")

            p <- ggsurvplot(fit1, data = curr_pheno_survival, pval = overallp, conf.int = conf.int, ggtheme = theme_minimal(), legend = "right", legend.title = feat, legend.labs = breakinfo, risk.table = FALSE)

            #Add description to feature, if applicable
            if (analysis != "LKT"){
                featdesc <- rowData(currobj)[feat, "Description"]
                featname <- paste(feat, featdesc)
            } else {
                featname <- feat
            }

            msgs <- c(maintit, featname, stattit)
            plotit <- paste0(msgs, collapse = "\n")

            if (!is.null(addtit)){
                plotit <- paste(plotit, addtit, sep = "\n")
            }

            p <- p + labs(title = plotit)

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
        }

    }#End loop for each subset

    if (return_plots){
        #Return plots, as nothing was printed
        return(gvec)
    }
}
