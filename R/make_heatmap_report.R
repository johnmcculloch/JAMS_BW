#' make_heatmap_report(report = "comparative", project = NULL, expvec = NULL, usefulexp = NULL, appendtofilename = NULL, glomby = NULL, featcutofftaxa = NULL, featcutofffunct = NULL, l2fchmtaxa = NULL, l2fchmfunc = NULL, samplesToKeep = NULL, featuresToKeep = NULL, applyfilters = NULL, absolutemaxfeats = NULL, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL, cluster_samples_per_heatmap = FALSE, cluster_features_per_heatmap = FALSE, PPM_normalize_to_bases_sequenced = FALSE, no_underscores = FALSE, hmasPA = FALSE, threshPA = 0, splitcolsby = NULL, splitcolsbycomp = FALSE, ordercolsby = NULL, textby = NULL, cluster_rows = TRUE, max_rows_in_heatmap = 50, secondaryheatmap = "GenomeCompleteness", minabscorrcoeff = NULL, variable_list = NULL, max_annot_legends = 4, includereaddata = TRUE, adjustpval = FALSE, scaled = FALSE, cdict = NULL, showonlypbelow = 0.05, makespreadsheets = FALSE, makeheatmaps = TRUE, maxnumheatmaps = NULL, numthreads = 4, nperm = 99, class_to_ignore = NULL, ...)
#'
#' Generates standard comparative heatmaps and spreadsheet based on relative abundance difference between categories for analyses named in expvectors.
#' @export

make_heatmap_report <- function(report = "comparative", project = NULL, expvec = NULL, usefulexp = NULL, appendtofilename = NULL, glomby = NULL, featcutofftaxa = NULL, featcutofffunct = NULL, l2fchmtaxa = NULL, l2fchmfunc = NULL, samplesToKeep = NULL, featuresToKeep = NULL, applyfilters = NULL, absolutemaxfeats = NULL, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL, cluster_samples_per_heatmap = FALSE, cluster_features_per_heatmap = FALSE, PPM_normalize_to_bases_sequenced = FALSE, no_underscores = FALSE, hmasPA = FALSE, threshPA = 0, splitcolsby = NULL, splitcolsbycomp = FALSE, ordercolsby = NULL, textby = NULL, cluster_rows = TRUE, max_rows_in_heatmap = 50, secondaryheatmap = "GenomeCompleteness", minabscorrcoeff = NULL, variable_list = NULL, max_annot_legends = 4, includereaddata = TRUE, adjustpval = FALSE, scaled = FALSE, cdict = NULL, showonlypbelow = 0.05, makespreadsheets = FALSE, makeheatmaps = TRUE, maxnumheatmaps = NULL, numthreads = 4, nperm = 99, class_to_ignore = NULL, ...){

    if (is.null(usefulexp)){
        usefulexp <- names(expvec)[!(names(expvec) %in% c("FeatType", "vfdb", "SFLD", "Coils", "Gene3D", "Phobius", "ProSitePatterns", "SMART", "resfinder", "ProDom"))]
    } else {
        usefulexp <- names(expvec)[(names(expvec) %in% usefulexp)]
    }
    expvec2 <- expvec[usefulexp]

    #if (is.null(applyfilters)){
    #    applyfilters <- "none"
    #}

    #Set variables
    variables_disc <- variable_list$discrete
    variables_bin <- variable_list$binary
    variables_cont <- variable_list$continuous
    variables_subs <- variable_list$subsettable
    variables_all <- unique(c(variables_disc, variables_subs, variables_bin, variables_cont))
    variables_discont <- variables_all[!(variables_all %in% variables_cont)]

    pt <- as.data.frame(colData(expvec2[[1]]))

    colcategories_list <- split(variables_all, ceiling(seq_along(variables_all) / max_annot_legends))

    if  (!is.null(class_to_ignore)){
        flog.info(paste("Metadata value(s)", paste0(class_to_ignore, collapse = ", "), "will be omitted." ))
    }

    if (includereaddata == TRUE){
        readscols <- (colnames(pt)[(colnames(pt) %in% c("GbNAHS", "PctAss"))])
        if (length(readscols) > 0){
            #Squeeze in reads into last list position if there is space, else do them separately
            lastpos <- colcategories_list[[length(colcategories_list)]]
            if (length(lastpos) < max_annot_legends){
                colcategories_list[[length(colcategories_list)]] <- append(lastpos, readscols)
            } else {
                colcategories_list[[length(colcategories_list) + 1]] <- readscols
            }
        }
    }

    if (report == "comparative"){
        plot_category_message <- "Relative Abundance Heatmaps"
        plotname <- "Feature_Heatmaps_Comparative"
        statname <- "Feature_Stats_Comparative"
        comparisons <- variables_discont
        plotmsg <- c("Difference in features", "between discrete categories")
        #Careful not to overplot
        if (is.null(showonlypbelow)){
            topcats <- absolutemaxfeats
        } else {
            topcats <- NULL
        }

    } else if (report == "exploratory"){
        plot_category_message <- "Relative Abundance Heatmaps"
        plotname <- "Feature_Heatmaps_Variance"
        statname <- "Feature_Stats_Variance"
        comparisons <- "bRobDingnag"
        plotmsg <- c("Features most variant", "across samples")
        topcats <- absolutemaxfeats
        if (is.null(topcats)){
            topcats <- 100
        }

    } else if (report == "PA"){
        plot_category_message <- "Relative Abundance Heatmaps"
        plotname <- "Feature_Heatmaps_PA_Fisher"
        statname <- "Feature_Stats_PA_Fisher"
        comparisons <- variables_bin
        plotmsg <- c("Features present or absent", "between binary categories")
        #Careful not to overplot
        if (is.null(showonlypbelow)){
            topcats <- absolutemaxfeats
        } else {
            topcats <- NULL
        }

    } else if (report == "correlation"){
        plot_category_message <- "Feature Correlation Heatmaps"
        plotname <- "Feature_Heatmaps_Pairwise_Correlation"
        statname <- "Feature_Heatmaps_Stats_Correlation"
        comparisons <- "bRobDingnag"
        plotmsg <- "Pairwise correlation of features"
    }

    #Set counters
    XLn <- 0
    HMn <- 0

    basepdffn <- paste("JAMS", project, plotname, sep = "_")
    if (!(is.null(appendtofilename))){
        basepdffn <- paste(basepdffn, appendtofilename, sep = "_")
    }
    pdffn <- paste(basepdffn, "pdf", sep=".")
    #Check if fn exists to avoid overwriting.
    ffn <- 1
    while (file.exists(pdffn)){
        pdffn <- paste(paste0(basepdffn, ffn), "pdf", sep=".")
        ffn <- ffn + 1
        if(ffn > 100){
            stop("I think you have enough plots already. May I suggest you try looking through them.")
        }
    }

    if (makeheatmaps == TRUE){
        pdf(pdffn, paper = "a4r")
        authors <- as.character(as.person(packageDescription("JAMS")$Author))
        intromessage <- c(packageDescription("JAMS")$Title, paste("JAMS version", packageVersion("JAMS")), authors, paste("Contact:", "john.mcculloch@nih.gov"), "National Cancer Institute", "National Institutes of Health", "Bethesda, MD, USA")
        plot.new()
        grid.table(intromessage, rows = NULL, cols = NULL, theme = ttheme_default(base_size = 15))
        plot.new()
        grid.table(c("JAMS", plot_category_message, project, plotmsg), rows = NULL, cols = NULL, theme = ttheme_default(base_size = 20))
    }

    #Generate XL files with comparisons throughout samples
    compvec <- NULL
    mycomp <- NULL

    #Loop round comparisons
    for (cmp in comparisons){
        #Loop round experiments
        for (a in 1:length(expvec2)){
            XLn <- XLn+1
            flog.info(paste("Calculating stats for", names(expvec2)[a]))
            if (names(expvec2)[a] == "LKT"){
                featcutoff <- featcutofftaxa
                minl2fc <- l2fchmtaxa
            } else {
                featcutoff <- featcutofffunct
                minl2fc <- l2fchmfunc
            }

            for (annotcnk in 1:length(colcategories_list)){
                colcategories <- colcategories_list[[annotcnk]]
                #Plot heatmap in tandem
                if (report == "comparative"){
                    flog.info(paste("Plotting", names(expvec2)[a], "relabund heatmaps of features most different between", cmp))
                    #check if the comparison is already in colcategories
                    if (!(cmp %in% colcategories)){
                        colcategories <- c(colcategories, cmp)
                    }
                    flog.info(paste("Using legend annotations:", paste0(colcategories, collapse = ", ")))
                    if (annotcnk < 2){
                        plot.new()
                        grid.table(c(names(expvec2)[a], "Features most different", paste("between", cmp), "No subsetting"), rows = NULL, cols = NULL, theme = ttheme_default(base_size = 20))
                    }

                    if (splitcolsbycomp) {
                        splitby <- cmp
                    } else {
                        splitby <- splitcolsby
                    }

                    mycomp <- plot_relabund_heatmap(ExpObj = expvec2[[a]], glomby = glomby, hmtype = report, hmasPA = hmasPA, compareby = cmp, ntop = topcats, invertbinaryorder = FALSE, colcategories = colcategories, applyfilters = applyfilters, minl2fc = minl2fc, featcutoff = featcutoff, samplesToKeep = samplesToKeep, featuresToKeep = featuresToKeep, GenomeCompletenessCutoff = GenomeCompletenessCutoff, PctFromCtgscutoff = PctFromCtgscutoff, ignoreunclassified = TRUE, adjustpval = adjustpval, showonlypbelow = showonlypbelow, showpval = TRUE, showl2fc = TRUE, cdict = cdict, maxnumheatmaps = maxnumheatmaps, numthreads = numthreads, nperm = nperm, returnstats = makespreadsheets, cluster_samples_per_heatmap = cluster_samples_per_heatmap, cluster_features_per_heatmap = cluster_features_per_heatmap, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, class_to_ignore = class_to_ignore, scaled = scaled, no_underscores = no_underscores, threshPA = threshPA, splitcolsby = splitby, ordercolsby = ordercolsby, textby = textby, cluster_rows = cluster_rows, max_rows_in_heatmap = max_rows_in_heatmap, secondaryheatmap = secondaryheatmap)

                    if (annotcnk < 2){
                        compvec <- append(compvec, mycomp, after = length(compvec))
                    }

                } else if (report == "exploratory"){
                    flog.info(paste("Plotting", names(expvec2)[a], "relabund heatmaps of features with highest variance."))
                    flog.info(paste("Using legend annotations:", paste0(colcategories, collapse = ", ")))

                    if (annotcnk < 2){
                        plot.new()
                        grid.table(c("Features with highest variance", names(expvec2)[a], "No subsetting"), rows = NULL, cols = NULL, theme = ttheme_default(base_size = 20))
                    }
                    mycomp <- plot_relabund_heatmap(ExpObj = expvec2[[a]], glomby = glomby, hmtype = report, hmasPA = hmasPA, compareby = cmp, ntop = topcats, colcategories = colcategories, applyfilters = applyfilters, featcutoff = featcutoff, samplesToKeep = samplesToKeep, featuresToKeep = featuresToKeep, GenomeCompletenessCutoff = GenomeCompletenessCutoff, PctFromCtgscutoff = PctFromCtgscutoff, ignoreunclassified = TRUE, cdict = cdict, maxnumheatmaps = maxnumheatmaps, returnstats = makespreadsheets, cluster_samples_per_heatmap = cluster_samples_per_heatmap, cluster_features_per_heatmap = cluster_features_per_heatmap, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, class_to_ignore = class_to_ignore, scaled = scaled, no_underscores = no_underscores, threshPA = threshPA, splitcolsby = splitcolsby, ordercolsby = ordercolsby, textby = textby, cluster_rows = cluster_rows, max_rows_in_heatmap = max_rows_in_heatmap, secondaryheatmap = secondaryheatmap)

                    if (annotcnk < 2){
                        compvec <- append(compvec, mycomp, after=length(compvec))
                    }

                } else if (report == "PA"){
                    flog.info(paste("Plotting", names(expvec2)[a], "relabund heatmaps of features most differentially present between", cmp))
                    #check if the comparison is already in colcategories

                    if (annotcnk < 2){
                        plot.new()
                        grid.table(c(names(expvec2)[a], "Features present or absent", paste("between", cmp), "No subsetting"), rows = NULL, cols = NULL, theme = ttheme_default(base_size = 20))
                    }

                    if (splitcolsbycomp) {
                        splitby <- cmp
                    } else {
                        splitby <- splitcolsby
                    }

                    mycomp <- plot_relabund_heatmap(ExpObj = expvec2[[a]], glomby = glomby, hmtype = report, hmasPA = hmasPA, compareby = cmp, ntop = topcats, invertbinaryorder = FALSE, colcategories = colcategories, applyfilters = applyfilters, minl2fc = minl2fc, featcutoff = featcutoff, samplesToKeep = samplesToKeep, featuresToKeep = featuresToKeep, GenomeCompletenessCutoff = GenomeCompletenessCutoff, PctFromCtgscutoff = PctFromCtgscutoff, ignoreunclassified = TRUE, adjustpval = adjustpval, showonlypbelow = showonlypbelow, showpval = TRUE, showl2fc = TRUE, cdict = cdict, maxnumheatmaps = maxnumheatmaps, numthreads = numthreads, nperm = nperm, returnstats = makespreadsheets, cluster_samples_per_heatmap = cluster_samples_per_heatmap, cluster_features_per_heatmap = cluster_features_per_heatmap, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, class_to_ignore = class_to_ignore, scaled = scaled, no_underscores = no_underscores, threshPA = threshPA, splitcolsby = splitby, ordercolsby = ordercolsby, textby = textby, cluster_rows = cluster_rows, max_rows_in_heatmap = max_rows_in_heatmap, secondaryheatmap = secondaryheatmap)

                    if (annotcnk < 2){
                        compvec <- append(compvec, mycomp, after=length(compvec))
                    }

                } else if (report == "correlation"){
                    if (annotcnk < 2){
                        plot.new()
                        grid.table(c(names(expvec2)[a], "Pairwise Correlation of Features", "No subsetting"), rows = NULL, cols = NULL, theme = ttheme_default(base_size = 20))

                        plot_correlation_heatmap(ExpObj = expvec2[[a]], glomby = glomby, stattype = "spearman", subsetby = NULL, applyfilters = applyfilters, featcutoff = featcutoff, samplesToKeep = samplesToKeep, featuresToKeep = featuresToKeep, GenomeCompletenessCutoff = GenomeCompletenessCutoff, PctFromCtgscutoff = PctFromCtgscutoff, ignoreunclassified = TRUE, ntopvar = 250, minabscorrcoeff = minabscorrcoeff, class_to_ignore = class_to_ignore)
                    }
                    compvec <- NULL
                }
                HMn <- HMn + 1
            }

            if (length(variables_subs) > 0){
                #If there is any subsettable data, then do that.
                validsubs <- NULL
                validsubs <- variables_subs[!(variables_subs %in% cmp)]
                for (vs in validsubs){
                    XLn <- XLn + 1

                    for (annotcnk in 1:length(colcategories_list)){
                        colcategories <- colcategories_list[[annotcnk]]
                        if (report == "comparative"){
                            flog.info(paste("Plotting", names(expvec2)[a], "relabund heatmaps of features most different between", cmp, "within", vs))
                            #check if the comparison is already in colcategories
                            if (!(cmp %in% colcategories)){
                                colcategories <- c(colcategories, cmp)
                            }
                            flog.info(paste("Using legend annotations:", paste0(colcategories, collapse = ", ")))

                            if (annotcnk < 2){
                                plot.new()
                                grid.table(c(names(expvec2)[a], "Features most different", paste("between", cmp), paste("Subset by", vs)), rows = NULL, cols = NULL, theme = ttheme_default(base_size = 20))
                            }

                            if (splitcolsbycomp) {
                                splitby <- cmp
                            } else {
                                splitby <- splitcolsby
                            }

                            mycomp <- plot_relabund_heatmap(ExpObj = expvec2[[a]], glomby = glomby, subsetby = vs, hmtype = report, hmasPA = hmasPA, compareby = cmp, ntop = topcats, invertbinaryorder = FALSE, colcategories = colcategories, applyfilters = applyfilters, minl2fc = minl2fc, featcutoff = featcutoff, samplesToKeep = samplesToKeep, featuresToKeep = featuresToKeep, GenomeCompletenessCutoff = GenomeCompletenessCutoff, PctFromCtgscutoff = PctFromCtgscutoff, ignoreunclassified = TRUE, adjustpval = adjustpval, showonlypbelow = showonlypbelow, showpval = TRUE, showl2fc = TRUE, cdict = cdict, maxnumheatmaps = maxnumheatmaps, numthreads = numthreads, nperm = nperm, returnstats = makespreadsheets, cluster_samples_per_heatmap = cluster_samples_per_heatmap, cluster_features_per_heatmap = cluster_features_per_heatmap, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, class_to_ignore = class_to_ignore, scaled = scaled, no_underscores = no_underscores, threshPA = threshPA, splitcolsby = splitby, ordercolsby = ordercolsby, textby = textby, cluster_rows = cluster_rows, max_rows_in_heatmap = max_rows_in_heatmap, secondaryheatmap = secondaryheatmap)

                            if (annotcnk < 2){
                                compvec <- append(compvec, mycomp, after = length(compvec))
                            }

                        } else if (report == "exploratory"){
                            #Explain what we are doing
                            flog.info(paste("Plotting", names(expvec2)[a], "relabund heatmaps of features with higest variance within", vs))
                            flog.info(paste("Using legend annotations:", paste0(colcategories, collapse = ", ")))

                            if (annotcnk < 2){
                                plot.new()
                                grid.table(c("Features with highest variance", names(expvec2)[a], "Subsetting by group"), rows = NULL, cols = NULL, theme = ttheme_default(base_size = 20))
                            }

                            mycomp <- plot_relabund_heatmap(ExpObj = expvec2[[a]], glomby = glomby, subsetby = vs, hmtype = report, hmasPA = hmasPA, compareby = cmp, ntop = topcats, colcategories = colcategories, applyfilters = applyfilters, featcutoff = featcutoff, samplesToKeep = samplesToKeep, featuresToKeep = featuresToKeep, GenomeCompletenessCutoff = GenomeCompletenessCutoff, PctFromCtgscutoff = PctFromCtgscutoff, ignoreunclassified = TRUE, cdict = cdict, maxnumheatmaps = maxnumheatmaps, returnstats = makespreadsheets, cluster_samples_per_heatmap = cluster_samples_per_heatmap, cluster_features_per_heatmap = cluster_features_per_heatmap, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, class_to_ignore = class_to_ignore, scaled = scaled, no_underscores = no_underscores, threshPA = threshPA, splitcolsby = splitcolsby, ordercolsby = ordercolsby, textby = textby, cluster_rows = cluster_rows, max_rows_in_heatmap = max_rows_in_heatmap, secondaryheatmap = secondaryheatmap)

                            if (annotcnk < 2){
                                compvec <- append(compvec, mycomp, after = length(compvec))
                            }

                        } else if (report == "PA"){
                            flog.info(paste("Plotting", names(expvec2)[a], "relabund heatmaps of features most differentially present between", cmp, "within", vs))
                            if (annotcnk < 2){
                                plot.new()
                                grid.table(c(names(expvec2)[a], "Features present or absent", paste("between", cmp), paste("Subset by", vs)), rows = NULL, cols = NULL, theme = ttheme_default(base_size = 20))
                            }

                            if (splitcolsbycomp) {
                                splitby <- cmp
                            } else {
                                splitby <- splitcolsby
                            }

                            mycomp <- plot_relabund_heatmap(ExpObj = expvec2[[a]], glomby = glomby, subsetby = vs, hmtype = report, hmasPA = hmasPA, compareby = cmp, ntop = topcats, invertbinaryorder = FALSE, colcategories = colcategories, applyfilters = applyfilters, minl2fc = minl2fc, featcutoff = featcutoff, samplesToKeep = samplesToKeep, featuresToKeep = featuresToKeep, GenomeCompletenessCutoff = GenomeCompletenessCutoff, PctFromCtgscutoff = PctFromCtgscutoff, ignoreunclassified = TRUE, adjustpval = adjustpval, showonlypbelow = showonlypbelow, showpval = TRUE, showl2fc = TRUE, cdict = cdict, maxnumheatmaps = maxnumheatmaps, numthreads = numthreads, nperm = nperm, returnstats = makespreadsheets, cluster_samples_per_heatmap = cluster_samples_per_heatmap, cluster_features_per_heatmap = cluster_features_per_heatmap, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, class_to_ignore = class_to_ignore, scaled = scaled, no_underscores = no_underscores, threshPA = threshPA, splitcolsby = splitby, ordercolsby = ordercolsby, textby = textby, cluster_rows = cluster_rows, max_rows_in_heatmap = max_rows_in_heatmap, secondaryheatmap = secondaryheatmap)

                            if (annotcnk < 2){
                                compvec <- append(compvec, mycomp, after = length(compvec))
                            }

                        } else if (report == "correlation"){
                            if (annotcnk < 2){
                                plot.new()
                                grid.table(c(names(expvec2)[a], "Pairwise Correlation of Features", paste("Subset by", vs)), rows = NULL, cols = NULL, theme = ttheme_default(base_size = 20))

                                plot_correlation_heatmap(ExpObj = expvec2[[a]], glomby = glomby, stattype = "spearman", subsetby = vs, applyfilters = applyfilters, featcutoff = featcutoff, samplesToKeep = samplesToKeep, featuresToKeep = featuresToKeep, GenomeCompletenessCutoff = GenomeCompletenessCutoff, PctFromCtgscutoff = PctFromCtgscutoff, ignoreunclassified = TRUE, ntopvar = 250, minabscorrcoeff = minabscorrcoeff, class_to_ignore = class_to_ignore)
                            }
                            compvec <- NULL
                        } #Type of report conditional
                        HMn <- HMn + 1
                    } #End loop for colcategories chunks
                } #End of for each valid subset loop
            } #If there are any subsets conditional
        } #End of for each experiment loop
    } #End of for each comparison loop

    if (makeheatmaps == TRUE){
        #Final message
        authors = as.character(as.person(packageDescription("JAMS")$Author))
        finalmessage = c(packageDescription("JAMS")$Title, paste("JAMS version", packageVersion("JAMS")), authors, paste("Contact:", "john.mcculloch@nih.gov"), "National Cancer Institute", "National Institutes of Health", "Bethesda, MD, USA")
        plot.new()
        grid.table(finalmessage, rows = NULL, cols = NULL, theme = ttheme_default(base_size = 15))
        dev.off()
    }

    if (makespreadsheets == TRUE){
        #Output statistics to system as spreadsheet.
        #flush out any null dfs in list
        compvec <- compvec[sapply(compvec, function(x){!(is.null(x))})]

        if (length(compvec) > 1){
            #Truncate names to 31 chars. Excel does not accept more. Sigh.
            compvecXL <- compvec
            names(compvecXL) <- sapply(1:length(compvecXL), function(x){ stringr::str_trunc(names(compvecXL)[x], 30) })

            baseXLfn <- paste("JAMS", project, statname, sep="_")
            if (!(is.null(appendtofilename))){
                baseXLfn <- paste(baseXLfn, appendtofilename, sep = "_")
            }
            XLfn <- paste(baseXLfn, "xlsx", sep=".")
            #Ensure no overwriting.
            xfn <- 1
            while (file.exists(XLfn)){
                XLfn <- paste(paste0(baseXLfn, xfn), "xlsx", sep = ".")
                xfn <- xfn + 1
                if(xfn > 100){
                    stop("I think you have enough tables already. May I suggest you try looking through them.")
                }
            }
            write.xlsx(compvecXL, file = XLfn, asTable = TRUE, rowNames = TRUE, colNames = TRUE, borders = "surrounding", colWidths = "auto")
        }
    }

    return(flog.info(paste("Comparative report complete with", XLn, "spreadsheet comparisons and", HMn, "heatmap comparisons")))
}
