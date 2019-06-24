#' make_heatmap_report(report=c("comparative", "exploratory", "PA"), project=NULL, expvec=NULL, usefulexp=NULL, featcutofftaxa=NULL, featcutofffunct=NULL, l2fchmtaxa=NULL, l2fchmfunc=NULL, applyfilters=c("stringent", "moderate", "none"), absolutemaxfeats=NULL, genomecompleteness=NULL, variable_list=NULL, adjustpval=TRUE, list.data=NULL, asPA=FALSE, mgSeqnorm=FALSE, cdict=NULL, showonlypbelow=NULL, makespreadsheets=FALSE, makeheatmaps=TRUE, maxnumheatmaps=NULL, ...)
#'
#' Generates standard comparative heatmaps and spreadsheet based on relative abundance difference between binary categories for analyses named in expvectors.
#' @export

make_heatmap_report <- function(report = "comparative", project = NULL, expvec = NULL, usefulexp = NULL, appendtofilename = NULL, featcutofftaxa = NULL, featcutofffunct = NULL, l2fchmtaxa = NULL, l2fchmfunc = NULL, samplesToKeep = NULL, featuresToKeep = NULL, applyfilters = NULL, absolutemaxfeats = NULL, genomecompleteness = NULL, variable_list = NULL, adjustpval = NULL, list.data = NULL, hmasPA = FALSE, mgSeqnorm = FALSE, cdict = NULL, showonlypbelow = 0.05, makespreadsheets = FALSE, makeheatmaps = TRUE, maxnumheatmaps = NULL, numthreads = 4, nperm = 99, ...){

    if (is.null(usefulexp)){
        usefulexp <- names(expvec)[!(names(expvec) %in% c("FeatType", "vfdb", "SFLD", "Coils", "Gene3D", "Phobius", "ProSitePatterns", "SMART", "resfinder", "ProDom"))]
    } else {
        usefulexp <- names(expvec)[(names(expvec) %in% usefulexp)]
    }
    expvec2 <- expvec[usefulexp]

    if (is.null(applyfilters)){
        applyfilters <- "none"
    }

    #Set variables
    variables_disc <- variable_list$discrete
    variables_bin <- variable_list$binary
    variables_cont <- variable_list$continuous
    variables_subs <- variable_list$subsettable
    variables_all <- unique(c(variables_disc, variables_subs, variables_bin, variables_cont))
    variables_discont <- variables_all[!(variables_all %in% variables_cont)]

    pt <- as.data.frame(pData(expvec2[[1]]))
    readscols <- (colnames(pt)[(colnames(pt) %in% c("GbNAHS", "PctAss"))])
    colcategories <- append(variables_all, readscols, after=length(variables_all))

    if (report == "comparative"){
        plotname <- "Feature_Heatmaps_Comparative"
        statname <- "Feature_Stats_Comparative"
        comparisons <- variables_discont
        plotmsg <- c("Difference in features", "between discrete categories")
    } else if (report == "exploratory"){
        plotname <- "Feature_Heatmaps_Variance"
        statname <- "Feature_Stats_Variance"
        comparisons <- "bRobDingnag"
        plotmsg <- c("Features most variant", "across samples")
    } else if (report == "PA"){
        plotname <- "Feature_Heatmaps_PA_Fisher"
        statname <- "Feature_Stats_PA_Fisher"
        comparisons <- variables_bin
        plotmsg <- c("Features present or absent", "between binary categories")
    }

    #Set counters
    XLn <- 0
    HMn <- 0

    basepdffn <- paste("JAMS", project, plotname, sep="_")
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
        grid.table(c("JAMS", "Relative Abundance Heatmaps", project, plotmsg), rows = NULL, cols = NULL, theme = ttheme_default(base_size = 20))
    }

    #Generate XL files with comparisons throughout samples
    compvec <- NULL
    mycomp <- NULL

    if (is.null(showonlypbelow)){
        topcats <- absolutemaxfeats
    } else {
        topcats <- NULL
    }

    #Loop round comparisons
    for (cmp in comparisons){
        #Loop round experiments
        for (a in 1:length(expvec2)){
            XLn <- XLn+1
            flog.info(paste("Calculating stats for", names(expvec2)[a]))
            if (makeheatmaps == TRUE){
                #Plot heatmap in tandem
                if (report == "comparative"){
                    flog.info(paste("Plotting", names(expvec2)[a], "relabund heatmaps of features most different between", cmp))
                    plot.new()
                    grid.table(c(names(expvec2)[a], "Features most different", paste("between", cmp), "No subsetting"), rows = NULL, cols = NULL, theme = ttheme_default(base_size = 20))

                    #Verify if there are enough samples in each category to go ahead

                    mycomp <- plot_relabund_heatmap(mgseqobj=expvec2[[a]], heatpalette="diverging", hmtype=report, hmasPA=FALSE, compareby=cmp, ntop=topcats, invertbinaryorder=FALSE, colcategories=colcategories, cluster_rows=TRUE, applyfilters=applyfilters, samplesToKeep=samplesToKeep, featuresToKeep=featuresToKeep, adjustpval=adjustpval, showonlypbelow=showonlypbelow, showpval=TRUE, showl2fc=TRUE, list.data=list.data, cdict=cdict, maxnumheatmaps=maxnumheatmaps, numthreads=numthreads, nperm=nperm, statsonlog=TRUE, returnstats=makespreadsheets)

                    compvec <- append(compvec, mycomp, after=length(compvec))

                } else if (report == "exploratory"){
                    flog.info(paste("Plotting", names(expvec2)[a], "relabund heatmaps of features with higest variance."))
                    plot.new()
                    grid.table(c("Features with highest variance", names(expvec2)[a], "No subsetting"), rows = NULL, cols = NULL, theme = ttheme_default(base_size = 20))

                    mycomp<-plot_relabund_heatmap(mgseqobj=expvec2[[a]],  heatpalette="diverging", hmtype=report, ntop=topcats, colcategories=colcategories, cluster_rows=TRUE, applyfilters=applyfilters, samplesToKeep=samplesToKeep, featuresToKeep=featuresToKeep, list.data=list.data, cdict=cdict, maxnumheatmaps=maxnumheatmaps, numthreads=numthreads, nperm=nperm, statsonlog=TRUE, returnstats=makespreadsheets)

                    compvec <- append(compvec, mycomp, after=length(compvec))

                } else if (report == "PA"){
                    plot.new()
                    grid.table(c(names(expvec2)[a], "Features present or absent", paste("between", cmp), "No subsetting"), rows = NULL, cols = NULL, theme = ttheme_default(base_size = 20))

                    mycomp<-plot_relabund_heatmap(mgseqobj=expvec2[[a]],  heatpalette="diverging", hmtype=report, hmasPA=TRUE, compareby=cmp, ntop=topcats, invertbinaryorder=FALSE, colcategories=colcategories, cluster_rows=TRUE, applyfilters=applyfilters, samplesToKeep=samplesToKeep, featuresToKeep=featuresToKeep, adjustpval=adjustpval, showonlypbelow=showonlypbelow, showpval=TRUE, showl2fc=TRUE, list.data=list.data, cdict=cdict, maxnumheatmaps=maxnumheatmaps, numthreads=numthreads, nperm=nperm, statsonlog=TRUE, returnstats=makespreadsheets)

                    compvec <- append(compvec, mycomp, after=length(compvec))
                }
            }
            HMn=HMn+1

            if(length(variables_subs) > 0){
                #If there is any subsettable data, then do that.
                validsubs<-NULL
                validsubs<-variables_subs[!(variables_subs %in% cmp)]
                for (vs in validsubs){
                    XLn=XLn+1

                    if(makeheatmaps==TRUE){
                        if(report == "comparative"){
                            plot.new()
                            grid.table(c(names(expvec2)[a], "Features most different", paste("between", cmp), paste("Subset by", vs)), rows = NULL, cols = NULL, theme = ttheme_default(base_size = 20))

                            mycomp<-plot_relabund_heatmap(mgseqobj=expvec2[[a]], subsetby=vs, heatpalette="diverging", hmtype=report, hmasPA=FALSE, compareby=cmp, ntop=topcats, invertbinaryorder=FALSE, colcategories=colcategories, cluster_rows=TRUE, applyfilters=applyfilters, adjustpval=adjustpval, samplesToKeep=samplesToKeep, featuresToKeep=featuresToKeep, showonlypbelow=showonlypbelow, showpval=TRUE, showl2fc=TRUE, list.data=list.data, cdict=cdict, maxnumheatmaps=maxnumheatmaps, numthreads=numthreads, nperm=nperm, statsonlog=TRUE, returnstats=makespreadsheets)

                            compvec <- append(compvec, mycomp, after=length(compvec))

                        } else if (report == "exploratory"){
                            #Explain what we are doing
                            plot.new()
                            grid.table(c("Features with highest variance", names(expvec2)[a], "Subsetting by group"), rows = NULL, cols = NULL, theme = ttheme_default(base_size = 20))

                            mycomp<-plot_relabund_heatmap(mgseqobj=expvec2[[a]], subsetby=vs, heatpalette="diverging", hmtype=report, ntop=topcats, colcategories=colcategories, cluster_rows=TRUE, applyfilters=applyfilters, samplesToKeep=samplesToKeep, featuresToKeep=featuresToKeep, list.data=list.data, cdict=cdict, maxnumheatmaps=maxnumheatmaps, numthreads=numthreads, nperm=nperm, statsonlog=TRUE, returnstats=makespreadsheets)

                            compvec <- append(compvec, mycomp, after=length(compvec))

                        } else if (report == "PA"){
                            plot.new()
                            grid.table(c(names(expvec2)[a], "Features present or absent", paste("between", cmp), paste("Subset by", vs)), rows = NULL, cols = NULL, theme = ttheme_default(base_size = 20))

                            mycomp<-plot_relabund_heatmap(mgseqobj=expvec2[[a]], subsetby=vs, heatpalette="diverging", hmtype=report, hmasPA=TRUE, compareby=cmp, ntop=topcats, invertbinaryorder=FALSE, colcategories=colcategories, cluster_rows=TRUE, applyfilters=applyfilters, samplesToKeep=samplesToKeep, featuresToKeep=featuresToKeep, adjustpval=adjustpval, showonlypbelow=showonlypbelow, showpval=TRUE, showl2fc=TRUE, list.data=list.data, cdict=cdict, maxnumheatmaps=maxnumheatmaps, numthreads=numthreads, nperm=nperm, statsonlog=TRUE, returnstats=makespreadsheets)

                            compvec <- append(compvec, mycomp, after=length(compvec))

                        } #Type of report conditional
                    } #If heatmaps are to be made conditional
                    HMn=HMn+1
                } #End of for each valid subset loop
            } #If there are any subsets conditional
        } #End of for each experiment loop
    } #End of for each comparison loop

    if(makeheatmaps==TRUE){
        #Final message
        authors=as.character(as.person(packageDescription("JAMS")$Author))
        finalmessage=c(packageDescription("JAMS")$Title, paste("JAMS version", packageVersion("JAMS")), authors, paste("Contact:", "john.mcculloch@nih.gov"), "National Cancer Institute", "National Institutes of Health", "Bethesda, MD, USA")
        plot.new()
        grid.table(finalmessage, rows = NULL, cols = NULL, theme = ttheme_default(base_size = 15))
        dev.off()
    }

    if(makespreadsheets==TRUE){
        #Output statistics to system as spreadsheet.
        #flush out any null dfs in list
        compvec<-compvec[sapply(compvec, function(x){!(is.null(x))})]
        #Truncate names to 31 chars. Excel does not accept more. Sigh.
        compvecXL<-compvec
        names(compvecXL)<-sapply(1:length(compvecXL), function(x){ stringr::str_trunc(names(compvecXL)[x], 30) })

        baseXLfn <- paste("JAMS", project, statname, sep="_")
        if (!(is.null(appendtofilename))){
            baseXLfn <- paste(baseXLfn, appendtofilename, sep = "_")
        }
        XLfn <- paste(baseXLfn, "xlsx", sep=".")
        #Ensure no overwriting.
        xfn <- 1
        while (file.exists(XLfn)){
            XLfn <- paste(paste0(baseXLfn, xfn), "xlsx", sep=".")
            xfn <- xfn + 1
            if(xfn > 100){
                stop("I think you have enough tables already. May I suggest you try looking through them.")
            }
        }

        write.xlsx(compvecXL, file = XLfn, asTable = TRUE, rowNames = TRUE, colNames = TRUE, borders = "surrounding", colWidths="auto")
    }

    return(flog.info(paste("Comparative report complete with", XLn, "spreadsheet comparisons and", HMn, "heatmap comparisons")))
}
