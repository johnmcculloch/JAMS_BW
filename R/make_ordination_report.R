#' make_ordination_report(project = NULL, algorithm = "PCA", expvec = NULL, usefulexp = NULL, mgSeqnorm = FALSE, genomecompleteness = NULL, variable_list = NULL, list.data = NULL, doreads = NULL, cdict = NULL, ellipse = "auto", samplesToKeep = NULL, featuresToKeep = NULL, ignoreunclassified=TRUE, appendtofilename = NULL, ...)
#'
#' Generates standard ordination plots for analyses named in expobjects.
#' @export

make_ordination_report <- function(project = NULL, algorithm = "PCA", expvec = NULL, usefulexp = NULL, mgSeqnorm = FALSE, genomecompleteness = NULL, variable_list = NULL, list.data = NULL, doreads = NULL, cdict = NULL, ellipse = "auto", samplesToKeep = NULL, featuresToKeep = NULL, ignoreunclassified = TRUE, appendtofilename = NULL, ...){

    if(is.null(usefulexp)){
        usefulexp <- names(expvec)[!(names(expvec) %in% c("SFLD", "Coils", "Gene3D", "Phobius", "ProSitePatterns", "SMART", "ProDom"))]
    } else {
        usefulexp <- names(expvec)[(names(expvec) %in% usefulexp)]
    }
    expvec2 <- expvec[usefulexp]

    #Set variables
    variables_disc <- variable_list$discrete
    variables_bin <- variable_list$binary
    variables_cont <- variable_list$continuous
    variables_subs <- variable_list$subsettable

    pt <- as.data.frame(pData(expvec2[[1]]))
    variables_all <- unique(c(variables_disc, variables_subs, variables_bin, variables_cont))
    readscols <- (colnames(pt)[(colnames(pt) %in% c("GbNAHS", "PctAss"))])
    if (length(readscols) > 0){
        if (is.null(doreads)){
            doreads=TRUE
        }
    }

    #Plot ordination plots in tandem to saving dataframes with stats.
    basepdffn <- paste("JAMS", project, "Ordination_plots", sep="_")
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

    pdf(pdffn, paper="a4r")

    plot.new()
    grid.table(c("JAMS", "Ordination plots", project), rows = NULL, cols = NULL, theme = ttheme_default(base_size = 25))

    #If project includes read data, plot sequencing depth information.
    if(doreads == TRUE){
        flog.info("Plotting ordination and marking by sequencing depth.")
        plot.new()
        grid.table(c("Check for sequencing depth bias"), rows = NULL, cols = NULL, theme = ttheme_default(base_size = 25))

        #Cycle through possible experiments
        for(e in 1:length(expvec2)){
            print(plot_Ordination(mgseqobj = expvec2[[e]], mgSeqnorm=mgSeqnorm, algorithm = algorithm, colourby = "GbNAHS", logtran = TRUE, transp = TRUE, permanova = FALSE, ellipse = FALSE, samplesToKeep=samplesToKeep, featuresToKeep=featuresToKeep))
        }
    }

    variables_all <- unique(c(variables_disc, variables_subs, variables_bin, variables_cont))
    variables_discontinuous <- variables_all[!(variables_all %in% variables_cont)]

    #plot ordination if there are any plottable variables
    if(length(variables_all) > 0){
        flog.info("Plotting ordination and marking by available variables.")
        plot.new()
        grid.table(c("Marking samples by metadata variables"), rows = NULL, cols = NULL, theme = ttheme_default(base_size = 25))
        for(c in 1:length(variables_all)){
            plot.new()
            grid.table(c("Marking samples by", variables_all[c]), rows = NULL, cols = NULL, theme = ttheme_default(base_size = 25))

            for(e in 1:length(expvec2)){
                flog.info(paste("Plotting", names(expvec2)[e], " and marking by", variables_all[c]))

                print(plot_Ordination(mgseqobj = expvec2[[e]], mgSeqnorm=mgSeqnorm, algorithm = algorithm, colourby = variables_all[c], shapeby=NULL, log2tran = TRUE, transp = TRUE, permanova = TRUE, ellipse = ellipse, cdict=cdict, samplesToKeep=samplesToKeep, featuresToKeep=featuresToKeep, ignoreunclassified = ignoreunclassified))
                #If there are any subsettable variables, subset by them.
                validsubs <- variables_subs[!(variables_subs %in% variables_all[c])]

                if (length(validsubs)>0){
                    for (vs in validsubs){
                        #Test if there is more than one class within the subset
                        discnamesvs <- unique(pt[, vs])
                        classtest <- NULL
                        for(k in 1:length(discnamesvs)){
                            #check if there is more than one class when subset
                            classtest[k] <- (length(unique(pt[, variables_all[c]][which(pt[, vs] == discnamesvs[k])])) > 1)
                        }

                        #Plot if there is more than one class in a subset, else skip, as it is pointless.
                        if (all(classtest)){
                            plot.new()
                            grid.table(c("Marking samples by", variables_all[c], "Subsetting samples by", vs), rows = NULL, cols = NULL, theme = ttheme_default(base_size = 25))
                            print(plot_Ordination(mgseqobj = expvec2[[e]], mgSeqnorm=mgSeqnorm, algorithm = algorithm, subsetby=vs, colourby = variables_all[c], log2tran = TRUE, transp = TRUE, permanova = TRUE, ellipse = ellipse, cdict=cdict, samplesToKeep=samplesToKeep, featuresToKeep=featuresToKeep, ignoreunclassified=ignoreunclassified))
                        } else {
                            flog.info(paste("Will not plot", variables_all[c], "within", vs, "because",  paste0(discnamesvs[!classtest], collapse=", " ),  "do(es) not have more than a single class."))
                        } #End conditional that there is more than a single class
                    } #End loop for each subset that is valid
                } #End conditional that there is at least one valid subset
            } #End loop for each experiment
        } #End loop for each variable to mark by
    } #End conditional if there is more than one variable to mark by

    dev.off()

    return(flog.info("Ordination report complete."))
}
