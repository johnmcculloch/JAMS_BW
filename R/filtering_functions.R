#' filter_experiment(ExpObj = NULL, featmaxatleastPPM = 0, featcutoff = c(0, 0), discard_SDoverMean_below = NULL, samplesToKeep = NULL, featuresToKeep = NULL, asPPM = TRUE, PPM_normalize_to_bases_sequenced = FALSE, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL)
#'
#' Filters a SummarizedExperiment object by several criteria.
#' @export

filter_experiment <- function(ExpObj = NULL, featmaxatleastPPM = 0, featcutoff = c(0, 0), discard_SDoverMean_below = NULL, samplesToKeep = NULL, featuresToKeep = NULL, asPPM = TRUE, PPM_normalize_to_bases_sequenced = FALSE, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL){

    #Get only samples you asked for
    if (!(is.null(samplesToKeep))){
        samplesToKeep <- samplesToKeep[samplesToKeep %in% colnames(ExpObj)]
        ExpObj <- ExpObj[, samplesToKeep]
    }

    #Get appropriate object with which to work
    if (as.character(class(ExpObj)[1]) == "SummarizedExperiment"){
        rawcts <- as.matrix(assays(ExpObj)$BaseCounts)
        if (PPM_normalize_to_bases_sequenced) {
            totbases <- metadata(ExpObj)$TotalBasesSequenced
        } else {
            totbases <- metadata(ExpObj)$TotalBasesSequencedinAnalysis
            #totbases <- t(data.frame(NumBases = colSums(rawcts)))
        }
    } else {
        stop("Object must be a SummarizedExperiment object.")
    }

    #If setting featmaxatleastPPM or featcutoff to anything other than the defaults, then return as PPM.
    if ((featmaxatleastPPM != 0) || featcutoff != c(0, 0)){
        asPPM <- TRUE
    }

    #Flush out empty rows
    ExpObj <- ExpObj[(rowSums(rawcts) > 0), ]
    #Flush out empty Samples
    emptysamples <- names(which(colSums(rawcts) == 0) == TRUE)
    if (length(emptysamples) > 0){
        flog.info(paste("Samples", paste0(emptysamples, collapse = ", "), "are empty and will be discarded."))
        validsamples <- names(which(colSums(rawcts) > 0) == TRUE)
        ExpObj <- ExpObj[ , validsamples]
    }

    countmat <- as.matrix(assays(ExpObj)$BaseCounts)
    countmat <- countmat[!is.na(row.names(countmat)),]

    #Transform to PPM if applicable
    if (asPPM == TRUE){
        #transform into PPM
        getPPM <- function(Sample = NULL){
            PPMs <- (countmat[ , Sample] / totbases["NumBases", Sample]) * 1000000
            return(PPMs)
        }

        countmat2 <- sapply(1:ncol(countmat), function(x){ getPPM(Sample = colnames(countmat)[x])} )
        colnames(countmat2) <- colnames(countmat)
        #Round to integer
        for (c in 1:ncol(countmat2)){
            countmat2[, c] <- round(countmat2[, c], 0)
        }
        countmatrix <- countmat2

        #Regenerate an experiment object
        assays(ExpObj)$BaseCounts <- countmatrix

    } else {
        countmatrix <- countmat
    }

    #Discard features which do not match certain criteria
    if (!(missing(featcutoff))){
        if(length(featcutoff) != 2){
            stop("Please specify the minimum PPM in what percentage of the samples you want with a numerical vector of size two. For example, featcutoff=c(2000,10) would discard features which are not at least 2000 PPM in at least 10% of samples.")
        }
        thresholdPPM <- featcutoff[1]
        sampcutoffpct <- min(featcutoff[2], 100)

        prop_above_threshold <- function(x) {
            proportionPassing <- length(which(countmatrix[x, ] >= thresholdPPM)) / ncol(countmatrix)
            return(proportionPassing)
        }

        featuresToKeep2 <- rownames(countmatrix)[which((sapply(1:nrow(countmatrix), function(x) { prop_above_threshold(x) })) > (sampcutoffpct/100))]

        ExpObj <- ExpObj[featuresToKeep2, ]
        countmatrix <- countmatrix[featuresToKeep2, ]

        cutoffmsg <- paste("Feature must be >", thresholdPPM, "PPM in at least", sampcutoffpct, "% of samples", sep=" ")
        flog.info(cutoffmsg)
    }

    if (featmaxatleastPPM > 0){
        featuresToKeep2 <- rownames(countmatrix)[which(rowMax(countmatrix) >= featmaxatleastPPM)]
        ExpObj <- ExpObj[featuresToKeep2, ]
        countmatrix <- countmatrix[featuresToKeep2, ]
        minPPMmsg <- paste("Highest feature must be >", featmaxatleastPPM, "PPM", sep=" ")
        flog.info(minPPMmsg)
    }

    if (!is.null(discard_SDoverMean_below)){
        dfm <- (rowSds(countmatrix) / rowMeans(countmatrix))
        featuresToKeepSDflt <- names(dfm[dfm > discard_SDoverMean_below])
        ExpObj <- ExpObj[featuresToKeepSDflt, ]
        countmatrix <- countmatrix[featuresToKeepSDflt, ]
        SDoverMeanmsg <- paste("Feature must have >", discard_SDoverMean_below, "SDs over mean", sep=" ")
        flog.info(SDoverMeanmsg)
    }

    if (all(c(("PctFromCtgs" %in% names(assays(ExpObj))), (!is.null(PctFromCtgscutoff))))){
        if(length(PctFromCtgscutoff) != 2){
            stop("Please specify the minimum percentage of information coming from contigs in what percentage of the samples you want with a numerical vector of size two. For example, PctFromCtgscutoff = c(90, 10) would discard features whose taxonomic information does not come at least 90% from contigs in at least 10% of samples.")
        }
        minPctFromCtgs <- PctFromCtgscutoff[1]
        sampcutoffpct <- min(PctFromCtgscutoff[2], 100)
        PctFromCtgsmatrix <- assays(ExpObj)[["PctFromCtgs"]]

        prop_above_minCtg <- function(x) {
            validsamples <- names(which(countmatrix[x, ] > 0))
            proportionPassing <- length(which(PctFromCtgsmatrix[x, validsamples] >= minPctFromCtgs)) / length(validsamples)
            return(proportionPassing)
        }

        featuresToKeep3 <- rownames(PctFromCtgsmatrix)[which((sapply(1:nrow(PctFromCtgsmatrix), function(x) { prop_above_minCtg(x) })) > (sampcutoffpct / 100))]

        ExpObj <- ExpObj[featuresToKeep3, ]
        countmatrix <- countmatrix[featuresToKeep3, ]
    }

    if (all(c(("GenomeCompleteness" %in% names(assays(ExpObj))), (!is.null(GenomeCompletenessCutoff))))){
        if(length(GenomeCompletenessCutoff) != 2){
            stop("Please specify the minimum Genome Completeness in what percentage of the samples you want with a numerical vector of size two. For example, GenomeCompletenessCutoff = c(10, 5) would discard features whose genome completeness is smaller than 10% from contigs in at least 5% of samples.")
        }
        minGenomeCompleteness <- (GenomeCompletenessCutoff[1] / 100)
        sampcutoffpct <- min(GenomeCompletenessCutoff[2], 100)
        GenomeCompletenessmatrix <- assays(ExpObj)[["GenomeCompleteness"]]

        prop_above_minGComp <- function(x) {
            validsamples <- names(which(countmatrix[x, ] > 0))
            proportionPassing <- length(which(GenomeCompletenessmatrix[x, validsamples] >= minGenomeCompleteness)) / length(validsamples)
            return(proportionPassing)
        }

        featuresToKeep4 <- rownames(GenomeCompletenessmatrix)[which((sapply(1:nrow(GenomeCompletenessmatrix), function(x) { prop_above_minGComp(x) })) > (sampcutoffpct / 100))]

        ExpObj <- ExpObj[featuresToKeep4, ]
        countmatrix <- countmatrix[featuresToKeep4, ]
    }

    #Get only features you asked for
    if (!(is.null(featuresToKeep))){
        featurespresent <- rownames(countmatrix)
        wantedfeatures <- featuresToKeep[featuresToKeep %in% featurespresent]
        #Check that you still have the features you're asking for
        if (length(wantedfeatures) > 0){
            ExpObj <- ExpObj[wantedfeatures, ]
        } else {
            flog.warn("There are no features you requested surviving the criteria.")
        }
    }

    return(ExpObj)
}


#' ExpObjVetting(ExpObj = NULL)
#'
#' Performs vetting of a SummarizedExperiment object for use in several functions
#' @export

ExpObjVetting <- function(ExpObj = NULL, samplesToKeep = NULL, featuresToKeep = NULL, glomby = NULL, variables_to_fix = NULL, class_to_ignore = NULL){

        #Get appropriate object to work with
        if (as.character(class(ExpObj)[1]) != "SummarizedExperiment"){
            stop("This function can only take a SummarizedExperiment object as input.")
        }
        obj <- ExpObj

        if (!(is.null(glomby))){
            obj <- agglomerate_features(ExpObj = obj, glomby = glomby)
        }

        #Exclude samples and features if specified
        if (!(is.null(samplesToKeep))){
            samplesToKeep <- samplesToKeep[samplesToKeep %in% colnames(obj)]
            obj <- obj[, samplesToKeep]
        }

        if (!(is.null(featuresToKeep))){
            featuresToKeep <- featuresToKeep[featuresToKeep %in% rownames(obj)]
            obj <- obj[featuresToKeep, ]
        }

        obj <- suppressWarnings(filter_sample_by_class_to_ignore(SEobj = obj, variables = variables_to_fix, class_to_ignore = class_to_ignore))

    return(obj)
}


#' filter_sample_by_class_to_ignore(SEobj = NULL, variables = NULL, class_to_ignore = NULL)
#'
#' Filters a SummarizedExperiment object by several criteria.
#' @export

filter_sample_by_class_to_ignore <- function(SEobj = NULL, variables = NULL, class_to_ignore = NULL){

    ptb <- as.data.frame(colData(SEobj))
    allmetadata <- metadata(SEobj)
    Samples <- rownames(ptb)

    if (!is.null(class_to_ignore)){
        #Start off with all of them
        valid_samples <- rownames(ptb)
        variables <- colnames(ptb)[colnames(ptb) %in% variables]
        #Keep omitting samples which do not fit the criteria, as long as variables are contained in current metadata
        if (all(c((!is.null(variables)), (length(variables) > 0)))){
            for (v in 1:length(variables)){
                valid_samples <-  valid_samples[valid_samples %in% (rownames(ptb)[!(ptb[ , variables[v]] %in% class_to_ignore)])]
            }
        }
        if (length(valid_samples) < 1){
            flog.warn("There are no samples matching the criteria. Returning original object with all samples")
        } else {
            omitted_samples <- rownames(ptb)[!(rownames(ptb) %in% valid_samples)]
            flog.info(paste("A total of", length(omitted_samples), "samples were omitted for containing", paste0(class_to_ignore, collapse = ", "), "within metadata variables", paste0(variables, collapse = ", ")))
            SEobj <- SEobj[ , valid_samples]
            ptb <- as.data.frame(colData(SEobj))
        }
    }

    #Fix categories of metadata within Experiment object
    for (colm in 1:ncol(ptb)){
        numtest <- length(which(is.na(as.numeric(ptb[, colm])) == TRUE))
        if (numtest == 0){
            ptb[, colm] <- as.numeric(ptb[, colm])
        }
    }

    ExpObj <- SummarizedExperiment(assays = assays(SEobj), rowData = rowData(SEobj), colData = ptb)
    metadata(ExpObj) <- allmetadata

    return(ExpObj)
}


#' declare_filtering_presets(analysis = NULL, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL, maxl2fc = NULL, minl2fc = NULL)
#'
#' Performs vetting of a SummarizedExperiment object for use in several functions
#' @export
declare_filtering_presets <- function(analysis = NULL, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL, maxl2fc = NULL, minl2fc = NULL, minabscorrcoeff = NULL){

    if ((analysis != "LKT") && (!(is.null(GenomeCompletenessCutoff)))){
        warning("Genome completeness only makes sense for taxa. Please choose a taxonomic (non functional) analysis.")
        GenomeCompletenessCutoff <- NULL
    }

    presetlist <- list()

    if (!is.null(applyfilters)){
        if (applyfilters == "stringent"){
            if (analysis == "LKT"){
                presetlist$featcutoff <- c(2000, 15)
                presetlist$GenomeCompletenessCutoff <- c(30, 10)
                presetlist$PctFromCtgscutoff <- c(50, 10)
                presetlist$minl2fc <- 2
            } else {
                presetlist$featcutoff <- c(50, 15)
                presetlist$minl2fc <- 2.5
            }
            presetlist$minabscorrcoeff <- 0.8
        } else if (applyfilters == "moderate"){
            if (analysis == "LKT"){
                presetlist$featcutoff <- c(250, 15)
                presetlist$GenomeCompletenessCutoff <- c(10, 5)
                presetlist$PctFromCtgscutoff <- c(25, 10)
                presetlist$minl2fc <- 1
            } else {
                presetlist$featcutoff <- c(10, 5)
                presetlist$minl2fc <- 1
            }
            presetlist$minabscorrcoeff <- 0.6
        } else if (applyfilters == "light"){
            if (analysis == "LKT"){
                presetlist$featcutoff <- c(50, 5)
                presetlist$GenomeCompletenessCutoff <- c(5, 5)
                presetlist$PctFromCtgscutoff <- c(25, 10)
                presetlist$minl2fc <- 1
            } else {
                presetlist$featcutoff <- c(5, 5)
                presetlist$minl2fc <- 1
            }
            presetlist$minabscorrcoeff <- 0.4
        }
    }

    #Replace with any values explicitly set by the user
    argstoset <- c("featcutoff", "GenomeCompletenessCutoff", "PctFromCtgscutoff", "maxl2fc", "minl2fc", "minabscorrcoeff")[!unlist(lapply(list(featcutoff, GenomeCompletenessCutoff, PctFromCtgscutoff, maxl2fc, minl2fc, minabscorrcoeff), is.null))]

    if (length(argstoset) > 0){
        for (ats in argstoset){
            presetlist[[ats]] <- get(ats)
        }
    }

    #Generate a filtration message
    presetlist$filtermsg <- NULL
    #Discard features which do not match certain criteria
    if (!(is.null(presetlist$featcutoff))){
        presetlist$thresholdPPM <- presetlist$featcutoff[1]
        presetlist$sampcutoffpct <- min(presetlist$featcutoff[2], 100)
        presetlist$filtermsg <- paste("Feature must be >", presetlist$thresholdPPM, "PPM in at least ", presetlist$sampcutoffpct, "% of samples", sep = "")
    } else {
        presetlist$filtermsg <- NULL
        presetlist$featcutoff <- c(0, 0)
    }

    if (!(is.null(presetlist$PctFromCtgscutoff))){
        presetlist$thresholdPctFromCtgs <- presetlist$PctFromCtgscutoff[1]
        presetlist$sampcutoffpctPctFromCtgs <- min(presetlist$PctFromCtgscutoff[2], 100)
        presetlist$filtermsg <- paste(presetlist$filtermsg, (paste("Taxonomy information must come from >", presetlist$thresholdPctFromCtgs, "% contigs in at least ", presetlist$sampcutoffpctPctFromCtgs, "% of samples", sep = "")), sep = "\n")
    }

    if (!(is.null(presetlist$GenomeCompletenessCutoff))){
        presetlist$thresholdGenomeCompleteness <- presetlist$GenomeCompletenessCutoff[1]
        presetlist$sampcutoffpctGenomeCompleteness <- min(presetlist$GenomeCompletenessCutoff[2], 100)
        presetlist$filtermsg <- paste(presetlist$filtermsg, (paste("Taxon genome completeness must be >", presetlist$thresholdGenomeCompleteness, "% in at least ", presetlist$sampcutoffpctGenomeCompleteness, "% of samples", sep = "")), sep = "\n")
    }

    return(presetlist)
}

#' filter_correlations(corrmat = NULL, mincorrelcoeff = NULL)
#' Given a pairwise correlation matrix, eliminate features which do not present an absolute correlation coefficient smaller than mincorrelcoeff with all other features other than itself.
#'
#' @export
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
