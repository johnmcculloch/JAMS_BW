#' filter_experiment(ExpObj = NULL, featmaxatleastPPM = 0, featcutoff = c(0, 0), samplesToKeep = NULL, featuresToKeep = NULL, asPPM = TRUE, PPM_normalize_to_bases_sequenced = FALSE, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL, mgSeqnorm = FALSE)
#'
#' Filters a SummarizedExperiment or metagenomeSeq object by several criteria.
#' @export

filter_experiment <- function(ExpObj = NULL, featmaxatleastPPM = 0, featcutoff = c(0, 0), samplesToKeep = NULL, featuresToKeep = NULL, asPPM = TRUE, PPM_normalize_to_bases_sequenced = FALSE, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL, mgSeqnorm = FALSE){

    #Get only samples you asked for
    if (!(is.null(samplesToKeep))){
        samplesToKeep <- samplesToKeep[samplesToKeep %in% colnames(ExpObj)]
        ExpObj <- ExpObj[, samplesToKeep]
    }

    #Get appropriate object with which to work
    if (as.character(class(ExpObj)) == "SummarizedExperiment"){
        rawcts <- as.matrix(assays(ExpObj)$BaseCounts)
        if (PPM_normalize_to_bases_sequenced) {
            totbases <- metadata(ExpObj)$TotalBasesSequenced
        } else {
            #totbases <- metadata(ExpObj)$TotalBasesSequencedinAnalysis
            totbases <- t(data.frame(NumBases = colSums(rawcts)))
        }

    } else if (as.character(class(ExpObj)) == "MRexperiment"){
        rawcts <- MRcounts(ExpObj, norm = mgSeqnorm, log = FALSE)
        totbases <- colSums(rawcts)

    } else {
        stop("Input must be either a SummarizedExperiment or a MetagenomeSeq object.")
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

    if (as.character(class(ExpObj)) == "SummarizedExperiment"){
        countmat <- as.matrix(assays(ExpObj)$BaseCounts)
    } else if (as.character(class(ExpObj)) == "MRexperiment"){
        countmat <- MRcounts(ExpObj, norm = mgSeqnorm, log = FALSE)
    }
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
        if (as.character(class(ExpObj)) == "SummarizedExperiment"){
            assays(ExpObj)$BaseCounts <- countmatrix

        } else if (as.character(class(ExpObj)) == "MRexperiment"){
            pheno2 <- pData(ExpObj)
            #just make sure
            pheno2 <- pheno2[colnames(countmatrix), ]
            ftt <- fData(ExpObj)
            ##Create a non-normalised class object in metagenomeseq (mgseq)
            phenotypeData <- AnnotatedDataFrame(pheno2)
            ttdata <- AnnotatedDataFrame(ftt)
            ExpObj <- newMRexperiment(countmatrix, phenoData = phenotypeData, featureData = ttdata)
        }
    } else {
        countmatrix <- countmat
    }

    #Discard features which do not match certain criteria
    if(!(missing(featcutoff))){
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

    #Allow for filtering by genomecompleteness or percent from contigs, if SummarizedExperiment
    if (as.character(class(ExpObj)) == "SummarizedExperiment"){

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

    }#End conditional that object is a SummarizedExperiment

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
