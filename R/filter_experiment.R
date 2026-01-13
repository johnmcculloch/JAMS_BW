#' filter_experiment(ExpObj = NULL, featmaxatleastPPM = 0, featcutoff = c(0, 0), discard_SDoverMean_below = NULL, samplesToKeep = NULL, featuresToKeep = NULL, normalization = "relabund", asPPM = TRUE, PPM_normalize_to_bases_sequenced = FALSE, flush_out_empty_samples = FALSE, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL)
#'
#' Filters a SummarizedExperiment object by several criteria.
#' @export

filter_experiment <- function(ExpObj = NULL, featmaxatleastPPM = 0, featcutoff = c(0, 0), discard_SDoverMean_below = NULL, samplesToKeep = NULL, featuresToKeep = NULL, normalization = "relabund", asPPM = TRUE, PPM_normalize_to_bases_sequenced = FALSE, flush_out_empty_samples = FALSE, GenomeCompletenessCutoff = NULL, PctFromCtgscutoff = NULL){

    #Get only samples you asked for
    if (!(is.null(samplesToKeep))){
        samplesToKeep <- samplesToKeep[samplesToKeep %in% colnames(ExpObj)]
        ExpObj <- ExpObj[, samplesToKeep]
        #Also fix total counts vector in ExpObj metadata
        for (mdTU in c("TotalBasesSequenced", "TotalBasesSequencedinAnalysis")){
            MTU <- t(as.matrix(metadata(ExpObj)[[mdTU]]["NumBases", samplesToKeep]))
            rownames(MTU) <- "NumBases"
            colnames(MTU) <- samplesToKeep
            metadata(ExpObj)[[mdTU]] <- MTU
        }
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
        JAMSversion <- assays(ExpObj)$version
        if (is.null(JAMSversion)){
            JAMS2_SEobj <- TRUE
        } else {
            JAMS2_SEobj <- FALSE
        }
    } else {
        stop("Object must be a SummarizedExperiment object.")
    }

    #If setting featmaxatleastPPM or featcutoff to anything other than the defaults, then return as PPM.
    if (any(featcutoff != c(0, 0))){
        asPPM <- TRUE
    }

    #Flush out empty rows
    #Don't do this anymore
    #ExpObj <- ExpObj[(rowSums(rawcts) > 0), ]

    #Flush out empty Samples
    if (flush_out_empty_samples){
        emptysamples <- names(which(colSums(rawcts) == 0) == TRUE)
        if (length(emptysamples) > 0){
            flog.info(paste("Samples", paste0(emptysamples, collapse = ", "), "are empty and will be discarded."))
            validsamples <- names(which(colSums(rawcts) > 0) == TRUE)
            ExpObj <- ExpObj[ , validsamples]
        }
    }

    countmat <- as.matrix(assays(ExpObj)$BaseCounts)
    countmat <- countmat[!is.na(row.names(countmat)),]

    if (normalization == "compositions"){
        asPPM <- FALSE
    }

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

        if (normalization == "compositions"){
            countmat2 <- sapply(1:ncol(countmat), function(x){ as.numeric(compositions::clr(countmat[,x])) })
            colnames(countmat2) <- colnames(countmat)
            rownames(countmat2) <- rownames(countmat)
            countmatrix <- countmat2
        } else {
            countmatrix <- countmat
        }
    }

    #Discard features which do not match certain criteria
    if (!(missing(featcutoff))){
        if (length(featcutoff) != 2){
            stop("Please specify the minimum PPM in what percentage of the samples you want with a numerical vector of size two. For example, featcutoff=c(2000,10) would discard features which are not at least 2000 PPM in at least 10% of samples.")
        }
        thresholdPPM <- featcutoff[1]
        sampcutoffpct <- min(featcutoff[2], 100)

        prop_above_threshold <- function(countmatrix = NULL, feat = NULL, thresholdPPM = NULL) {
            proportionPassing <- length(which(countmatrix[feat, ] >= thresholdPPM)) / ncol(countmatrix)
            return(proportionPassing)
        }

        featuresToKeep2 <- rownames(countmatrix)[which((sapply(1:nrow(countmatrix), function(x) { prop_above_threshold(countmatrix = countmatrix, feat = x, thresholdPPM = thresholdPPM) })) > (sampcutoffpct/100))]

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

    #Filter by Genome Completeness
    if (all(c(("GenomeCompleteness" %in% names(assays(ExpObj))), (!is.null(GenomeCompletenessCutoff))))){
        if(length(GenomeCompletenessCutoff) != 2){
            stop("Please specify the minimum Genome Completeness in what percentage of the samples you want with a numerical vector of size two. For example, GenomeCompletenessCutoff = c(10, 5) would discard features whose genome completeness is smaller than 10% from contigs in at least 5% of samples.")
        }
        if (JAMS2_SEobj) {
            minGenomeCompleteness <- GenomeCompletenessCutoff[1]
        } else {
            minGenomeCompleteness <- (GenomeCompletenessCutoff[1] / 100)
        }

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

    #filter out unwanted feature/sample stratification if SummarizedExperiment object contains that kind of information
    #This is dealt with quite differently whether it is a JAMS1 or JAMS2 SEobj. Will keep functionality for JAMS1 for now for backwards compatibility, although this is now deprecated.

    if (JAMS2_SEobj) {
        if (all(c("allfeaturesbytaxa_matrix", "allfeaturesbytaxa_GeneCounts_matrix") %in% names(metadata(ExpObj)))){
            #No need to eliminate unwanted samples from SparseIndex matrix, because this is done when subsetting the entire ExpObj at the samplesToKeep phase above.
            #Obtain valid index row numbers to prune sparse matrices.
            index_matrix <- as.data.frame(assays(ExpObj)[["SparseIndex"]])
            rowsinterestdf <- as.data.frame(index_matrix) %>% rownames_to_column("Accession") %>% pivot_longer(cols = -Accession, names_to = "Sample", values_to = "RowNumber")
            rowsinterestdf <- as.data.frame(rowsinterestdf)
            #Ensure correct class
            rowsinterestdf$RowNumber <- as.integer(rowsinterestdf$RowNumber)
            #Eliminate unavailable information.
            rowsinterestdf <- rowsinterestdf[!is.na(rowsinterestdf$RowNumber), ]

            #Process sparse martices
            allfeaturesbytaxa_matrix <- metadata(ExpObj)[["allfeaturesbytaxa_matrix"]][rowsinterestdf$RowNumber, , drop = FALSE]
            #Prune empty taxa, i.e. taxa which have 0 counts for all of the features of interest
            LKTskeep <- Matrix::colSums(allfeaturesbytaxa_matrix) != 0
            allfeaturesbytaxa_matrix <- allfeaturesbytaxa_matrix[, LKTskeep, drop = FALSE]

            #Prune the genecounts matrix accordingly
            allfeaturesbytaxa_GeneCounts_matrix <- metadata(ExpObj)[["allfeaturesbytaxa_GeneCounts_matrix"]][rowsinterestdf$RowNumber, , drop = FALSE]
            allfeaturesbytaxa_GeneCounts_matrix <- allfeaturesbytaxa_GeneCounts_matrix[, colnames(allfeaturesbytaxa_matrix), drop = FALSE]

            #Reset row numbers for SparseIndex and the sparse matrix
            rowsinterestdf$RowNumber <- 1:nrow(rowsinterestdf)
            rownames(allfeaturesbytaxa_matrix) <- 1:nrow(allfeaturesbytaxa_matrix)
            rownames(allfeaturesbytaxa_GeneCounts_matrix) <- 1:nrow(allfeaturesbytaxa_GeneCounts_matrix)

            #Commit to SEobj metadata
            metadata(ExpObj)[["allfeaturesbytaxa_matrix"]] <- allfeaturesbytaxa_matrix
            metadata(ExpObj)[["allfeaturesbytaxa_GeneCounts_matrix"]] <- allfeaturesbytaxa_GeneCounts_matrix

            #Regenerate SparseIndex assay and insert back into SummarizedExperiment object
            SparseIndex <- rowsinterestdf %>% group_by(Sample, Accession) %>% tidyr::pivot_wider(names_from = Sample, values_from = RowNumber, values_fill = NA)
            SparseIndex <- as.data.frame(SparseIndex)
            rownames(SparseIndex) <- SparseIndex$Accession
            SparseIndex$Accession <- NULL
            #Add empty features if necessary (some may have been pruned)
            emptyFeatures <- rownames(ExpObj)[!(rownames(ExpObj) %in% rownames(SparseIndex))]
            if (length(emptyFeatures) > 0){
                complementarycts <- matrix(ncol = ncol(ExpObj), nrow = length(emptyFeatures), data = 0)
                complementarycts <- as.data.frame(complementarycts, stringsAsFactors = FALSE)
                colnames(complementarycts) <- colnames(ExpObj)
                rownames(complementarycts) <- emptyFeatures
                SparseIndex <- rbind(SparseIndex, complementarycts)
            }
            SparseIndex <- SparseIndex[rownames(ExpObj), colnames(ExpObj)]
            assays(ExpObj)[["SparseIndex"]] <- as.matrix(SparseIndex)
        }

    } else {

        #Treat it as a JAMS1 style stratification object
        if (all(c("allfeaturesbytaxa_matrix", "allfeaturesbytaxa_GeneCounts_matrix") %in% names(metadata(ExpObj)))){
            #Obtain current index and matrix
            allfeaturesbytaxa_index <- metadata(ExpObj)$allfeaturesbytaxa_index
            allfeaturesbytaxa_matrix <- metadata(ExpObj)$allfeaturesbytaxa_matrix

            #Eliminate from index irrelevant samples
            allfeaturesbytaxa_index <- subset(allfeaturesbytaxa_index, Sample %in% colnames(ExpObj))

            #Eliminate from index irrelevant features
            allfeaturesbytaxa_index <- subset(allfeaturesbytaxa_index, Accession %in% rownames(ExpObj))

            #Subset sparse matrix to contain only relevant rows
            allfeaturesbytaxa_matrix <- allfeaturesbytaxa_matrix[allfeaturesbytaxa_index$RowNumber, ]

            #Rename rows with consecutive numbers
            allfeaturesbytaxa_index$RowNumber <- 1:nrow(allfeaturesbytaxa_index)
            rownames(allfeaturesbytaxa_index) <- allfeaturesbytaxa_index$RowNumber
            rownames(allfeaturesbytaxa_matrix) <- allfeaturesbytaxa_index$RowNumber

            #Eliminate irrelevant columns from sparse matrix
            matcSums <- Matrix::colSums(allfeaturesbytaxa_matrix)
            matcSums <- matcSums[matcSums > 0]
            taxaToKeep <- names(matcSums)
            allfeaturesbytaxa_matrix <- allfeaturesbytaxa_matrix[ , taxaToKeep]

            #Stick back into ExpObj
            metadata(ExpObj)$allfeaturesbytaxa_index <- allfeaturesbytaxa_index
            metadata(ExpObj)$allfeaturesbytaxa_matrix <- allfeaturesbytaxa_matrix
        }
    }

    return(ExpObj)
}
