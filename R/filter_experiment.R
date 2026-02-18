#' filter_experiment(SEobj = NULL, samplesToKeep = NULL, featuresToKeep = NULL, featmaxatleastPPM = 0, featcutoff = c(0, 0), GenomeCompletenessCutoff = NULL, GenomeContaminationCutoff = NULL, applyfilters = NULL, discard_SDoverMean_below = NULL, normalization = NULL, PPM_normalize_to_bases_sequenced = TRUE, flush_out_empty_samples = FALSE, clr_pseudocount = 1, give_info = TRUE)
#'
#' Filters a SummarizedExperiment object by several criteria.
#' @export

filter_experiment <- function(SEobj = NULL, samplesToKeep = NULL, featuresToKeep = NULL, featmaxatleastPPM = 0, featcutoff = c(0, 0), GenomeCompletenessCutoff = NULL, GenomeContaminationCutoff = NULL, applyfilters = NULL, discard_SDoverMean_below = NULL, normalization = NULL, PPM_normalize_to_bases_sequenced = TRUE, flush_out_empty_samples = FALSE, clr_pseudocount = 1, give_info = TRUE){

    #Check that, if using normalization, the method is a valid one.
    if (!is.null(normalization)){
        #Consider synonym
        if (normalization[1] == FALSE){
            normalization <- NULL
        }
        #Check validity of method
        if (!any(normalization %in% c("relabund", "clr"))){
          flog.warm('Normalization must be set to relative abundance with "relabund", or to CLR transform with "clr"')
          stop("Invalid parameter for normalization. Set to NULL for raw counts.")
        } else {
            normalization <- normalization[normalization %in% c("relabund", "clr")]
        }
    }

    #Get appropriate object with which to work
    if (as.character(class(SEobj)[1]) == "SummarizedExperiment"){
        #Get only samples you asked for
        if (!(is.null(samplesToKeep))){
            samplesToKeep <- colnames(SEobj)[colnames(SEobj) %in% samplesToKeep]
            #Stop if there is nothing left. You will be amazed as to what people ask for.
            if (length(samplesToKeep) < 2){
                flog.warn("There are less than 2 samples in samplesToKeep matching original samples in the input SEobj.")
                stop("Impossible to make a SummarizedExperiment object with less than 2 samples. Check samplesToKeep.")
            }
            SEobj <- SEobj[, samplesToKeep]

            #Also prune total counts vector in SEobj metadata
            for (mdTU in c("TotalBasesSequenced", "TotalBasesSequencedinAnalysis")){
                metadata(SEobj)[[mdTU]] <- metadata(SEobj)[[mdTU]][ , colnames(SEobj), drop = FALSE]
            }
        }

        #Get only features you asked for
        if (!(is.null(featuresToKeep))){
            featuresToKeep <- rownames(SEobj)[rownames(SEobj) %in% featuresToKeep]
            #Stop if there is nothing left. You will be amazed as to what people ask for.
            if (length(featuresToKeep) < 2){
                flog.warn("There are less than 2 features in featuresToKeep matching original features in the input SEobj.")
                stop("Impossible to make a SummarizedExperiment object with less than 2 features. Check featuresToKeep.")
            }
            SEobj <- SEobj[featuresToKeep, ]
        }

        if (PPM_normalize_to_bases_sequenced) {
            totbases <- metadata(SEobj)$TotalBasesSequenced
        } else {
            totbases <- metadata(SEobj)$TotalBasesSequencedinAnalysis
        }

        JAMSversion <- metadata(SEobj)$version

        if (is.null(JAMSversion)){
            JAMS2_SEobj <- FALSE
        } else {
            JAMS2_SEobj <- TRUE
        }

    } else {
        stop("SEobj must be a SummarizedExperiment object.")
    }

    #Declare useful function
    prop_above_threshold <- function(nummatix = NULL, feat = NULL, numthreshold = NULL) {
        proportionPassing <- length(which(nummatix[feat, ] >= numthreshold)) / ncol(nummatix)

        return(proportionPassing)
    }

    #Declare filtering presets
    presetlist <- declare_filtering_presets(analysis = metadata(SEobj)$analysis, applyfilters = applyfilters, featcutoff = featcutoff, GenomeCompletenessCutoff = GenomeCompletenessCutoff)

    #Flush out empty Samples
    if (flush_out_empty_samples){
        emptysamples <- names(which(colSums(rawcts) == 0) == TRUE)
        if (length(emptysamples) > 0){
            flog.info(paste("Samples", paste0(emptysamples, collapse = ", "), "are empty and will be discarded."))
            validsamples <- names(which(colSums(rawcts) > 0) == TRUE)
            SEobj <- SEobj[ , validsamples]
        }
    }

    countmat_raw <- as.matrix(assays(SEobj)$BaseCounts)
    countmat_raw <- countmat_raw[!is.na(row.names(countmat_raw)),]

    #Obtain a relabund matrix if relabund in normalization, or featcutoff is not equal to c(0,0)
    if (any(c(("relabund" %in% normalization), (all(presetlist$featcutoff == c(0,0))), !is.null(discard_SDoverMean_below)))){
        #transform into PPM
        getPPM <- function(countmat_raw = NULL, Sample = NULL){
            PPMs <- round(((countmat_raw[ , Sample] / totbases["NumBases", Sample]) * 1000000), 0)

            return(PPMs)
        }
        PPMmatrix <- sapply(1:ncol(countmat_raw), function(x){ getPPM(countmat_raw = countmat_raw, Sample = colnames(countmat_raw)[x])} )
        colnames(PPMmatrix) <- colnames(countmat_raw)

        #Commit to SEobj
        assays(SEobj)$PPM <- PPMmatrix
    }

    #Obtain a CLR transform matrix if requested.
    if ("clr" %in% normalization) {
        if (!is.numeric(clr_pseudocount) || length(clr_pseudocount) != 1 || clr_pseudocount <= 0){
            stop("clr_pseudocount must be a single positive numeric value (e.g., 1 or 0.5).")
        }

        getCLR <- function(countmat_raw = NULL, Sample = NULL, clr_pseudocount = NULL){
            x <- (as.numeric(countmat_raw[, Sample]) + clr_pseudocount)
            CLRs <- as.numeric(compositions::clr(x))

            return(CLRs)
        }
        CLRmatrix <- sapply(1:ncol(countmat_raw), function(x){ getCLR(countmat_raw = countmat_raw, Sample = colnames(countmat_raw)[x], clr_pseudocount = clr_pseudocount)} )

        colnames(CLRmatrix) <- colnames(countmat_raw)
        rownames(CLRmatrix) <- rownames(countmat_raw)

        #Commit to SEobj
        assays(SEobj)$CLR <- CLRmatrix
    }

    #Discard features which do not match certain criteria
    if (!(all(presetlist$featcutoff == c(0,0)))){
        if (length(presetlist$featcutoff) != 2){
            stop("Please specify the minimum PPM in what percentage of the samples you want with a numerical vector of size two. For example, featcutoff = c(2000,10) would discard features which are not at least 2000 PPM in at least 10% of samples.")
        }
        thresholdPPM <- presetlist$featcutoff[1]
        sampcutoffpct <- min(presetlist$featcutoff[2], 100)
        #New rule. If sample percent threshold is 0, consider a single sample. i.e. at least thresholdPPM in a single sample.
        if (sampcutoffpct == 0){
            #Calculate what percentage is exactly equal to a single sample.
            flog.warn("Percentage of samples for filtration (prevalence) was set to 0%. This would lead to no filtering, so instead will remove features whose PPM does not cross the threshold in at least one sample. To turn off all PPM filtration, instead, set featcutoff to c(0,0).")
            sampcutoffpct <- (1 / ncol(PPMmatrix)) * 100
        }

        featuresToKeep2 <- rownames(countmat_raw)[which((sapply(1:nrow(countmat_raw), function(x) { prop_above_threshold(nummatix = PPMmatrix, feat = x, numthreshold = thresholdPPM) })) > (sampcutoffpct / 100))]

        SEobj <- SEobj[featuresToKeep2, ]
        countmat_raw <- countmat_raw[featuresToKeep2, ]

        cutoffmsg <- paste("Feature must be >", thresholdPPM, "PPM in at least", sampcutoffpct, "% of samples", sep = " ")
        if (give_info){
            flog.info(cutoffmsg)
        }
    }

    if (!is.null(discard_SDoverMean_below)){
        dfm <- (rowSds(PPMmatrix) / rowMeans(PPMmatrix))
        featuresToKeepSDflt <- names(dfm[dfm > discard_SDoverMean_below])
        SEobj <- SEobj[featuresToKeepSDflt, ]
        countmat_raw <- countmat_raw[featuresToKeepSDflt, ]
        SDoverMeanmsg <- paste("Feature must have >", discard_SDoverMean_below, "SDs over mean", sep=" ")
        flog.info(SDoverMeanmsg)
    }

    #Filter by Genome Completeness
    if (all(c(("GenomeCompleteness" %in% names(assays(SEobj))), (!(all(presetlist$GenomeCompletenessCutoff == c(0,0))))))){
        if (length(presetlist$GenomeCompletenessCutoff) != 2){
            stop("Please specify the minimum Genome Completeness in what percentage of the samples you want with a numerical vector of size two. For example, GenomeCompletenessCutoff = c(10, 5) would discard features whose genome completeness is smaller than 10% from contigs in at least 5% of samples.")
        }
        if (JAMS2_SEobj) {
            minGenomeCompleteness <- presetlist$GenomeCompletenessCutoff[1]
        } else {
            minGenomeCompleteness <- (presetlist$GenomeCompletenessCutoff[1] / 100)
        }

        GenomeCompletenessmatrix <- assays(SEobj)[["GenomeCompleteness"]]

        sampcutoffpct <- min(presetlist$GenomeCompletenessCutoff[2], 100)
        #New rule. If sample percent threshold is 0, consider a single sample. i.e. at least thresholdPPM in a single sample.
        if (sampcutoffpct == 0){
            #Calculate what percentage is exactly equal to a single sample.
            flog.warn("Percentage of samples for filtration (prevalence) was set to 0%. This would lead to no filtering, so instead will remove features whose GenomeCompleteness does not cross the threshold in at least one sample. To turn off all GenomeCompleteness filtration, instead, set GenomeCompletenessCutoff to c(0,0).")
            sampcutoffpct <- (1 / ncol(GenomeCompletenessmatrix)) * 100
        }

        featuresToKeep4 <- rownames(GenomeCompletenessmatrix)[which((sapply(1:nrow(GenomeCompletenessmatrix), function(x) { prop_above_threshold(nummatix = GenomeCompletenessmatrix, feat = x, numthreshold = minGenomeCompleteness) })) >= (sampcutoffpct / 100))]

        cutoffmsg <- paste0("Feature must have a genome completeness >= ", minGenomeCompleteness, "% in at least ", sampcutoffpct, "% of samples", sep = " ")
        if (give_info){
            flog.info(cutoffmsg)
        }

        SEobj <- SEobj[featuresToKeep4, ]
        countmat_raw <- countmat_raw[featuresToKeep4, ]
    }

    #filter out unwanted feature/sample stratification if SummarizedExperiment object contains that kind of information
    #This is dealt with quite differently whether it is a JAMS1 or JAMS2 SEobj. Will keep functionality for JAMS1 for now for backwards compatibility, although this is now deprecated.

    if (JAMS2_SEobj) {
        if (all(c("allfeaturesbytaxa_matrix", "allfeaturesbytaxa_GeneCounts_matrix") %in% names(metadata(SEobj)))){
            #No need to eliminate unwanted samples from SparseIndex matrix, because this is done when subsetting the entire SEobj at the samplesToKeep phase above.
            #Obtain valid index row numbers to prune sparse matrices.
            index_matrix <- as.data.frame(assays(SEobj)[["SparseIndex"]])
            rowsinterestdf <- as.data.frame(index_matrix) %>% rownames_to_column("Accession") %>% pivot_longer(cols = -Accession, names_to = "Sample", values_to = "RowNumber")
            rowsinterestdf <- as.data.frame(rowsinterestdf)
            #Ensure correct class
            rowsinterestdf$RowNumber <- as.integer(rowsinterestdf$RowNumber)
            #Eliminate unavailable information.
            rowsinterestdf <- rowsinterestdf[!is.na(rowsinterestdf$RowNumber), ]

            #Process sparse martices
            allfeaturesbytaxa_matrix <- metadata(SEobj)[["allfeaturesbytaxa_matrix"]][rowsinterestdf$RowNumber, , drop = FALSE]
            #Prune empty taxa, i.e. taxa which have 0 counts for all of the features of interest
            LKTskeep <- Matrix::colSums(allfeaturesbytaxa_matrix) != 0
            allfeaturesbytaxa_matrix <- allfeaturesbytaxa_matrix[, LKTskeep, drop = FALSE]

            #Prune the genecounts matrix accordingly
            allfeaturesbytaxa_GeneCounts_matrix <- metadata(SEobj)[["allfeaturesbytaxa_GeneCounts_matrix"]][rowsinterestdf$RowNumber, , drop = FALSE]
            allfeaturesbytaxa_GeneCounts_matrix <- allfeaturesbytaxa_GeneCounts_matrix[, colnames(allfeaturesbytaxa_matrix), drop = FALSE]

            #Reset row numbers for SparseIndex and the sparse matrix
            rowsinterestdf$RowNumber <- 1:nrow(rowsinterestdf)
            rownames(allfeaturesbytaxa_matrix) <- 1:nrow(allfeaturesbytaxa_matrix)
            rownames(allfeaturesbytaxa_GeneCounts_matrix) <- 1:nrow(allfeaturesbytaxa_GeneCounts_matrix)

            #Commit to SEobj metadata
            metadata(SEobj)[["allfeaturesbytaxa_matrix"]] <- allfeaturesbytaxa_matrix
            metadata(SEobj)[["allfeaturesbytaxa_GeneCounts_matrix"]] <- allfeaturesbytaxa_GeneCounts_matrix

            #Regenerate SparseIndex assay and insert back into SummarizedExperiment object
            SparseIndex <- rowsinterestdf %>% group_by(Sample, Accession) %>% tidyr::pivot_wider(names_from = Sample, values_from = RowNumber, values_fill = NA)
            SparseIndex <- as.data.frame(SparseIndex)
            rownames(SparseIndex) <- SparseIndex$Accession
            SparseIndex$Accession <- NULL
            #Add empty features if necessary (some may have been pruned)
            emptyFeatures <- rownames(SEobj)[!(rownames(SEobj) %in% rownames(SparseIndex))]
            if (length(emptyFeatures) > 0){
                complementarycts <- matrix(ncol = ncol(SEobj), nrow = length(emptyFeatures), data = 0)
                complementarycts <- as.data.frame(complementarycts, stringsAsFactors = FALSE)
                colnames(complementarycts) <- colnames(SEobj)
                rownames(complementarycts) <- emptyFeatures
                SparseIndex <- rbind(SparseIndex, complementarycts)
            }
            SparseIndex <- SparseIndex[rownames(SEobj), colnames(SEobj)]
            assays(SEobj)[["SparseIndex"]] <- as.matrix(SparseIndex)
        }

    } else {

        #Treat it as a JAMS1 style stratification object
        if (all(c("allfeaturesbytaxa_matrix", "allfeaturesbytaxa_GeneCounts_matrix") %in% names(metadata(SEobj)))){
            #Obtain current index and matrix
            allfeaturesbytaxa_index <- metadata(SEobj)$allfeaturesbytaxa_index
            allfeaturesbytaxa_matrix <- metadata(SEobj)$allfeaturesbytaxa_matrix

            #Eliminate from index irrelevant samples
            allfeaturesbytaxa_index <- subset(allfeaturesbytaxa_index, Sample %in% colnames(SEobj))

            #Eliminate from index irrelevant features
            allfeaturesbytaxa_index <- subset(allfeaturesbytaxa_index, Accession %in% rownames(SEobj))

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

            #Stick back into SEobj
            metadata(SEobj)$allfeaturesbytaxa_index <- allfeaturesbytaxa_index
            metadata(SEobj)$allfeaturesbytaxa_matrix <- allfeaturesbytaxa_matrix
        }
    }

    return(SEobj)
}
