#' filter_experiment(SEobj = NULL, samplesToKeep = NULL, featuresToKeep = NULL, featmaxatleastPPM = 0, featcutoff = c(0, 0), GenomeCompletenessCutoff = NULL, GenomeContaminationCutoff = NULL, applyfilters = NULL, discard_SDoverMean_below = NULL, normalization = NULL, PPM_normalize_to_bases_sequenced = TRUE, only_allow_CSBs = FALSE, flush_out_empty_samples = FALSE, clr_pseudocount = 1, give_info = TRUE)
#'
#' Filters a SummarizedExperiment object by several criteria.
#' @export

filter_experiment <- function(SEobj = NULL, samplesToKeep = NULL, featuresToKeep = NULL, featmaxatleastPPM = 0, featcutoff = c(0, 0), GenomeCompletenessCutoff = NULL, GenomeContaminationCutoff = NULL, applyfilters = NULL, discard_SDoverMean_below = NULL, normalization = NULL, PPM_normalize_to_bases_sequenced = TRUE, only_allow_CSBs = FALSE, flush_out_empty_samples = FALSE, clr_pseudocount = 1, give_info = TRUE){

    #Check that, if using normalization, the method is a valid one.
    if (!is.null(normalization)){
        #Check validity of method
        if (!any(normalization %in% c("relabund", "clr"))){
          flog.warn('Normalization must be set to relative abundance with "relabund", or to CLR transform with "clr"')
          stop("Invalid parameter for normalization. Set to NULL for raw counts.")
        } else {
            normalization <- normalization[normalization %in% c("relabund", "clr")]
        }
    }

    #Get appropriate object with which to work
    if (as.character(class(SEobj)[1]) == "SummarizedExperiment"){
        #Reject deprecated JAMS1 objects up front. This used to be a silent branch that rescaled
        #the genome completeness cutoff by 100 for objects with no version stamp.
        assert_JAMS2_SEobj(SEobj = SEobj, objname = "SEobj", calling_function = "filter_experiment")

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

        #If requesting only CSBs, check that the request is valid and that there are any.
        if (only_allow_CSBs){
            #Check this is a ConsolidatedGenomeBin SEobj
            if (metadata(SEobj)$analysis != "ConsolidatedGenomeBin"){
                flog.warn("Restriction of taxonomic features to CSBs only makes sense for a ConsolidatedGenomeBin analysis type JAMS2 SummarizedExperiment object.")
                flog.warn(paste("Ignoring request to restrict features to CSBs, this is a", metadata(SEobj)$analysis, "type SEobj."))
            } else {
                CSBsToKeep <- rownames(SEobj)[grep("^CSB_", rownames(SEobj))]
                #Stop if there is nothing left. You will be amazed as to what people use as inputs.
                if (length(CSBsToKeep) < 2){
                    flog.warn("There are less than 2 features in which are valid CSBs in the input SEobj.")
                    stop("Impossible to make a SummarizedExperiment object with less than 2 features. Check featuresToKeep.")
                }
                PctCSBs <- paste0(round(((length(CSBsToKeep) / nrow(SEobj)) * 100), 1), "%")
                if (give_info){
                    flog.warn(paste("Filtering out features which are not CSBs (Consolidated Species Bins), i.e. features which were not binned by MetaBAT2 in JAMSalpha."))
                    flog.info(paste("CSBs kept are", PctCSBs, "of the number of features on the original SEobj."))
                }
                SEobj <- SEobj[CSBsToKeep, ]
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
        emptysamples <- names(which(colSums(assays(SEobj)$BaseCounts) == 0) == TRUE)
        if (length(emptysamples) > 0){
            flog.info(paste("Samples", paste0(emptysamples, collapse = ", "), "are empty and will be discarded."))
            validsamples <- names(which(colSums(assays(SEobj)$BaseCounts) > 0) == TRUE)
            SEobj <- SEobj[ , validsamples]
        }
    }

    countmat_raw <- as.matrix(assays(SEobj)$BaseCounts)
    countmat_raw <- countmat_raw[!is.na(row.names(countmat_raw)),]

    #Obtain a relabund matrix if relabund in normalization, or featcutoff is not equal to c(0,0)
    if (any(c(("relabund" %in% normalization), (all(presetlist$featcutoff != c(0,0))), !is.null(discard_SDoverMean_below)))){
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
            #Centred log-ratio, identical to as.numeric(compositions::clr(x)) for a strictly positive vector.
            CLRs <- as.numeric(log(x) - mean(log(x)))

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
        dfm <- (jams_rowSds(PPMmatrix) / rowMeans(PPMmatrix))
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
        #GenomeCompleteness is stored on a 0-100 percent scale (CheckM2 style) in every JAMS2
        #object, so the cutoff is used as supplied. JAMS1 objects stored a 0-1 fraction and were
        #silently rescaled here; they are now rejected by assert_JAMS2_SEobj above.
        minGenomeCompleteness <- presetlist$GenomeCompletenessCutoff[1]

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

    #NOTE ON FEATURE-BY-TAXON STRATIFICATION.
    #This function deliberately does NOT prune the stratification sparse matrices in
    #metadata(SEobj) (allfeaturesbytaxa_matrix and allfeaturesbytaxa_GeneCounts_matrix). It does
    #not need to: the row pointers live in the SparseIndex ASSAY, which SummarizedExperiment
    #subsets along with the object, so every lookup made by retrieve_features_by_taxa stays
    #correct after any filtration performed above. The matrices simply keep rows and taxon
    #columns that nothing points at any more, which costs memory but never correctness.
    #A large commented-out block that used to prune and renumber them (plus a parallel JAMS1
    #branch keyed on the long-form allfeaturesbytaxa_index) was removed in JAMS 2.3.1. If
    #pruning is ever reinstated for memory reasons, note that it must renumber SparseIndex and
    #both matrices consistently, and that merge_JAMS_SEobj does exactly that renumbering already.


    return(SEobj)
}
