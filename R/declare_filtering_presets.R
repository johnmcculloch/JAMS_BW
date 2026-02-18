#' declare_filtering_presets(analysis = NULL, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, maxl2fc = NULL, minl2fc = NULL)
#'
#' Performs vetting of a SummarizedExperiment object for use in several functions
#' @export

declare_filtering_presets <- function(analysis = NULL, is16S = FALSE, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, maxl2fc = NULL, minl2fc = NULL, minabscorrcoeff = NULL){

    taxonomic_spaces <- c("LKT", "Contig_LKT", "ConsolidatedGenomeBin", "MB2bin", "16S")

    if ((!analysis %in% taxonomic_spaces) && (!(is.null(GenomeCompletenessCutoff)))){
        warning("Genome completeness only makes sense for taxa. Please choose a taxonomic (non functional) analysis.")
        GenomeCompletenessCutoff <- NULL
    }

    presetlist <- list()

    if (!is.null(applyfilters)){
        if (applyfilters == "stringent"){
            if (analysis %in% taxonomic_spaces){
                presetlist$featcutoff <- c(2000, 15)
                presetlist$GenomeCompletenessCutoff <- c(30, 10)
                presetlist$minl2fc <- 2
            } else {
                presetlist$featcutoff <- c(50, 15)
                presetlist$minl2fc <- 2.5
            }
            presetlist$minabscorrcoeff <- 0.8
        } else if (applyfilters == "moderate"){
            if (analysis %in% taxonomic_spaces){
                presetlist$featcutoff <- c(250, 15)
                presetlist$GenomeCompletenessCutoff <- c(10, 5)
                presetlist$minl2fc <- 1
            } else {
                presetlist$featcutoff <- c(5, 5)
                presetlist$minl2fc <- 1
            }
            presetlist$minabscorrcoeff <- 0.6
        } else if (applyfilters == "light"){
            if (analysis %in% taxonomic_spaces){
                presetlist$featcutoff <- c(50, 5)
                presetlist$GenomeCompletenessCutoff <- c(5, 5)
                presetlist$minl2fc <- 1
            } else {
                presetlist$featcutoff <- c(0, 0)
                presetlist$minl2fc <- 1
            }
            presetlist$minabscorrcoeff <- 0.4
        }
    }

    #Replace with any values explicitly set by the user
    argstoset <- c("featcutoff", "GenomeCompletenessCutoff", "maxl2fc", "minl2fc", "minabscorrcoeff")[!unlist(lapply(list(featcutoff, GenomeCompletenessCutoff, maxl2fc, minl2fc, minabscorrcoeff), is.null))]

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

    if (!(is.null(presetlist$GenomeCompletenessCutoff))){
        if (!is16S){
            presetlist$thresholdGenomeCompleteness <- presetlist$GenomeCompletenessCutoff[1]
            presetlist$sampcutoffpctGenomeCompleteness <- min(presetlist$GenomeCompletenessCutoff[2], 100)
            presetlist$filtermsg <- paste(presetlist$filtermsg, (paste("Taxon genome completeness must be >", presetlist$thresholdGenomeCompleteness, "% in at least ", presetlist$sampcutoffpctGenomeCompleteness, "% of samples", sep = "")), sep = "\n")
        }
    }

    return(presetlist)
}
