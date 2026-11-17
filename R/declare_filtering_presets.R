#' declare_filtering_presets(analysis = NULL, is16S = FALSE, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, ...)
#'
#' Resolves the JAMS filtration presets into the concrete cutoffs used by several functions.
#'
#' A filtration preset decides one thing only: whether a feature survives into the analysis at
#' all. It therefore sets featcutoff and, for taxonomic spaces, GenomeCompletenessCutoff, and
#' nothing else. Effect-size thresholds such as minl2fc, maxl2fc and minabscorrcoeff are NOT
#' filtration criteria. They decide how a surviving feature is judged or drawn, which is an
#' entirely separate question, and they are consequently set by their own arguments in the
#' functions which offer them, each with its own sensible default. Prior to JAMS 2.3.2 the
#' presets did set minl2fc and minabscorrcoeff, which coupled "is this feature real" to "is
#' this feature interesting" and made a plain applyfilters string silently change the
#' statistical verdict of a plot. That coupling has been removed.
#'
#' @export

declare_filtering_presets <- function(analysis = NULL, is16S = FALSE, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, ...){

    taxonomic_spaces <- c("LKT", "Contig_LKT", "ConsolidatedGenomeBin", "MB2bin", "16S")
    #Resolved once. A bare `analysis %in% taxonomic_spaces` is logical(0) when analysis is NULL,
    #and `if (logical(0))` is a hard error, so every branch below would abort on a NULL analysis.
    is_taxonomic <- ((!is.null(analysis)) && (length(analysis) == 1) && (analysis %in% taxonomic_spaces))

    #Named, rather than silently swallowed, so that a stale caller is visible in the log
    #instead of quietly having its effect-size threshold ignored.
    deadargs <- names(list(...))
    if (is.null(deadargs)){
        deadargs <- character(0)
    }
    deadargs <- deadargs[nzchar(deadargs)]
    if (length(deadargs) > 0){
        effectargs <- intersect(deadargs, c("minl2fc", "maxl2fc", "minabscorrcoeff"))
        if (length(effectargs) > 0){
            flog.warn(paste0("declare_filtering_presets no longer handles effect-size thresholds, so the following argument(s) were ignored: ", paste0(effectargs, collapse = ", "), ". Filtration presets now govern feature survival only. Set these on the plotting function itself, where they have their own defaults (minl2fc = 1, minabscorrcoeff = 0.4)."))
        }
        otherargs <- setdiff(deadargs, effectargs)
        if (length(otherargs) > 0){
            flog.warn(paste0("declare_filtering_presets does not recognize the following argument(s), which were ignored: ", paste0(otherargs, collapse = ", "), "."))
        }
    }

    #analysis must be tested for NULL first: `NULL %in% x` is logical(0), and `!logical(0) && TRUE`
    #aborts with "invalid 'x' type in 'x && y'" on newer R rather than falling through.
    if ((!is_taxonomic) && (!(is.null(GenomeCompletenessCutoff)))){
        flog.warn("Genome completeness only makes sense for taxa. Please choose a taxonomic (non functional) analysis.")
        GenomeCompletenessCutoff <- NULL
    }

    presetlist <- list()

    if (!is.null(applyfilters)){
        if (applyfilters == "stringent"){
            if (is_taxonomic){
                presetlist$featcutoff <- c(2000, 15)
                presetlist$GenomeCompletenessCutoff <- c(30, 10)
            } else {
                presetlist$featcutoff <- c(50, 15)
            }
        } else if (applyfilters == "moderate"){
            if (is_taxonomic){
                presetlist$featcutoff <- c(250, 15)
                presetlist$GenomeCompletenessCutoff <- c(10, 5)
            } else {
                presetlist$featcutoff <- c(5, 5)
            }
        } else if (applyfilters == "light"){
            if (is_taxonomic){
                presetlist$featcutoff <- c(50, 5)
                presetlist$GenomeCompletenessCutoff <- c(5, 5)
            } else {
                presetlist$featcutoff <- c(0, 0)
            }
        }
    }

    #Replace with any values explicitly set by the user
    argstoset <- c("featcutoff", "GenomeCompletenessCutoff")[!unlist(lapply(list(featcutoff, GenomeCompletenessCutoff), is.null))]

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
