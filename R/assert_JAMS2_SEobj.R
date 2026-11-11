#' assert_JAMS2_SEobj(SEobj = NULL, objname = "SEobj", require_stratification = FALSE, calling_function = NULL)
#'
#' Validates that an object is a JAMS2-generation JAMS-style SummarizedExperiment object, and
#' stops with an informative error if it is not.
#'
#' WHY THIS EXISTS. Up to JAMS 2.3.0 several downstream functions carried silent fallback paths
#' for JAMS1-era objects: filter_experiment rescaled GenomeCompletenessCutoff by a factor of 100,
#' and retrieve_features_by_taxa / constrain_JAMS_SEobj_to_wanted_taxa read a long-form
#' "allfeaturesbytaxa_index" metadata slot that JAMS2 objects do not have. Those branches were a
#' standing hazard: merge_JAMS_SEobj used to drop metadata$version, which silently reclassified a
#' perfectly good JAMS2 object as JAMS1 and then quietly divided the genome completeness cutoff
#' by 100 - a wrong answer rather than a crash. JAMS1 objects are long deprecated and nobody
#' should be analysing with them, so as of JAMS 2.3.1 the fallbacks are gone and every entry
#' point asserts the generation up front instead.
#'
#' The test is the presence of metadata(SEobj)$version, which is stamped by
#' make_SummarizedExperiments() and make_SummarizedExperiment_from_tables(). Where the string is
#' parseable, a major version below 2 is also rejected.
#'
#' @param SEobj The object to validate.
#' @param objname Name of the argument in the calling function, used to make the error message
#'   point at the right thing. Default "SEobj".
#' @param require_stratification If TRUE, additionally require that the object carries
#'   feature-by-taxon stratification (the JAMSbeta -u payload): a SparseIndex assay and an
#'   allfeaturesbytaxa_matrix metadata slot. Default FALSE.
#' @param calling_function Optional name of the calling function, prepended to the error.
#' @return Invisibly, the version string found on the object.
#' @export

assert_JAMS2_SEobj <- function(SEobj = NULL, objname = "SEobj", require_stratification = FALSE, calling_function = NULL){

    prefix <- if (!is.null(calling_function)) paste0(calling_function, ": ") else ""

    if (!inherits(SEobj, "SummarizedExperiment")){
        stop(paste0(prefix, objname, " must be a SummarizedExperiment object, but an object of class ", paste0(class(SEobj), collapse = "/"), " was supplied."))
    }

    JAMSversion <- metadata(SEobj)$version

    if (is.null(JAMSversion) || (length(JAMSversion) < 1) || all(is.na(JAMSversion))){
        stop(paste0(prefix, objname, " does not carry a JAMS version stamp in its metadata, which means it was built by JAMS 1.x (or by hand without one). Support for JAMS1 SummarizedExperiment objects was removed in JAMS 2.3.1, because their taxon stratification and genome completeness are stored on incompatible scales and silently produce wrong results. Please rebuild the object with a current JAMSbeta (or JAMS16, for amplicon data). If you are certain this object was in fact built by JAMS 2.x and has merely lost its stamp, you can restore it with: metadata(SEobj)$version <- paste0(\"JAMS v\", packageVersion(\"JAMS\"))"))
    }

    JAMSversion <- as.character(JAMSversion)[1]

    #Parse the leading integer of the version string, e.g. "JAMS v2.3.0" -> 2. If the string is
    #not parseable we accept it: the presence of the stamp at all is the primary test, and a
    #future reformatting of verstr must not start rejecting valid objects.
    majorversion <- suppressWarnings(as.integer(sub("^[^0-9]*([0-9]+).*$", "\\1", JAMSversion)))

    if (!is.na(majorversion) && (majorversion < 2)){
        stop(paste0(prefix, objname, " was built by ", JAMSversion, ". Support for JAMS1 SummarizedExperiment objects was removed in JAMS 2.3.1, because their taxon stratification and genome completeness are stored on incompatible scales and silently produce wrong results. Please rebuild the object from the .jams files with a current JAMSbeta."))
    }

    if (require_stratification){
        has_index <- "SparseIndex" %in% names(assays(SEobj))
        has_matrix <- "allfeaturesbytaxa_matrix" %in% names(metadata(SEobj))
        if (!all(c(has_index, has_matrix))){
            missingbits <- c("SparseIndex assay", "allfeaturesbytaxa_matrix metadata")[!c(has_index, has_matrix)]
            stop(paste0(prefix, objname, " does not carry feature-by-taxon stratification data (missing: ", paste0(missingbits, collapse = ", "), "). This information is only generated when JAMSbeta is run with the -u / --stratify_functions_by_taxon option, and only for functional analysis spaces."))
        }
    }

    return(invisible(JAMSversion))
}
