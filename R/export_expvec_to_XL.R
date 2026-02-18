#' Exports counts and featuredata in a SummarizedExperiment vector into a single spreadsheet.

#' @param expvec list of JAMS-style SummarizedExperiment objects.

#' @param expvec_analysis_spaces String specifying the name(s) of the JAMS-style SummarizedExperiment object(s) to be exported into Excel. If NULL, will include all analysis contained in the expvec SummarizedExperiment list. Default is NULL.

#' @param SEobj single JAMS-style SummarizedExperiment object. If passed, parameter expvec will be ignored.

#' @param filename String specifying the Excel file name. If NULL (the default), the filename will be automatically generated. Be mindful that in this case, the filename might be long and very descriptive.

#' @param featcutoff Requires a numeric vector of length 2 for specifying how to filter out features by relative abundance. The first value of the vector specifies the minimum relative abundance in Parts per Million (PPM) and the second value is the percentage of samples which must have at least that relative abundance. Thus, passing c(250, 10) to featcutoff would filter out any feature which does not have at least 250 PPM (= 0.025 percent) of relative abundance in at least 10 percent of all samples being plot. Please note that when using the subsetby option (q.v.) to automatically plot multiple plots of sample subsets, the featcutoff parameters are applied within the subset. The default is c(0, 0), meaning no feature is filtered. If NULL is passed, then the value defaults to c(0, 0). See also applyfilters for a shorthand way of applying multiple filtration settings.

#' @param GenomeCompletenessCutoff Requires a numeric vector of length 2 for specifying how to filter out features by genome completeness. This is, of course, only applicble for taxonomic shotgun SummarizedExperiment objects. When passed to non-taxonomic shotgun SummarizedExperiment objects, GenomeCompletenessCutoff will be ignored. The first value of the vector specifies the minimum genome completeness in percentage  and the second value is the percentage of samples which must have at least that genome completeness. Thus, passing c(50, 5) to GenomeCompletenessCutoff would filter out any taxonomic feature which does not have at least 50 percent of genome completeness in at least 5 percent of all samples being plot. Please note that when using the subsetby option (q.v.) to automatically plot multiple plots of sample subsets, the GenomeCompletenessCutoff parameters are applied within the subset. The default is c(0, 0), meaning no feature is filtered. If NULL is passed, then the value defaults to c(0, 0). See also applyfilters for a shorthand way of applying multiple filtration settings.

#' @param applyfilters Optional string specifying filtration setting "combos", used as a shorthand for setting the featcutoff, GenomeCompletenessCutoff, minl2fc and minabscorrcoeff arguments in JAMS plotting functions. If NULL, none of these arguments are set if not specified. Permissible values for applyfilters are "light", "moderate" or "stringent". The actual values vary whether the SummarizedExperiment object is taxonomical (LKT) or not. For a taxonomical SummarizedExperiment object, using "light" will set featcutoff=c(50, 5), GenomeCompletenessCutoff=c(5, 5), minl2fc=1, minabscorrcoeff=0.4; using "moderate" will set featcutoff=c(250, 15), GenomeCompletenessCutoff=c(10, 5), minl2fc=1, minabscorrcoeff=0.6; and using "stringent" will set featcutoff=c(2000, 15), GenomeCompletenessCutoff=c(30, 10), minl2fc=2, minabscorrcoeff=0.8. For non-taxonomical (i.e. functional) SummarizedExperiment objects, using "light" will set featcutoff=c(0, 0), minl2fc=1, minabscorrcoeff=0.4; using "moderate" will set featcutoff=c(5, 5), minl2fc=1, minabscorrcoeff=0.6; and using "stringent" will set featcutoff=c(50, 15), minl2fc=2.5, minabscorrcoeff=0.8. When using applyfilters, one can still set one or more of featcutoff, GenomeCompletenessCutoff, minl2fc and minabscorrcoeff, which will then take the user set value in lieu of those set by the applyfilters shorthand. Default is light.

#' @param PPM_normalize_to_bases_sequenced Requires a logical value. Non-filtered JAMS feature counts tables (the BaseCounts assay within SummarizedExperiment objects) always includes unclassified taxonomical features (for taxonomical SummarizedExperiment objects) or unknown/unattributed functional features (for non-taxonomical SummarizedExperiment objects), so the relative abundance for each feature (see normalization) will be calculated in Parts per Million (PPM) by dividing the number of bases covering each feature by the sum of each sample column **previous to any filtration**. Relative abundances are thus representative of the entirety of the genomic content for taxonomical objects, whereas for non-taxonomical objects, strictly speaking, it is the abundance of each feature relative to only the coding regions present in the metagenome, even if these are annotationally unatributed. In other words, intergenic regions are not taken into account. In order to relative-abundance-normalize a **non-taxonomical** SummarizedExperiment object with the total genomic sequencing content, including non-coding regions, set PPM_normalize_to_bases_sequenced = TRUE. Default is FALSE.

#' @param includemetadata Requires a logical value. If TRUE, information within colData of your inputted expvec, which contains all information pertaining to the samples, will be exported into a "Metadata" tab within the Excel sheet. If FALSE, will omit this information in the Excel. Default is TRUE.

#' @param returncounts Requires a logical value. If TRUE, will print the information to be exported to Excel onto the terminal (which can be directed into a variable). If the analysis is taxonomical, this information will include metadata (if desired), the relative abundance of each feature for each taxonomic feature in Parts Per Million(PPM) for each of the samples, the Genome Completeness for each taxonomic feature in each of the samples, the Percentage from Contigs for each taxonomic feature in each of the samples, and a featuretable containing the feature rownames and any associated taxonomic information. If the analysis is functional, this will include metadata (if desired), the relative abundance of each functional feature in the analysis in Parts Per Million(PPM) for each of the samples, the Genecounts consisting of the number of features annotated in the analysis, and a featuretable containing the feature rownames, an accession, and any addtional information on the feature. If FALSE (the default), the information described above will be exported to Excel only.

#' @export

export_expvec_to_XL <- function(expvec = NULL, expvec_analysis_spaces = NULL, SEobj = NULL, filename = NULL, samplesToKeep = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, applyfilters = NULL, normalization = c("relabund", "clr"), PPM_normalize_to_bases_sequenced = TRUE, includemetadata = TRUE, returncounts = FALSE, ...){

    if (all(c(is.null(expvec), is.null(SEobj)))) {
        stop("You must choose as input for exporting either a list of SummarizedExperiment objects through the expvec argument, or a single SEobj via the SEobj argument.")
    }

    if (!is.null(SEobj)){
        expvec2 <- list()
        expvec2[[1]] <- SEobj
        names(expvec2)[1] <- metadata(SEobj)$analysis
    } else {
        if (!(is.null(expvec_analysis_spaces))){
            usefulexp <- names(expvec)[(names(expvec) %in% expvec_analysis_spaces)]
            expvec2 <- expvec[usefulexp]
        } else {
            expvec2 <- expvec
        }
    }

    #Set very descriptive filename if not passed to argument
    if (is.null(filename)){
        filename <- generate_filename(title = paste0(c("Counts_from_SummarizedExperiments", names(expvec2), "basecounts", normalization), collapse = "_"), add_date = TRUE, suffix = "xlsx")
    }

    pt <- as.data.frame(colData(expvec2[[1]]))

    countvec <- NULL
    countvec <- list()
    cvn <- 1
    if (includemetadata == TRUE){
        flog.info("Exporting metadata.")
        countvec[[cvn]] <- pt
        names(countvec)[cvn] <- "Metadata"
        cvn <- cvn + 1
    }

    if (is.null(featcutoff)){
        featcutoff <- c(0,0)
    }

    if (is.null(GenomeCompletenessCutoff)){
        GenomeCompletenessCutoff <- c(0,0)
    }

    expvec3 <- lapply(names(expvec2), function (x) { filter_experiment(SEobj = expvec2[[x]], samplesToKeep = samplesToKeep, featuresToKeep = NULL, featcutoff = featcutoff, GenomeCompletenessCutoff = GenomeCompletenessCutoff, applyfilters = applyfilters, normalization = normalization, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, flush_out_empty_samples = FALSE, clr_pseudocount = 1, give_info = FALSE) } )
    names(expvec3) <- names(expvec2)

    #Obtain list of possible assay names
    valid_countstables <- c("BaseCounts", "PPM", "CLR", "GenomeCompleteness", "GenomeContamination", "GeneCounts", "GeneLengths")
    valid_countstables <- valid_countstables[valid_countstables %in% unique(unname(unlist(sapply(names(expvec3), function (x) { names(assays(expvec3[[x]])) }))))]

    #Get tables
    for (analnumb in 1:length(expvec3)){
        analysis <- metadata(expvec3[[analnumb]])$analysis
        flog.info(paste("Exporting", analysis))

        #Export counts tables present in SEobj
        present_countstables <- valid_countstables[valid_countstables %in% names(assays(expvec3[[analnumb]]))]
        for (ct in present_countstables){
            cts <- assays(expvec3[[analnumb]])[[ct]]
            ctsname <- paste(analysis, ct, sep="_")
            countvec[[cvn]] <- as.data.frame(cts)
            names(countvec)[cvn] <- ctsname
            cvn <- cvn + 1
        }

        #Export feature tables present in SEobj
        feattbl <- as.data.frame(rowData(expvec3[[analnumb]]))
        ftsname <- paste(analysis, "featuretable", sep="_")
        countvec[[cvn]] <- feattbl
        names(countvec)[cvn] <- ftsname
        cvn <- cvn + 1
    }

    countvec <- countvec[sapply(countvec, function(x){ !(is.null(x)) })]
    #abbreviate Consolidated Genome Bin to CGB because of Excel's constraints
    names(countvec) <- gsub("ConsolidatedGenomeBin", "CGB", names(countvec))

    #Cap tab names to 31 characters because of excel
    if (any(nchar(names(countvec)) > 30)){
        flog.info("Capping tab names to 31 characters because of Microsoft Excel.")
        newnames <- sapply(names(countvec), function(x){ stringr::str_trunc(x, 29) })
        newnames <- make.unique(newnames)
        names(countvec) <- newnames
    }

    flog.info(paste("Saving spreadsheet as", filename))
    write.xlsx(countvec, file = filename, asTable = TRUE, rowNames = TRUE, colNames = TRUE, borders = "all", colWidths = "auto")
    if (returncounts){
        return(countvec)
    }
}
