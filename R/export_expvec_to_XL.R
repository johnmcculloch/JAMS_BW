#' Exports counts and featuredata in a SummarizedExperiment vector into a single spreadsheet.

#' @param expvec list of JAMS-style SummarizedExperiment objects.

#' @param usefulexp String specifying the name(s) of the JAMS-style SummarizedExperiment object(s) to be exported into Excel. If NULL, will include all analysis contained in the expvec SummarizedExperiment list. Default is NULL.

#' @param filename String specifying the Excel file name. If NULL (the default), the filename will be automatically assigned.

#' @param featcutoff Requires a numeric vector of length 2 for specifying how to filter out features by relative abundance. The first value of the vector specifies the minimum relative abundance in Parts per Million (PPM) and the second value is the percentage of samples which must have at least that relative abundance. Thus, passing c(250, 10) to featcutoff would filter out any feature which does not have at least 250 PPM (= 0.025 percent) of relative abundance in at least 10 percent of all samples being plot. Please note that when using the subsetby option (q.v.) to automatically plot multiple plots of sample subsets, the featcutoff parameters are applied within the subset. The default is c(0, 0), meaning no feature is filtered. If NULL is passed, then the value defaults to c(0, 0). See also applyfilters for a shorthand way of applying multiple filtration settings.

#' @param GenomeCompletenessCutoff Requires a numeric vector of length 2 for specifying how to filter out features by genome completeness. This is, of course, only applicble for taxonomic shotgun SummarizedExperiment objects. When passed to non-taxonomic shotgun SummarizedExperiment objects, GenomeCompletenessCutoff will be ignored. The first value of the vector specifies the minimum genome completeness in percentage  and the second value is the percentage of samples which must have at least that genome completeness. Thus, passing c(50, 5) to GenomeCompletenessCutoff would filter out any taxonomic feature which does not have at least 50 percent of genome completeness in at least 5 percent of all samples being plot. Please note that when using the subsetby option (q.v.) to automatically plot multiple plots of sample subsets, the GenomeCompletenessCutoff parameters are applied within the subset. The default is c(0, 0), meaning no feature is filtered. If NULL is passed, then the value defaults to c(0, 0). See also applyfilters for a shorthand way of applying multiple filtration settings.

#' @param applyfilters Optional string specifying filtration setting "combos", used as a shorthand for setting the featcutoff, GenomeCompletenessCutoff, minl2fc and minabscorrcoeff arguments in JAMS plotting functions. If NULL, none of these arguments are set if not specified. Permissible values for applyfilters are "light", "moderate" or "stringent". The actual values vary whether the SummarizedExperiment object is taxonomical (LKT) or not. For a taxonomical SummarizedExperiment object, using "light" will set featcutoff=c(50, 5), GenomeCompletenessCutoff=c(5, 5), minl2fc=1, minabscorrcoeff=0.4; using "moderate" will set featcutoff=c(250, 15), GenomeCompletenessCutoff=c(10, 5), minl2fc=1, minabscorrcoeff=0.6; and using "stringent" will set featcutoff=c(2000, 15), GenomeCompletenessCutoff=c(30, 10), minl2fc=2, minabscorrcoeff=0.8. For non-taxonomical (i.e. functional) SummarizedExperiment objects, using "light" will set featcutoff=c(0, 0), minl2fc=1, minabscorrcoeff=0.4; using "moderate" will set featcutoff=c(5, 5), minl2fc=1, minabscorrcoeff=0.6; and using "stringent" will set featcutoff=c(50, 15), minl2fc=2.5, minabscorrcoeff=0.8. When using applyfilters, one can still set one or more of featcutoff, GenomeCompletenessCutoff, minl2fc and minabscorrcoeff, which will then take the user set value in lieu of those set by the applyfilters shorthand. Default is light.

#' @param PPM_normalize_to_bases_sequenced Requires a logical value. Non-filtered JAMS feature counts tables (the BaseCounts assay within SummarizedExperiment objects) always includes unclassified taxonomical features (for taxonomical SummarizedExperiment objects) or unknown/unattributed functional features (for non-taxonomical SummarizedExperiment objects), so the relative abundance for each feature (see normalization) will be calculated in Parts per Million (PPM) by dividing the number of bases covering each feature by the sum of each sample column **previous to any filtration**. Relative abundances are thus representative of the entirety of the genomic content for taxonomical objects, whereas for non-taxonomical objects, strictly speaking, it is the abundance of each feature relative to only the coding regions present in the metagenome, even if these are annotationally unatributed. In other words, intergenic regions are not taken into account. In order to relative-abundance-normalize a **non-taxonomical** SummarizedExperiment object with the total genomic sequencing content, including non-coding regions, set PPM_normalize_to_bases_sequenced = TRUE. Default is FALSE.

#' @param includemetadata Requires a logical value. If TRUE, information within colData of your inputted expvec, which contains all information pertaining to the samples, will be exported into a "Metadata" tab within the Excel sheet. If FALSE, will omit this information in the Excel. Default is TRUE.

#' @param returncounts Requires a logical value. If TRUE, will print the information to be exported to Excel onto the terminal (which can be directed into a variable). If the analysis is taxonomical, this information will include metadata (if desired), the relative abundance of each feature for each taxonomic feature in Parts Per Million(PPM) for each of the samples, the Genome Completeness for each taxonomic feature in each of the samples, the Percentage from Contigs for each taxonomic feature in each of the samples, and a featuretable containing the feature rownames and any associated taxonomic information. If the analysis is functional, this will include metadata (if desired), the relative abundance of each functional feature in the analysis in Parts Per Million(PPM) for each of the samples, the Genecounts consisting of the number of features annotated in the analysis, and a featuretable containing the feature rownames, an accession, and any addtional information on the feature. If FALSE (the default), the information described above will be exported to Excel only.

#' @export

export_expvec_to_XL <- function(expvec = NULL, usefulexp = NULL, filename = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, applyfilters = NULL, asPPM = TRUE, PPM_normalize_to_bases_sequenced = FALSE, includemetadata = TRUE, returncounts = FALSE,...){

   #Hardwire PctFromCtgscutoff, as this should never be used without filtering because of huge amounts of false positives when evaluating taxonomic information from unassembled reads. The use of classifying unassembled reads is deprecated in JAMS and the default is to NOT classify unassembled reads, so this is usually not an issue.
   PctFromCtgscutoff <- c(70, 50)

  #Adding a default filename.
  if ((is.null(filename)& !(is.null(usefulexp)))){
      filename <- paste0("JAMSexported", usefulexp, ".xlsx")
      flog.info("Customize file name with filename option.")
  } else if ((is.null(filename)& is.null(usefulexp))){
    filename <- "JAMSexportedexpvec.xlsx"
    flog.info("Customize file name with filename option, you utter fool")
  } else{
    filename <- filename
  }

    if (!(is.null(usefulexp))){
        usefulexp <- names(expvec)[(names(expvec) %in% usefulexp)]
        expvec2 <- expvec[usefulexp]
    } else {
        expvec2 <- expvec
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

    if (any(c((asPPM == TRUE), !is.null(applyfilters), !is.null(featcutoff), !is.null(GenomeCompletenessCutoff), !is.null(PctFromCtgscutoff)))){
        flog.info("Counts will be exported as PPM")
        countunits <- "PPM"
    } else {
        flog.info("Counts will be exported as raw number of bases")
        countunits <- "NumBases"
    }

    #Get counts
    for (x in 1:length(expvec2)){
        analysis <- metadata(expvec2[[x]])$analysis
        flog.info(paste("Exporting", analysis))

        presetlist <- declare_filtering_presets(analysis = analysis, applyfilters = applyfilters, featcutoff = featcutoff, GenomeCompletenessCutoff = GenomeCompletenessCutoff, PctFromCtgscutoff = PctFromCtgscutoff)

        exp_filt <- filter_experiment(ExpObj = expvec2[[x]], featcutoff = presetlist$featcutoff, asPPM = asPPM, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, GenomeCompletenessCutoff = presetlist$GenomeCompletenessCutoff, PctFromCtgscutoff = presetlist$PctFromCtgscutoff)

        cts <- assays(exp_filt)$BaseCounts
        ctsname <- paste(names(expvec2)[x], countunits, sep="_")
        countvec[[cvn]] <- as.data.frame(cts)
        names(countvec)[cvn] <- ctsname
        cvn <- cvn + 1

        if ("GenomeCompleteness" %in% names(assays((expvec2)[[x]]))){
            cts <- assays(exp_filt)$GenomeCompleteness
            ctsname <- paste(names(expvec2)[x], "GnmCompl", sep="_")
            countvec[[cvn]] <- as.data.frame(cts)
            names(countvec)[cvn] <- ctsname
            cvn <- cvn + 1
        }

        if ("PctFromCtgs" %in% names(assays((expvec2)[[x]]))){
            cts <- assays(exp_filt)$PctFromCtgs
            ctsname <- paste(names(expvec2)[x], "PctFromCtgs", sep="_")
            countvec[[cvn]] <- as.data.frame(cts)
            names(countvec)[cvn] <- ctsname
            cvn <- cvn + 1
        }

        if ("GeneCounts" %in% names(assays((expvec2)[[x]]))){
            cts <- assays(exp_filt)$GeneCounts
            ctsname <- paste(names(expvec2)[x], "GeneCounts", sep="_")
            countvec[[cvn]] <- as.data.frame(cts)
            names(countvec)[cvn] <- ctsname
            cvn <- cvn + 1
        }

        feattbl <- as.data.frame(rowData(expvec2[[x]]))
        ftsname <- paste(names(expvec2)[x], "featuretable", sep="_")
        countvec[[cvn]] <- feattbl
        names(countvec)[cvn] <- ftsname
        cvn <- cvn + 1
    }

    countvec <- countvec[sapply(countvec, function(x){ !(is.null(x)) })]
    #Cap tab names to 31 characters because of excel
    if(any(nchar(names(countvec)) > 30)){
        flog.info("Capping tab names to 31characters because of Microsoft Excel.")
        newnames <- sapply(names(countvec), function(x){ stringr::str_trunc(x, 28) })
        newnames <- make.unique(newnames)
        names(countvec) <- newnames
    }

    flog.info(paste("Saving spreadsheet as", filename))
    write.xlsx(countvec, file = filename, asTable = TRUE, rowNames = TRUE, colNames = TRUE, borders = "all", colWidths = "auto")
    if (returncounts){
        return(countvec)
    }
}
