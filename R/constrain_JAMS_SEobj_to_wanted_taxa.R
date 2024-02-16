#' Constrains a functional (non-taxonomic) JAMS-style SummarizedExperiment object which has feature-by-taxon stratification information (see JAMSbeta -u option) into a SummarizedExperiment object with BaseCounts and GeneCounts only pertaining to the wanted LKTs.

#' @param SEobj a functional (non-taxonomic) JAMS-style SummarizedExperiment object which has feature-by-taxon stratification information (see JAMSbeta -u option)

#' @param wanted_LKTs vector of LKTs to constrain the functional features of the original SEobj into

#' @export

constrain_JAMS_SEobj_to_wanted_taxa <- function(SEobj = NULL, wanted_LKTs = NULL){

    #Determine the characteristics of ExpObj
    characterize_ExpObj <- function(ExpObj = NULL){
        ExpObj_facts <- list()
        ExpObj_facts$class <- class(ExpObj)[1]
        ExpObj_facts$analysis <- metadata(ExpObj)$analysis
        ExpObj_facts$assays <- sort(names(assays(ExpObj)))
        ExpObj_facts$metadata <- sort(names(metadata(ExpObj)))
        ExpObj_facts$dimensions <- dim(ExpObj)

        return(ExpObj_facts)
    }

    SEobj_facts <- characterize_ExpObj(SEobj)
    #Rebuild a sample-by-feature matrix with only the LKT.
    #First, the BaseCounts
    row_index <- metadata(SEobj)$allfeaturesbytaxa_index

    relevant_matrix <- as.matrix(metadata(SEobj)$allfeaturesbytaxa_matrix[ , wanted_LKTs])
    relevant_matrix <- as.matrix(rowSums(relevant_matrix))
    colnames(relevant_matrix)[1] <- "Count"
    relevant_matrix <- as.data.frame(relevant_matrix)
    relevant_matrix <- cbind(row_index, relevant_matrix)
    relevant_matrix$RowNumber <- NULL

    cts <- relevant_matrix %>% pivot_wider(names_from = Sample, values_from = Count, values_fill = 0)
    cts <- as.data.frame(cts)
    cts[is.na(cts)] <- 0
    rownames(cts) <- cts$Accession
    cts$Accession <- NULL

    #Deal with samples with NO results for this analysis
    emptySamples <- colnames(SEobj)[!(colnames(SEobj) %in% colnames(cts))]
    if (length(emptySamples) > 0){
        complementarycts <- matrix(ncol = length(emptySamples), nrow = nrow(cts), data = 0)
        complementarycts <- as.data.frame(complementarycts, stringsAsFactors = FALSE)
        colnames(complementarycts) <- emptySamples
        rownames(complementarycts) <- rownames(cts)
        cts <- cbind(cts, complementarycts)
    }

    cts <- as.matrix(cts)

    #Second, the GeneCounts
    if ("allfeaturesbytaxa_GeneCounts_matrix" %in% SEobj_facts$metadata){

        row_index_GC <- metadata(SEobj)$allfeaturesbytaxa_GeneCounts_index
        relevant_matrix_GC <- as.matrix(metadata(SEobj)$allfeaturesbytaxa_GeneCounts_matrix[ , wanted_LKTs])
        relevant_matrix_GC <- as.matrix(rowSums(relevant_matrix_GC))
        colnames(relevant_matrix_GC)[1] <- "Count"
        relevant_matrix_GC <- as.data.frame(relevant_matrix_GC)
        relevant_matrix_GC <- cbind(row_index_GC, relevant_matrix_GC)
        relevant_matrix_GC$RowNumber <- NULL

        featcts <- relevant_matrix_GC %>% pivot_wider(names_from = Sample, values_from = Count, values_fill = 0)
        featcts <- as.data.frame(featcts)
        featcts[is.na(featcts)] <- 0

        rownames(featcts) <- featcts$Accession
        featcts$Accession <- NULL

        #Deal with samples with NO results for this analysis
        emptySamples <- colnames(SEobj)[!(colnames(SEobj) %in% colnames(featcts))]
        if (length(emptySamples) > 0){
            complementarycts <- matrix(ncol = length(emptySamples), nrow = nrow(featcts), data = 0)
            complementarycts <- as.data.frame(complementarycts, stringsAsFactors = FALSE)
            colnames(complementarycts) <- emptySamples
            rownames(complementarycts) <- rownames(featcts)
            featcts <- cbind(featcts, complementarycts)
        }
        featcts$Accession <- NULL
        featcts <- as.matrix(featcts)
    }

    #Ensure features are the same
    SEobj <- SEobj[intersect(rownames(SEobj), rownames(cts)), ]
    cts <- cts[rownames(SEobj), ]
    featcts <- featcts[rownames(SEobj), ]

    assays <- list()
    assays$BaseCounts <- cts
    if (!is.null(featcts)){
        assays$GeneCounts <- featcts
    }

    phenoanal <- colData(SEobj)
    ftt <- rowData(SEobj)
    constrained_SEobj <- SummarizedExperiment(assays = assays, rowData = as.matrix(ftt), colData = as.matrix(phenoanal))

    #Put all the original bells and whistles back
    relevant_md <- SEobj_facts$metadata[!(SEobj_facts$metadata %in% c("allfeaturesbytaxa_GeneCounts_index", "allfeaturesbytaxa_GeneCounts_matrix", "allfeaturesbytaxa_index","allfeaturesbytaxa_matrix"))]
    for (md in relevant_md){
        metadata(constrained_SEobj)[[md]] <- metadata(SEobj)[[md]]
    }

    return(constrained_SEobj)
}
