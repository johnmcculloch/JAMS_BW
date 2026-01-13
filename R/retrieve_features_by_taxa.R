#' retrieve_features_by_taxa(FuncExpObj = NULL, glomby = NULL, assay_for_matrix = "BaseCounts", wantedfeatures = NULL, wantedsamples = NULL, asPPM = TRUE, PPM_normalize_to_bases_sequenced = FALSE, PPMthreshold = 0)
#'
#' Returns a long form data frame of stratification by taxa of the relative abundance or number of bases wanted of functional features in wanted samples, given allfeaturesbytaxa_matrix and allfeaturesbytaxa_index metadata present in a JAMS SummarizedExperiment functional object.
#' @export

retrieve_features_by_taxa <- function(FuncExpObj = NULL, glomby = NULL, assay_for_matrix = "BaseCounts", wantedfeatures = NULL, wantedsamples = NULL, asPPM = TRUE, PPM_normalize_to_bases_sequenced = FALSE, PPMthreshold = 0){

    if (assay_for_matrix == "GeneCounts"){
        sparsematrix_name <- "allfeaturesbytaxa_GeneCounts_matrix"
        sparsematrix_index_name <- "allfeaturesbytaxa_GeneCounts_index"
    } else {
        sparsematrix_name <- "allfeaturesbytaxa_matrix"
        sparsematrix_index_name <- "allfeaturesbytaxa_index"
    }

    #Default to all samples and/or if wantedsamples/wantedfeatures is not passed
    if (is.null(wantedsamples)){
        wantedsamples <- colnames(FuncExpObj)
    } else {
        wantedsamples <- wantedsamples[wantedsamples %in% colnames(FuncExpObj)]
    }

    if (is.null(wantedfeatures)){
        wantedfeatures <- rownames(FuncExpObj)
    } else {
        wantedfeatures <- wantedfeatures[wantedfeatures %in% rownames(FuncExpObj)]
    }

    #Determine first if FuncExpObj is a JAMS version > 2.1 object
    if ("SparseIndex" %in% names(assays(FuncExpObj))){
        JAMS2_object <- TRUE
        index_matrix <- as.data.frame(assays(FuncExpObj)[["SparseIndex"]])
        #Subset to only wanted features and samples
        index_matrix <- index_matrix[wantedfeatures, wantedsamples, drop = FALSE]
        rowsinterestdf <- as.data.frame(index_matrix) %>% rownames_to_column("Accession") %>% pivot_longer(cols = -Accession, names_to = "Sample", values_to = "RowNumber")
        rowsinterestdf <- as.data.frame(rowsinterestdf)

    } else {
        JAMS2_object <- FALSE
        allfeaturesbytaxa_index <- metadata(FuncExpObj)[[sparsematrix_index_name]]
        #Get appropriate rows
        rowsinterestdf <- subset(allfeaturesbytaxa_index, Sample %in% wantedsamples)
        rowsinterestdf <- subset(rowsinterestdf, Accession %in% wantedfeatures)
    }

    #Ensure correct class
    rowsinterestdf$RowNumber <- as.integer(rowsinterestdf$RowNumber)
    #Eliminate unavailable information.
    rowsinterestdf <- rowsinterestdf[!is.na(rowsinterestdf$RowNumber), ]

    allfeaturesbytaxa_interest <- as.matrix(metadata(FuncExpObj)[[sparsematrix_name]][rowsinterestdf$RowNumber, , drop = FALSE])
    #Prune empty taxa, i.e. taxa which have 0 counts for all of the features of interest
    allfeaturesbytaxa_interest <- as.data.frame(allfeaturesbytaxa_interest[, which(colSums(allfeaturesbytaxa_interest) != 0)])
    allfeaturesbytaxa_interest$RowNumber <- as.integer(rownames(allfeaturesbytaxa_interest))
    allfeaturesbytaxa_interest <- left_join(allfeaturesbytaxa_interest, rowsinterestdf, by = "RowNumber")
    allfeaturesbytaxa_interest$RowNumber <- NULL
    allfeaturesbytaxa_interest <- allfeaturesbytaxa_interest[ , c("Sample", "Accession", (sort(colnames(allfeaturesbytaxa_interest)[which(!colnames(allfeaturesbytaxa_interest) %in% c("Sample", "Accession"))])))]

    #Agglomerate if applicable
    if (!is.null(glomby)){
        #Check this feature is available before proceeding
        if (!JAMS2_object){
            stop("Stratification of functions by taxonomic levels higher than LKT is only available for SummarizedExperiment objects made with JAMS version > 2.1.0 ")
        }
        #Obtain taxonomy table
        stt <- as.data.frame(metadata(FuncExpObj)$tt)
        LKTs_to_glom <- (sort(colnames(allfeaturesbytaxa_interest)[which(!colnames(allfeaturesbytaxa_interest) %in% c("Sample", "Accession", "Ultra_low_abundance_LKTs"))]))
        stt <- stt[which(LKTs_to_glom %in% stt$LKT), c(glomby, "LKT")]
        #Fill in unclassifieds
        stt[which(is.na(stt[ , glomby])), glomby] <- "Unclassified"
        #extract matrix to glom
        LKTs_to_glom_mat <- allfeaturesbytaxa_interest[ , c("Sample", "Accession", LKTs_to_glom)]
        #Filter LKTs present in stt
        common_LKTs <- intersect(colnames(LKTs_to_glom_mat), stt$LKT)
        #Keep only those LKT columns
        LKTs_filtered <- LKTs_to_glom_mat[ , common_LKTs, drop = FALSE]
        # Reshape LKTs_to_glom_mat to long format
        LKTs_to_glom_long <- LKTs_filtered %>%  mutate(row_id = row_number()) %>%  pivot_longer(cols = -row_id, names_to = "LKT", values_to = "NumBases")
        # Join with stt to get glomby level for each LKT
        LKTs_to_glom_long <- LKTs_to_glom_long %>% left_join(stt, by = "LKT")
        #Rename relevant glomby column to make things more practical
        colnames(LKTs_to_glom_long)[which(colnames(LKTs_to_glom_long) == glomby)] <- "Glomby"
        #Aggregate counts per row_id × Glomby
        Glomby_sum <- LKTs_to_glom_long %>% group_by(row_id, Glomby) %>% summarise(NumBases = sum(NumBases, na.rm = TRUE), .groups = "drop")
        #Spread back to wide format (one row per original row)
        Glomby_wide <- Glomby_sum %>% pivot_wider(names_from = Glomby, values_from = NumBases, values_fill = 0)
        Glomby_wide <- as.data.frame(Glomby_wide)
        rownames(Glomby_wide) <- Glomby_wide$row_id
        Glomby_wide$row_id <- NULL
        Glom_taxa <- colnames(Glomby_wide)
        Glomby_wide[ , "Sample"] <- allfeaturesbytaxa_interest[rownames(Glomby_wide), "Sample"]
        Glomby_wide[ , "Accession"] <- allfeaturesbytaxa_interest[rownames(Glomby_wide), "Accession"]
        #Rearrange order
        allfeaturesbytaxa_interest <- Glomby_wide[ , c("Sample", "Accession", Glom_taxa)]
    }

    if (asPPM){
        taxsplit <- allfeaturesbytaxa_interest
        #Transform to PPM
        if (PPM_normalize_to_bases_sequenced == TRUE){
            totbases <- "TotalBasesSequenced"
        } else {
            totbases <- "TotalBasesSequencedinAnalysis"
        }
        numbases2sampl <- as.data.frame(t(metadata(FuncExpObj)[[totbases]]))
        numbases2sampl$Sample <- rownames(numbases2sampl)
        taxsplit <- left_join(taxsplit, numbases2sampl, by = "Sample")
        LKTcolumns <- colnames(taxsplit)[!(colnames(taxsplit) %in% c("Sample", "Accession", "NumBases"))]

        #Transform to PPM
        for (colm in LKTcolumns){
            taxsplit[ , colm] <- round(((taxsplit[ , colm] / taxsplit$NumBases) * 1000000), 0)
        }

        #Denoise
        LKTsMaxima <- sapply(LKTcolumns, function(x) {max(taxsplit[ , x])})
        LKTsToKeep <- names(which(LKTsMaxima > PPMthreshold))
        sample2metadata <- as.data.frame(colData(FuncExpObj))
        nonSamplecolms <- colnames(sample2metadata)[colnames(sample2metadata) != "Sample"]
        taxsplit <- left_join(taxsplit, sample2metadata, by = "Sample")
        taxsplit <- taxsplit[ , c("Sample", "Accession", "NumBases", nonSamplecolms, LKTsToKeep)]
        taxsplit$NumBases <- NULL

        return(taxsplit)

    } else {

        return(allfeaturesbytaxa_interest)

    }
}
