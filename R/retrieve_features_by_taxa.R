#' retrieve_features_by_taxa(FuncExpObj = NULL, assay_for_matrix = "BaseCounts", wantedfeatures = NULL, wantedsamples = NULL, asPPM = TRUE, PPM_normalize_to_bases_sequenced = FALSE, PPMthreshold = 0)
#'
#' Returns a long form data frame of stratification by taxa of the relative abundance or number of bases wanted of functional features in wanted samples, given allfeaturesbytaxa_matrix and allfeaturesbytaxa_index metadata present in a JAMS SummarizedExperiment functional object.
#' @export

retrieve_features_by_taxa <- function(FuncExpObj = NULL, assay_for_matrix = "BaseCounts", wantedfeatures = NULL, wantedsamples = NULL, asPPM = TRUE, PPM_normalize_to_bases_sequenced = FALSE, PPMthreshold = 0){

    if (assay_for_matrix == "GeneCounts"){
        sparsematrix_name <- "allfeaturesbytaxa_GeneCounts_matrix"
        sparsematrix_index_name <- "allfeaturesbytaxa_GeneCounts_index"
    } else {
        sparsematrix_name <- "allfeaturesbytaxa_matrix"
        sparsematrix_index_name <- "allfeaturesbytaxa_index"
    }

    allfeaturesbytaxa_matrix <- metadata(FuncExpObj)[[sparsematrix_name]]
    allfeaturesbytaxa_index <- metadata(FuncExpObj)[[sparsematrix_index_name]]
    curr_pt <- colData(FuncExpObj)

    #Get appropriate rows
    rowsinterestdf <- subset(allfeaturesbytaxa_index, Sample %in% wantedsamples)
    rowsinterestdf <- subset(rowsinterestdf, Accession %in% wantedfeatures)

    allfeaturesbytaxa_interest <- as.matrix(allfeaturesbytaxa_matrix[rowsinterestdf$RowNumber, ])

    if (length(rowsinterestdf$RowNumber) == 1){
        allfeaturesbytaxa_interest <- t(allfeaturesbytaxa_interest)
    }

    allfeaturesbytaxa_interest <- as.data.frame(allfeaturesbytaxa_interest[, which(colSums(allfeaturesbytaxa_interest) != 0)])
    allfeaturesbytaxa_interest$RowNumber <- as.numeric(rownames(allfeaturesbytaxa_interest))

    allfeaturesbytaxa_interest <- left_join(allfeaturesbytaxa_interest, rowsinterestdf, by = "RowNumber")
    allfeaturesbytaxa_interest$RowNumber <- NULL
    allfeaturesbytaxa_interest <- allfeaturesbytaxa_interest[ , c("Sample", "Accession", (sort(colnames(allfeaturesbytaxa_interest)[which(!colnames(allfeaturesbytaxa_interest) %in% c("Sample", "Accession"))])))]

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
        for(colm in LKTcolumns){
            taxsplit[ , colm] <- round(((taxsplit[ , colm] / taxsplit$NumBases) * 1000000), 0)
        }

        #Denoise
        LKTsMaxima <-sapply(LKTcolumns, function(x) {max(taxsplit[ , x])})
        LKTsToKeep <- names(which(LKTsMaxima > PPMthreshold))
        sample2metadata <- as.data.frame(curr_pt)
        nonSamplecolms <- colnames(sample2metadata)[colnames(sample2metadata) != "Sample"]
        taxsplit <- left_join(taxsplit, sample2metadata, by = "Sample")
        taxsplit <- taxsplit[ , c("Sample", "Accession", "NumBases", nonSamplecolms, LKTsToKeep)]
        taxsplit$NumBases <- NULL

        return(taxsplit)

    } else {

        return(allfeaturesbytaxa_interest)

    }
}
