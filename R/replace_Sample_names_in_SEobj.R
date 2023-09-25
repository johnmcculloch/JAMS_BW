#' replace_Sample_names_in_SEobj(curr_SEobj = NULL, Sample_name_mapping = NULL, oldname_column = "Sample", newname_column = NULL)
#'
#' Renames samples in JAMS SummarizedExperiment objects safely.
#' Use with care and at your own risk

#' @export

replace_Sample_names_in_SEobj <- function(curr_SEobj = NULL, Sample_name_mapping = NULL, oldname_column = "Sample", newname_column = NULL){

    #Fix pheno
    curr_pt <- as.data.frame(colData(curr_SEobj))
    samples_present <- intersect(Sample_name_mapping[ , oldname_column], rownames(curr_pt))
    curr_pt <- curr_pt[samples_present, ]
    curr_pt_newnames <- Sample_name_mapping[match(rownames(curr_pt), Sample_name_mapping[ , oldname_column]), newname_column]
    rownames(curr_pt) <- curr_pt_newnames
    if ("Sample" %in% colnames(curr_pt)){
        curr_pt[ , "Sample"] <- rownames(curr_pt)
    }

    #Get feature table - nothing to change
    ft <- rowData(curr_SEobj)

    #Fix matrices
    matrices_have <- names(assays(curr_SEobj))[names(assays(curr_SEobj)) %in% c("BaseCounts", "GeneCounts", "PctFromCtgs", "GenomeCompleteness")]

    matrices_have_list <- list()
    for (mx in matrices_have){
        curr_mx <- as.data.frame(assays(curr_SEobj)[[mx]])
        samples_present <- intersect(Sample_name_mapping[ , oldname_column], colnames(curr_mx))
        curr_mx <- curr_mx[ , samples_present]
        curr_mx_newnames <- Sample_name_mapping[match(colnames(curr_mx), Sample_name_mapping[ , oldname_column]), newname_column]
        colnames(curr_mx) <- curr_mx_newnames
        matrices_have_list[[mx]] <- as.matrix(curr_mx)
    }

    new_SEobj <- SummarizedExperiment(assays = matrices_have_list, rowData = ft, colData = curr_pt)

    #Fix metadata
    metadata_have <- names(metadata(curr_SEobj))[names(metadata(curr_SEobj)) %in% c("TotalBasesSequenced", "TotalBasesSequencedinAnalysis")]
    for (md in metadata_have){
        curr_md <- as.data.frame(metadata(curr_SEobj)[[md]])
        samples_present <- intersect(Sample_name_mapping[ , oldname_column], colnames(curr_md))
        curr_md <- curr_md["NumBases" , samples_present]
        curr_md_newnames <- Sample_name_mapping[match(colnames(curr_md), Sample_name_mapping[ , oldname_column]), newname_column]
        colnames(curr_md) <-  curr_md_newnames
        metadata(new_SEobj)[[md]] <- curr_md
    }
    metadata(new_SEobj)[["analysis"]] <- metadata(curr_SEobj)[["analysis"]]

    return(new_SEobj)
}
