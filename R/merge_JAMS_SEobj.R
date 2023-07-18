#' Safely merges two compatible JAMS-style SummarizedExperiment objects into a single JAMS-style SummarizedExperiment object. SE objects must belong to the same analysis type, and have no Sample names in common. Ideally they should have been generated with the same version of JAMS using the same JAMS database. If the latter condition is not met, the consequences are unforseen and you are ill-advised to proceed.

#' @param ExpObj1 a JAMS-style SummarizedExperiment object

#' @param ExpObj2 a JAMS-style SummarizedExperiment object

#' @export

merge_JAMS_SEobj <- function(ExpObj1 = NULL, ExpObj2 = NULL){

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

    ExpObj1_facts <- characterize_ExpObj(ExpObj1)
    ExpObj2_facts <- characterize_ExpObj(ExpObj2)
    Samples_in_common <- intersect(colnames(ExpObj1), colnames(ExpObj2))
    Features_union <- union(rownames(ExpObj1), rownames(ExpObj2))

    #Check we're merging the same kind of ExpObj
    if (ExpObj1_facts$analysis != ExpObj1_facts$analysis){
        stop(flog.warn(paste("ExpObj1 is a", ExpObj1_facts$analysis, "SummarizedExperimentE object, whereas ExpObj2 is a", ExpObj2_facts$analysis, "SummarizedExperiment object, and therefore cannot be merged. Please use same analysis kind as input for merging SummarizedExperiment objects.")))
    }

    #Check we're not merging ExpObjs with the same samples
    if (length(Samples_in_common) > 0){
        stop(flog.warn(paste("ExpObj1 is and ExpObj2 have samples in common. Please use SummarizedExperiment objects with different sample names for merging.\nSamples in common are:", paste0(Samples_in_common, collapse = ", "))))
    }

    #Merge count matrices
    merged_assays <- list()
    for (ass in intersect(ExpObj1_facts$assays, ExpObj2_facts$assays)){
        flog.info(paste("Merging", ass))
        cts1 <- as.data.frame(assays(ExpObj1)[[ass]])
        cts2 <- as.data.frame(assays(ExpObj2)[[ass]])
        cts1$FeAtUrEnAmE <- rownames(cts1)
        cts2$FeAtUrEnAmE <- rownames(cts2)
        merged_cts <- base::merge(cts1, cts2, by = "FeAtUrEnAmE", all = TRUE)
        rownames(merged_cts) <- merged_cts$FeAtUrEnAmE
        merged_cts$FeAtUrEnAmE <- NULL
        merged_cts[is.na(merged_cts)] <- 0
        merged_cts <- as.matrix(merged_cts)
        merged_assays[[ass]] <- merged_cts
    }

    feat_order <- rownames(merged_assays[["BaseCounts"]])[order(rowSums(merged_assays[["BaseCounts"]]), decreasing = TRUE)]
    for (ass in intersect(ExpObj1_facts$assays, ExpObj2_facts$assays)){
        merged_assays[[ass]] <- merged_assays[[ass]][feat_order, ]
    }

    #Merge feature table(s)
    flog.info("Merging feature tables")
    ftt1 <- as.data.frame(rowData(ExpObj1))
    ftt2 <- as.data.frame(rowData(ExpObj2))
    ftt1$FeAtUrEnAmE <- rownames(ftt1)
    ftt2$FeAtUrEnAmE <- rownames(ftt2)
    #Not optimal, because there could be annotation conflict, but the user is responsible for merging compatible SE objects, so will keep it like this for the time being.
    common_cols <- intersect(colnames(ftt1), colnames(ftt2))
    ftt1 <- ftt1[ , common_cols]
    ftt2 <- ftt2[ , common_cols]
    merged_ftt <- rbind(ftt1, ftt2)
    merged_ftt <- merged_ftt[!duplicated(merged_ftt), ]
    #This should be the non-redundant feature table, but check and warn if there are discrepancies
    if (is.redundant(merged_ftt$FeAtUrEnAmE)){
        numdupes <- length(which(duplicated(merged_ftt$FeAtUrEnAmE) == TRUE))
        flog.warn(paste("ATTENTION: a total of", numdupes, "out of", nrow(merged_ftt), "unique feature table rows are divergent. This can be due to different annotation of features between the merging SummarizedExperiment objects. Making the feature names unique."))
        merged_ftt <- merged_ftt[!duplicated(merged_ftt$FeAtUrEnAmE), ]
    }
    rownames(merged_ftt) <- merged_ftt$FeAtUrEnAmE
    merged_ftt$FeAtUrEnAmE <- NULL
    #Should be the same, but enforce anyway
    merged_ftt <- merged_ftt[(rownames(merged_ftt) %in% feat_order), ]
    merged_ftt <- merged_ftt[feat_order, ]

    #Merge the metadata
    flog.info("Merging phenotable")
    pheno1 <- as.data.frame(colData(ExpObj1))
    pheno2 <- as.data.frame(colData(ExpObj2))
    pheno_merged <- base::merge(pheno1, pheno2, all = TRUE)
    rownames(pheno_merged) <- pheno_merged$Sample
    #Fill in the gaps with JAMS-style NA
    for (colm in 1:ncol(pheno_merged)){
        pheno_merged[is.na(pheno_merged[ , colm]), colm] <- "N_A"
    }
    pheno_merged <- pheno_merged[colnames(merged_assays[["BaseCounts"]]), ]

    #Merge the remaning bells and whistles (SEobj "metadata")
    SEobjMD_merged <- list()
    for (SEobjMD in intersect(ExpObj1_facts$metadata, ExpObj2_facts$metadata)){
        if (SEobjMD == "analysis"){
            #should be the same, but I am neurotic.
            SEobjMD_merged[[SEobjMD]] <- unique(metadata(ExpObj1)$analysis, metadata(ExpObj1)$analysis)[1]
        } else if (SEobjMD %in% c("TotalBasesSequenced", "TotalBasesSequencedinAnalysis")) {
            SEobjMD_merged[[SEobjMD]] <- cbind(metadata(ExpObj1)[[SEobjMD]], metadata(ExpObj2)[[SEobjMD]])
            #Ensure order, because, again, I am neurotic
            SEobjMD_merged[[SEobjMD]] <- t(SEobjMD_merged[[SEobjMD]]["NumBases", colnames(merged_assays[["BaseCounts"]])])
            rownames(SEobjMD_merged[[SEobjMD]]) <- "NumBases"
        } else if (SEobjMD == "ctable") {
            SEobjMD_merged[[SEobjMD]] <- base::merge(metadata(ExpObj1)[[SEobjMD]], metadata(ExpObj2)[[SEobjMD]], all = TRUE)
            SEobjMD_merged[[SEobjMD]] <- SEobjMD_merged[[SEobjMD]][!duplicated(SEobjMD_merged[[SEobjMD]]$Name), ]
            rownames(SEobjMD_merged[[SEobjMD]]) <- SEobjMD_merged[[SEobjMD]]$Name
        } else {
            SEobjMD_merged[[SEobjMD]] <- NULL
        }
    }

    #Build the merged SummarizedExperiment object
    merged_SEobj <- SummarizedExperiment(assays = merged_assays, rowData = as.matrix(merged_ftt), colData = as.matrix(pheno_merged), metadata = SEobjMD_merged)

    return(merged_SEobj)
}
