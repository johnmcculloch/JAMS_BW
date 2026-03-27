#' agglomerate_features(ExpObj = NULL, glomby = NULL)
#'
#' Agglomerates features in a JAMS SummarizedExperiment object safely
#' @export

agglomerate_features <- function(ExpObj = NULL, glomby = NULL){

    #Get appropriate object to work with
    if (as.character(class(ExpObj)[1]) != "SummarizedExperiment"){
        stop("This function can only take a SummarizedExperiment object as input.")
    }

    terminal_taxonomic_levels <- c("Contig_LKT", "MB2bin", "ConsolidatedGenomeBin", "LKT")
    #Test for silly requests
    if (glomby %in% terminal_taxonomic_levels){
        flog.warn(paste("No need to agglomerate by", glomby, "as it is a terminal level. Returning original SEobj."))

        return(ExpObj)
    }

    #Find out what kind of an object it is
    analysis <- metadata(ExpObj)$analysis
    if (analysis %in% c("Contig_LKT", "MB2bin", "ConsolidatedGenomeBin", "16S")){
        JAMS2TaxonomicObj <- TRUE
    } else {
        JAMS2TaxonomicObj <- FALSE
    }

    #Get feature table
    ftt <- as.data.frame(rowData(ExpObj))
    if (!(glomby %in% colnames(ftt))){
        stop(paste("Unable to agglomerate by", glomby, "because this category was not found in the feature table of the SummarizedExperiment object."))
    }

    pheno_original <- colData(ExpObj)
    taxonomic_levels <- c("Gram", "Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "IS1", "LKT")
    taxonomic_levels_tags <- c("", "d__", "k__", "p__", "c__", "o__", "f__", "g__", "s__", "is1__", "LKT__")
    #Deal with pure unclassifieds in Kingdom to Species
    if (glomby %in% c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")){
        #Reallocate LKT to tag__Unclassified to avoid glomming different LKTs at that level into the same glom
        appropriate_tag <- taxonomic_levels_tags[which(taxonomic_levels == glomby)]
        ftt[which(ftt[ , glomby] == paste0(appropriate_tag, "Unclassified")), glomby] <- ftt[which(ftt[ , glomby] == paste0(appropriate_tag, "Unclassified")), "LKT"]
    }

    assays <- list()
    #Aggregate counts by summing
    cts <- as.data.frame(assays(ExpObj)$BaseCounts)
    cts$Feats <- rownames(cts)
    feats2glomby_feats <- data.frame(Feats = rownames(ftt), Glomby_feats = as.character(ftt[ , glomby]), stringsAsFactors = FALSE)
    #Account for completely dark matter (i.e. LKT__Unclassified) not having a taxonomic path.
    if ("LKT__Unclassified" %in% feats2glomby_feats$Feats){
        #Attribute LKT__Unclassified to that glomby
        feats2glomby_feats[which(feats2glomby_feats$Feats == "LKT__Unclassified"), "Glomby_feats"] <- "LKT__Unclassified"
    }
    cts <- left_join(cts, feats2glomby_feats, by = "Feats")
    cts$Feats <- NULL

    glom_cts <- aggregate(. ~ Glomby_feats, data = cts, FUN = sum)
    rownames(glom_cts) <- glom_cts$Glomby_feats
    glom_cts$Glomby_feats <- NULL
    #Reorder from more to less prevalent
    glom_cts <- glom_cts[order(rowSums(glom_cts), decreasing = TRUE), ]
    featureorder <- rownames(glom_cts)
    sampleorder <- rownames(pheno_original)
    #Check everything is in the same order
    glom_cts <- glom_cts[, sampleorder]
    assays[["BaseCounts"]] <- glom_cts

    for (CurrAssay in c("GenomeCompleteness", "GenomeContamination", "GeneCounts", "GeneLengths")){
        if (CurrAssay %in% as.character(names(assays(ExpObj)))){
            #Aggregate counts by summing
            gcdf <- as.data.frame(assays(ExpObj)[[CurrAssay]])
            gcdf$Feats <- rownames(gcdf)
            gcdf <- left_join(gcdf, feats2glomby_feats, by = "Feats")
            gcdf$Feats <- NULL
            gcdf <- aggregate(. ~ Glomby_feats, data = gcdf, FUN = sum)
            rownames(gcdf) <- gcdf$Glomby_feats
            gcdf$Glomby_feats <- NULL
            gcdf <- gcdf[featureorder, sampleorder]
            assays[[CurrAssay]] <- gcdf
        }
    }

    #Iron out the feature table
    #Determine which version of JAMS the SEobj was made by, ensure backwards compatibility.
    if (JAMS2TaxonomicObj){
        #Get a novel feature table
        glom_ftt <- ftt
        if (glomby %in% taxonomic_levels){
            taxonomic_levels_have <- taxonomic_levels[taxonomic_levels %in% colnames(glom_ftt)]
            taxonomic_levels_to_keep <- taxonomic_levels_have[1:which(taxonomic_levels_have == glomby)]
            glom_ftt <- glom_ftt[ , taxonomic_levels_to_keep]

            #Again, account for completely dark matter (i.e. LKT__Unclassified) not having a taxonomic path.
            if ("LKT__Unclassified" %in% rownames(glom_ftt)){
                #Attribute LKT__Unclassified to that glomby
                glom_ftt[which(rownames(glom_ftt) == "LKT__Unclassified"), taxonomic_levels_to_keep] <- "LKT__Unclassified"
            }
            glom_ftt <- glom_ftt[ , taxonomic_levels_to_keep]
            #Try de-duplicating the entire taxonomy table data frame. Depending on how the taxonomy was built, there may be more than one alternate taxonomic path leading to a node. All was done to avoid this with JAMSbuildk2db, but you never know.
            glom_ftt <- glom_ftt[!duplicated(glom_ftt), ]
            #Check whether featureorder is the same length as glom_ftt rows and if they are 100% contained within
            if ((nrow(glom_ftt) == length(featureorder)) & (all(featureorder %in% glom_ftt[ , glomby]))){
                rownames(glom_ftt) <- glom_ftt[ , glomby]
                glom_ftt <- glom_ftt[featureorder, ]
            } else {
                #There may be some alternate taxonomic paths leading to the same terminal (glomby) leaf. de-duplicate the glomby and keep the first in each case.
                flog.warn(paste("Alternate taxonomic paths leading to the same node at", glomby, "detected. Keeping the first path in each instance."))
                glom_ftt <- glom_ftt[!duplicated(glom_ftt[ , glomby]), ]
                rownames(glom_ftt) <- glom_ftt[ , glomby]
                glom_ftt <- glom_ftt[featureorder, ]
            }
        }
    } else {
        #treat as before
        #Get classes for novel features
        glomby_feats <- unique(ftt[ , glomby])
        #Get a novel feature table
        glom_ftt <- ftt[(!duplicated(ftt[, glomby])), 1:(which(colnames(ftt) == glomby)), drop = FALSE]
        rownames(glom_ftt) <- glom_ftt[, glomby]
    }

    #Rebuild the SummarizedExperiment object
    assays <- assays[sapply(1:length(assays), function (x) { !is.null(assays[[x]]) })]

    glomExpObj <- SummarizedExperiment(assays = assays, rowData = glom_ftt, colData = pheno_original)
    metadata(glomExpObj) <- metadata(ExpObj)

    return(glomExpObj)
}
