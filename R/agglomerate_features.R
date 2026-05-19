#' agglomerate_features(ExpObj = NULL, glomby = NULL)
#'
#' Agglomerates features in a JAMS SummarizedExperiment object safely
#' @export

agglomerate_features <- function(ExpObj = NULL, glomby = NULL){

    #Get appropriate object to work with
    if (as.character(class(ExpObj)[1]) != "SummarizedExperiment"){
        stop("This function can only take a SummarizedExperiment object as input.")
    }

    terminal_taxonomic_levels <- c("Contig_LKT", "MB2bin", "ConsolidatedGenomeBin", "LKT", "is1_")
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
    #Add taxid if taxonomic
    if (JAMS2TaxonomicObj){
        ftt$Taxid <- sapply(ftt$LKT, function (x) { extract_NCBI_taxid_from_featname(Taxon = x) } )
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
    feats2glomby_feats <- data.frame(Feats = rownames(ftt), Glomby_feats = as.character(ftt[ , glomby]), Taxid = ftt[ , "Taxid"], stringsAsFactors = FALSE)
    #Account for completely dark matter (i.e. LKT__Unclassified) not having a taxonomic path.
    if ("LKT__Unclassified" %in% feats2glomby_feats$Feats){
        #Attribute LKT__Unclassified to that glomby
        feats2glomby_feats[which(feats2glomby_feats$Feats == "LKT__Unclassified"), "Glomby_feats"] <- "LKT__Unclassified"
    }
    rownames(feats2glomby_feats) <- feats2glomby_feats$Feats

    #Update from JAMS_BW version 2.2.2 onwards: aggregate by taxid, rather than name. This is preferred, because taxid on names of unclassified species are not lower than species level. So, CSB__g__12345_Bacterium_FC1 and CSB__g__12345_Bacterium_FC2 will have the same taxid. FCs are tags for infraspecies level, so, there should be no FCs when agglomerating by Species or higher level.
    #Split feats2glomby_feats into those glommed at the desired level and those still remaining CSBs or LKTs
    CSB_tags <- paste("^CSB", taxonomic_levels_tags[2:(which(taxonomic_levels == glomby) - 1)], sep = "__")
    LKT_tags <- paste("^LKT", taxonomic_levels_tags[2:(which(taxonomic_levels == glomby) - 1)], sep = "__")
    general_expressions <- paste0(c(CSB_tags, LKT_tags), collapse = "|")

    ### look for CSBs or LKTs which are NOT below species. i.e. LKT__s__ and LKT__is1__ are safe.
    feats2glomby_feats_to_resolve <- feats2glomby_feats[grep(general_expressions, feats2glomby_feats$Glomby_feats), ]

    #Remove those from feats2glomby_feats
    feats2glomby_feats <- feats2glomby_feats[setdiff(rownames(feats2glomby_feats), rownames(feats2glomby_feats_to_resolve)), ]
    #Fix glommed names
    #use JAMStaxtable
    data(JAMStaxtable)
    feats2glomby_feats_to_resolve <- left_join(feats2glomby_feats_to_resolve, JAMStaxtable[ , c("Taxid", "LKT")], by = "Taxid") 
    #Commit the LKT annotation to Glomby_feats
    feats2glomby_feats_to_resolve$Glomby_feats <- feats2glomby_feats_to_resolve$LKT
    feats2glomby_feats_to_resolve$LKT <- NULL

    #Reconstruct feats2glomby_feats
    feats2glomby_feats <- rbind(feats2glomby_feats, feats2glomby_feats_to_resolve)
    feats2glomby_feats$Taxid <- NULL

    cts <- left_join(cts, feats2glomby_feats, by = "Feats")
    cts$Feats <- NULL

    #glom_cts <- aggregate(. ~ Glomby_feats, data = cts, FUN = sum, na.rm = TRUE)
    glom_cts <- rowsum(cts[, colnames(cts) != "Glomby_feats"], group = cts$Glomby_feats)
    #Rownames are set automatically after aggregation with rowsum
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
            #Rownames are set automatically after aggregation with rowsum
            gcdf <- rowsum(gcdf[, colnames(gcdf) != "Glomby_feats"], group = gcdf$Glomby_feats)
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

            rownames(glom_ftt) <- 1:nrow(glom_ftt)
            rownames(feats2glomby_feats) <- 1:nrow(feats2glomby_feats)

            #Mark current features as leaves
            colnames(feats2glomby_feats)[which(colnames(feats2glomby_feats) == "Feats")] <- "LKT"
            #Mark Glomby_feats as the current agglomerating level
            colnames(feats2glomby_feats)[which(colnames(feats2glomby_feats) == "Glomby_feats")] <- glomby

            #Purge the stale information column from glom_ftt
            glom_ftt[ , glomby] <- NULL
            glom_ftt$Taxid <- NULL
            glom_ftt <- left_join(glom_ftt, feats2glomby_feats, by = "LKT")

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
