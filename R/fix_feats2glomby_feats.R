#' fix_feats2glomby_feats(feats2glomby_feats = NULL, glomby = NULL)
#'
#' Provides the correct algorithm for agglomerating taxonomic objects. This function is used by agglomerate_features and retrieve_features_by_taxa. No not attempt to use it outside of this context.
#' @export

fix_feats2glomby_feats <- function(feats2glomby_feats = NULL, glomby = NULL){

    feats2glomby_feats <- as.data.frame(feats2glomby_feats)

    #Account for completely dark matter (i.e. LKT__Unclassified) not having a taxonomic path.
    if ("LKT__Unclassified" %in% feats2glomby_feats[ , "Feats"]){
        #Attribute LKT__Unclassified to that glomby
        feats2glomby_feats[which(feats2glomby_feats[ , "Feats"] == "LKT__Unclassified"), "Glomby_feats"] <- "LKT__Unclassified"
    }
    rownames(feats2glomby_feats) <- feats2glomby_feats[ , "Feats"]

    #Update from JAMS_BW version 2.2.2 onwards: aggregate by taxid, rather than name. This is preferred, because taxid on names of unclassified species are not lower than species level. So, CSB__g__12345_Bacterium_FC1 and CSB__g__12345_Bacterium_FC2 will have the same taxid. FCs are tags for infraspecies level, so, there should be no FCs when agglomerating by Species or higher level.
    #Split feats2glomby_feats into those glommed at the desired level and those still remaining CSBs or LKTs
    taxonomic_levels <- c("Gram", "Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "IS1", "LKT")
    taxonomic_levels_tags <- c("", "d__", "k__", "p__", "c__", "o__", "f__", "g__", "s__", "is1__", "LKT__")

    CSB_tags <- paste("^CSB", taxonomic_levels_tags[2:(which(taxonomic_levels == glomby) - 1)], sep = "__")
    LKT_tags <- paste("^LKT", taxonomic_levels_tags[2:(which(taxonomic_levels == glomby) - 1)], sep = "__")
    general_expressions <- paste0(c(CSB_tags, LKT_tags), collapse = "|")

    ### look for CSBs or LKTs which are NOT below species. i.e. LKT__s__ and LKT__is1__ are safe.
    feats2glomby_feats_to_resolve <- feats2glomby_feats[grep(general_expressions, feats2glomby_feats[ , "Glomby_feats"]), ]

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

    return(feats2glomby_feats)

}
