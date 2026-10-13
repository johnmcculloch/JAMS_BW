#' rename_MAG_Accession(MAG_Accession = NULL, Cluster_Number = NULL, mag2species = NULL)
#'
#' Builds a Contextualized/Consolidated Species Bin (CSB) name for a MAG given its functional
#' cluster number. Used internally by contextualize_taxonomy (q.v.).
#'
#' Hoisted to package top-level (see cluster_strains for the reason). It previously took the
#' entire BinsDF to look up a Species label for infraspecies-level MAGs; it now takes a tiny
#' named lookup vector (mag2species) instead, so nothing large travels to parallel workers.
#'
#' @param MAG_Accession The "Sample§ConsolidatedGenomeBin" key for the MAG.
#' @param Cluster_Number The functional cluster number assigned by cluster_strains.
#' @param mag2species Named character vector mapping MAG_Accession -> Species label. Only
#'   consulted for infraspecies-level (is1) MAGs, which are promoted to their Species.
#' @export

rename_MAG_Accession <- function(MAG_Accession = NULL, Cluster_Number = NULL, mag2species = NULL){

    split_elements <- unlist(strsplit(MAG_Accession, split = "_"))
    #Get earliest position matching a taxonomic tag
    taxtag_pos <- which(split_elements %in% c("d", "k", "p", "c", "o", "f", "g", "s", "is1"))[1]

    #Keep current taxid, unless infraspecies (is1), in which case look up the species label.
    if (split_elements[taxtag_pos] != "is1"){
        #rejoin everything downstream from that
        CSB <- paste(split_elements[taxtag_pos:length(split_elements)], collapse = "_")
    } else {
        #MAG_Accession is at the infraspecies level, so look up the species label.
        CSB <- mag2species[[MAG_Accession]]
    }

    #Tack on the Functional Cluster k-Number (FuCk)
    CSB <- paste(CSB, paste0("FC", Cluster_Number), sep = "_")
    #Add Contextualized Species Bin (CSB) tag
    CSB <- paste("CSB", CSB, sep = "__")

    return(CSB)
}
