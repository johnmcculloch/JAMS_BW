#' cluster_strains(genes_df = NULL, normalize_length = FALSE, distmethod = "bray", cutoff = NULL)
#'
#' Clusters MAGs of a single working taxid into functional clades by Bray-Curtis
#' dissimilarity of their per-Product gene-length profiles, hierarchical average-linkage
#' clustering, and a dissimilarity cut. Used internally by contextualize_taxonomy (q.v.).
#'
#' Hoisted to package top-level (was a closure inside contextualize_taxonomy) so that when
#' this is called from a parallel worker, R does NOT serialize the enclosing call frame
#' (which held the whole list.data) along with it. It now receives everything it needs as
#' explicit arguments.
#'
#' @param genes_df Long data frame with columns Accession (the Product/functional accession),
#'   MAG_Accession, and LengthDNA (and, if normalize_length = TRUE, ProportionLengthDNA).
#' @param normalize_length Logical. If TRUE, cluster on ProportionLengthDNA rather than raw LengthDNA.
#' @param distmethod Dissimilarity index passed to vegan::vegdist. Default "bray".
#' @param cutoff Height at which to cut the dendrogram (stats::cutree h=).
#' @export

cluster_strains <- function(genes_df = NULL, normalize_length = FALSE, distmethod = "bray", cutoff = NULL){
    #Aggregate lengths of features
    if (normalize_length){
        length_sum_df <- genes_df %>% group_by(Accession, MAG_Accession) %>% summarise(Total_LengthDNA = sum(ProportionLengthDNA, na.rm = TRUE), .groups = "drop")
    } else {
        length_sum_df <- genes_df %>% group_by(Accession, MAG_Accession) %>% summarise(Total_LengthDNA = sum(LengthDNA, na.rm = TRUE), .groups = "drop")
    }
    length_sum_df <- length_sum_df %>% tidyr::pivot_wider(names_from = MAG_Accession, values_from = Total_LengthDNA, values_fill = 0)
    length_sum_df <- as.data.frame(length_sum_df)
    rownames(length_sum_df) <- length_sum_df$Accession
    length_sum_df$Accession <- NULL

    #Compute pairwise distance between MAGs (rows after transpose)
    d_length_sum <- vegdist(t(as.matrix(length_sum_df)), method = distmethod)

    #Cluster and bin
    curr_hc <- stats::hclust(d_length_sum, method = "average")
    hc_entity_clusters <- stats::cutree(curr_hc, h = cutoff)

    #Add annotation
    curr_strain_df <- as.data.frame(hc_entity_clusters)
    colnames(curr_strain_df)[1] <- "Cluster_Number"
    curr_strain_df$MAG_Accession <- rownames(curr_strain_df)

    return(curr_strain_df)
}
