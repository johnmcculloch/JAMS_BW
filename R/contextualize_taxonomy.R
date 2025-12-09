#' contextualize_taxonomy(LKTdosesall = LKTdosesall, list.data = list.data, normalize_length = FALSE, dissimilarity_cutoff = 0.15)
#'
#' This is an internal function used exclusively within the make_SummarizedExperiments function to cluster taxonomic entities belonging to the same species taxid into functional clades. Do not attempt to use this out of this context.
#' @export

contextualize_taxonomy <- function(LKTdosesall = LKTdosesall, list.data = list.data, normalize_length = FALSE, dissimilarity_cutoff = 0.15){

    data(JAMStaxtable)

    #n.b. There should be no duplicate ConsoildatedGenomeBin WITHIN a sample. It seems a few low level, bad quality ones may be duplicated in some samples. I don't know why, but meanwhile will aggregate these into a single CGB for the time being.
    LKTdosesall$MAG_Accession <- paste(LKTdosesall$Sample, LKTdosesall$ConsolidatedGenomeBin, sep = "§")
    if (any(duplicated(LKTdosesall$MAG_Accession))){
        dupes <- LKTdosesall$MAG_Accession[duplicated(LKTdosesall$MAG_Accession)]
        #extract dupes from LKTdosesall
        LKTdosesdupes <- LKTdosesall[which(LKTdosesall$MAG_Accession %in% dupes), ]
        LKTdosesall <- LKTdosesall[which(!(LKTdosesall$MAG_Accession %in% dupes)), ]
        LKTdosesdupes <- LKTdosesdupes %>% group_by(MAG_Accession) %>% summarise(Sample = Sample[1], ConsolidatedGenomeBin = ConsolidatedGenomeBin[1], Completeness = sum(Completeness), Contamination = sum(Contamination), Completeness_Model_Used = Completeness_Model_Used[1], NumBases = sum(NumBases), Taxid = Taxid[1], NCBI_taxonomic_rank = NCBI_taxonomic_rank[1], Domain = Domain[1], Kingdom = Kingdom[1], Phylum = Phylum[1], Class = Class[1], Order = Order[1], Family = Family[1], Genus = Genus[1], Species = Species[1], IS1 = IS1[1], LKT = LKT[1], Gram = Gram[1], RefScore = RefScore[1], PPM = sum(PPM), MAG_Accession = MAG_Accession[1])
        LKTdosesdupes <- as.data.frame(LKTdosesdupes)
        LKTdosesall <- rbind(LKTdosesall, LKTdosesdupes)
    }
    rownames(LKTdosesall) <- LKTdosesall$MAG_Accession

    #Define useful functions
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

        #Compute pairise distance
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

    rename_MAG_Accession <- function(MAG_Accession = NULL, Cluster_Number = NULL, BinsDF = NULL){

        split_elements <- unlist(strsplit(MAG_Accession, split = "_"))
        #Get earliest position matching a taxonomic tag
        taxtag_pos <- which(split_elements %in% c("d", "k", "p", "c", "o", "f", "g", "s", "is1"))[1]

        #Keep current taxid, unless infraspecies (is1), in which case look up the species taxid from BinsDF
        if (split_elements[taxtag_pos] != "is1"){
            #rejoin everything downstream from that
            CSB <- paste(split_elements[taxtag_pos:length(split_elements)], collapse = "_")
        } else {
            #MAG_Accession is at the infraspecies level, so look up the species level taxid from BinsDF
            CSB <- BinsDF[which(BinsDF$MAG_Accession == MAG_Accession), "Species"]
        }

        #Tack on the Functional Cluster k-Number (FuCk)
        CSB <- paste(CSB, paste0("FC", Cluster_Number), sep = "_")
        #Add Contextualized Species Bin (CSB) tag
        CSB <- paste("CSB", CSB, sep = "__")

        return(CSB)
    }

    #Rate all available MAG bins
    BinsDF <- rate_bin_quality(completeness_df = LKTdosesall, HQ_completeness_threshold = 90, HQ_contamination_threshold = 5, MHQ_completeness_threshold = 70, MHQ_contamination_threshold = 10, High_contamination_threshold = 20)

    #Consider only HQ and MHQ bins for taxonomic contextualization
    BinsDF <- subset(BinsDF, Quality %in% c("HQ", "MHQ"))

    #Return no change if nothing worthwile to contextualize was found.
    if (nrow(BinsDF) == 0){
        return(LKTdosesall)
    }

    BinsDF$MAG_Accession <- paste(BinsDF$Sample, BinsDF$ConsolidatedGenomeBin, sep = "§")

    #Set the Working Taxid for which entities will be clustered by function
    BinsDF$WorkingTaxid <- BinsDF$Taxid

    #Calculate Species taxid for taxonomic ranks which are downstream of Species
    if (any(BinsDF$NCBI_taxonomic_rank %in% c("strain", "subspecies"))){
        data(JAMStaxtable)
        downstream_taxids <- unique(BinsDF[which(BinsDF$NCBI_taxonomic_rank %in% c("strain", "subspecies")), "Taxid"])

        Taxid2WorkingTaxid <- data.frame(Taxid = downstream_taxids, WorkingTaxid = unname(sapply(JAMStaxtable[downstream_taxids, "Species"], function (x) { extract_NCBI_taxid_from_featname(Taxon = x) } )))
        for (tid in downstream_taxids){
            BinsDF[which(BinsDF$WorkingTaxid == tid), "WorkingTaxid"] <- Taxid2WorkingTaxid[which(Taxid2WorkingTaxid$Taxid == tid), "WorkingTaxid"]
        }
    }

    strain_df <- NULL
    #Only cluster or deconvolute working taxids which have > 1 MAG in it. Otherwise, it defeats the purpose.
    WorkingTaxids_to_decon <- names(which(table(BinsDF$WorkingTaxid) > 1))
    if (length(WorkingTaxids_to_decon) > 1){
        flog.info(paste("There are", length(WorkingTaxids_to_decon), "species-level taxonomical entities to deconvolute"))
        flog.info(paste("Entities will be clustered within bins presenting a distance >", dissimilarity_cutoff))

        #Loop round and resolve.
        for (WT in WorkingTaxids_to_decon){
            SampleEntities_df <- subset(BinsDF, WorkingTaxid == WT)[ , c("Sample", "ConsolidatedGenomeBin", "MAG_Accession")]
            #Build a data frame from which to obtain aggregate product (or functional feature) length sums.

            curr_genes_df <- NULL
            for (rn in 1:nrow(SampleEntities_df)){
                curr_func_df <- subset(list.data[[paste(SampleEntities_df$Sample[rn], "featuredata", sep = "_")]], ConsolidatedGenomeBin == SampleEntities_df$ConsolidatedGenomeBin[rn])[ , c("Feature", "LengthDNA", "Product", "ConsolidatedGenomeBin")]
                curr_func_df$LengthDNA <- as.numeric(curr_func_df$LengthDNA)
                colnames(curr_func_df)[which(colnames(curr_func_df) == "Product")] <- "Accession"
                curr_func_df$ProportionLengthDNA <- curr_func_df$LengthDNA / sum(curr_func_df$LengthDNA)
                curr_func_df$Sample <- SampleEntities_df$Sample[rn]
                curr_func_df$MAG_Accession <- SampleEntities_df$MAG_Accession[rn]
                curr_genes_df <- rbind(curr_genes_df, curr_func_df)
                curr_func_df <- NULL
            }

            curr_strain_df <- cluster_strains(genes_df = curr_genes_df, normalize_length = normalize_length, cutoff = dissimilarity_cutoff)
            curr_strain_df$ContextualizedSpecies <- sapply(1:nrow(curr_strain_df), function (x) { rename_MAG_Accession(MAG_Accession = curr_strain_df[x, "MAG_Accession"], Cluster_Number = curr_strain_df[x, "Cluster_Number"], BinsDF = BinsDF) } )
            strain_df <- rbind(strain_df, curr_strain_df)
            curr_strain_df <- NULL
        }
    }

    #Deal with singleton working taxids
    WorkingTaxids_singletons <- names(which(table(BinsDF$WorkingTaxid) == 1))
    if (length(WorkingTaxids_singletons) != 0){
        strain_df_singletons <- BinsDF[which(BinsDF$WorkingTaxid %in% WorkingTaxids_singletons), "MAG_Accession", drop = FALSE]
        strain_df_singletons$Cluster_Number <- 1
        strain_df_singletons$ContextualizedSpecies <- sapply(1:nrow(strain_df_singletons), function (x) { rename_MAG_Accession(MAG_Accession = strain_df_singletons[x, "MAG_Accession"], Cluster_Number = strain_df_singletons[x, "Cluster_Number"], BinsDF = BinsDF) } )
        strain_df_singletons <- strain_df_singletons[ , c("Cluster_Number", "MAG_Accession", "ContextualizedSpecies")]
        strain_df <- rbind(strain_df, strain_df_singletons)
    }

    if (nrow(strain_df) != 0){
        #Write changes to LKTdosesall
        #Replace the LKT value in LKTdosesall with the updated CGB
        LKTdosesall[rownames(strain_df), "LKT"] <- strain_df$ContextualizedSpecies
        rownames(LKTdosesall) <- 1:nrow(LKTdosesall)
    }

    return(LKTdosesall)
}