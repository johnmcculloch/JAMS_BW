#' retrieve_features_by_taxa(FuncExpObj = NULL, glomby = NULL, only_allow_CSBs = FALSE, LKTsToKeep = NULL, assay_for_matrix = "BaseCounts", wantedfeatures = NULL, wantedsamples = NULL, asPPM = TRUE, PPM_normalize_to_bases_sequenced = TRUE, PPMthreshold = 0, append_metatada = TRUE, include_samples_with_zero = TRUE, return_taxonomy_table = FALSE)
#'
#' Returns a long form data frame of stratification by taxa of the relative abundance or number of bases wanted of functional features in wanted samples, given allfeaturesbytaxa_matrix and allfeaturesbytaxa_index metadata present in a JAMS SummarizedExperiment functional object.
#' @export

retrieve_features_by_taxa <- function(FuncExpObj = NULL, glomby = NULL, only_allow_CSBs = FALSE, LKTsToKeep = NULL, assay_for_matrix = "BaseCounts", wantedfeatures = NULL, wantedsamples = NULL, asPPM = TRUE, PPM_normalize_to_bases_sequenced = TRUE, PPMthreshold = 0, append_metatada = TRUE, include_samples_with_zero = TRUE, return_taxonomy_table = FALSE){

    if (assay_for_matrix == "GeneCounts"){
        sparsematrix_name <- "allfeaturesbytaxa_GeneCounts_matrix"
        sparsematrix_index_name <- "allfeaturesbytaxa_GeneCounts_index"
        asPPM <- FALSE
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
    rowsinterestdf$RowNumber <- as.character(rowsinterestdf$RowNumber)
    #Eliminate unavailable information.
    rowsinterestdf <- rowsinterestdf[!is.na(rowsinterestdf$RowNumber), , drop = FALSE]
    #Safer to go by named rows than numeric, as it will throw an error if the name does not exist.
    allfeaturesbytaxa_interest <- as.matrix(metadata(FuncExpObj)[[sparsematrix_name]][as.character(rowsinterestdf$RowNumber), , drop = FALSE])
    #Prune empty taxa, i.e. taxa which have 0 counts for all of the features of interest
    allfeaturesbytaxa_interest <- as.data.frame(allfeaturesbytaxa_interest[, which(colSums(allfeaturesbytaxa_interest) != 0), drop = FALSE])
    allfeaturesbytaxa_interest$RowNumber <- as.character(rownames(allfeaturesbytaxa_interest))
    allfeaturesbytaxa_interest <- left_join(allfeaturesbytaxa_interest, rowsinterestdf, by = "RowNumber")
    allfeaturesbytaxa_interest$RowNumber <- NULL
    allfeaturesbytaxa_interest <- allfeaturesbytaxa_interest[ , c("Sample", "Accession", (sort(colnames(allfeaturesbytaxa_interest)[which(!colnames(allfeaturesbytaxa_interest) %in% c("Sample", "Accession"))])))]
    rownames(allfeaturesbytaxa_interest) <- paste(allfeaturesbytaxa_interest$Sample, allfeaturesbytaxa_interest$Accession, sep = "§")

    if (include_samples_with_zero){
        SampleAccession_tally <- expand.grid(list(Sample = wantedsamples, Accession = unique(allfeaturesbytaxa_interest$Accession)), KEEP.OUT.ATTRS = TRUE, stringsAsFactors = FALSE)
        SampleAccession_tally$SampleAccession <- paste(SampleAccession_tally$Sample, SampleAccession_tally$Accession, sep = "§")
        rownames(SampleAccession_tally) <- SampleAccession_tally$SampleAccession
        missingSApairs <- SampleAccession_tally$SampleAccession[!SampleAccession_tally$SampleAccession %in% rownames(allfeaturesbytaxa_interest)]
        #Only append if we're actually missing anything
        if (length(missingSApairs) > 0){
            suppl_df <- SampleAccession_tally[missingSApairs, c("Sample", "Accession")]
            Taxoncols <- colnames(allfeaturesbytaxa_interest)[!(colnames(allfeaturesbytaxa_interest) %in% c("Sample", "Accession"))]
            suppl_df[ , Taxoncols] <- 0
            suppl_df <- suppl_df[ , colnames(allfeaturesbytaxa_interest)]
            allfeaturesbytaxa_interest <- rbind(allfeaturesbytaxa_interest, suppl_df)
        }
    }

    Taxoncols <- colnames(allfeaturesbytaxa_interest)[!(colnames(allfeaturesbytaxa_interest) %in% c("Sample", "Accession"))]

    if (only_allow_CSBs){
        CSBsToKeep <- Taxoncols[grep("^CSB_", Taxoncols)]
        #Stop if there is nothing left. You will be amazed as to what people use as inputs.
        if (length(CSBsToKeep) < 1){
            flog.warn("There are less than 1 feature in which are valid CSBs in the input SummarizedExperiment object's taxonomic stratification sparse matrix.")
            stop("There is no information left to stratify with. Check your filtration parameters or set only_allow_CSBs - FALSE")
        } else {
            flog.warn(paste("Filtering out features which are not CSBs (Consolidated Species Bins), i.e. features which were not binned by MetaBAT2 in JAMSalpha."))
            Taxoncols <- CSBsToKeep
            #Subset taxa to what was requested.
            allfeaturesbytaxa_interest <- allfeaturesbytaxa_interest[ , c("Sample", "Accession", Taxoncols)]
        }  
    }

    #Restrict to LKTs before glomming, if required from LKTsToKeep
    if (!is.null(LKTsToKeep)){
        Taxoncols <- Taxoncols[Taxoncols %in% LKTsToKeep]
        #Check there is at least one left
        if (length(Taxoncols) < 1){
            flog.warn("No LKTs passed to LKTsToKeep were found within the taxonomic stratification matrix. Check your inputs.")
            stop("Aborting")
        }
        #Subset taxa to what was requested.
        allfeaturesbytaxa_interest <- allfeaturesbytaxa_interest[ , c("Sample", "Accession", Taxoncols)]
    }

    #Agglomerate if applicable
    if ((!is.null(glomby) && (glomby != "LKT"))){
        #Check this feature is available before proceeding
        if (!JAMS2_object){
            stop("Stratification of functions by taxonomic levels higher than LKT is only available for SummarizedExperiment objects made with JAMS version > 2.1.0 ")
        }
        #Obtain taxonomy table
        stt <- as.data.frame(metadata(FuncExpObj)$tt)
        LKTs_to_glom <- (sort(colnames(allfeaturesbytaxa_interest)[which(!colnames(allfeaturesbytaxa_interest) %in% c("Sample", "Accession", "Ultra_low_abundance_LKTs"))]))
        stt <- stt[which(stt$LKT %in% LKTs_to_glom), c(glomby, "LKT")]
        #Fill in unclassifieds
        stt[which(is.na(stt[ , glomby])), glomby] <- "Unclassified"
        stt$Taxid <- sapply(stt$LKT, function (x) { extract_NCBI_taxid_from_featname(Taxon = x) } )

        #Fix level unclassifieds (as opposed to complete Unclassifieds)
        taxonomic_levels <- c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "IS1", "LKT")
        taxonomic_levels_tags <- c("d__", "k__", "p__", "c__", "o__", "f__", "g__", "s__", "is1__", "LKT__")
        #Deal with unclassifieds in Kingdom to Species
        if (glomby %in% c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")){
            #Reallocate LKT to tag__Unclassified to avoid glomming different LKTs at that level into the same glom
            appropriate_tag <- taxonomic_levels_tags[which(taxonomic_levels == glomby)]
            stt[which(stt[ , glomby] == paste0(appropriate_tag, "Unclassified")), glomby] <- stt[which(stt[ , glomby] == paste0(appropriate_tag, "Unclassified")), "LKT"]
        }

        #extract matrix to glom
        LKTs_to_glom_mat <- allfeaturesbytaxa_interest[ , c("Sample", "Accession", LKTs_to_glom)]
        #Filter LKTs present in stt
        common_LKTs <- intersect(colnames(LKTs_to_glom_mat), stt$LKT)

        LKTs_filtered <- LKTs_to_glom_mat[ , c("Sample", "Accession", common_LKTs), drop = FALSE]
        LKTs_filtered$SampleAccession <- paste(LKTs_filtered$Sample, LKTs_filtered$Accession, sep = "-")
        SAsafeLookup <- LKTs_filtered[ , c("SampleAccession", "Sample", "Accession")]
        LKTs_filtered$Sample <- NULL
        LKTs_filtered$Accession <- NULL

        LKTs_filtered <- LKTs_filtered[ , c("SampleAccession", common_LKTs), drop = FALSE]
        # Reshape LKTs_to_glom_mat to long format
        LKTs_to_glom_long <- LKTs_filtered %>% pivot_longer(cols = -SampleAccession, names_to = "LKT", values_to = "NumBases")
        # Join with stt to get glomby level for each LKT
        LKTs_to_glom_long <- left_join(LKTs_to_glom_long, stt, by = "LKT")
        LKTs_to_glom_long <- as.data.frame(LKTs_to_glom_long)

        feats2glomby_feats <- LKTs_to_glom_long[ , c("LKT", glomby, "Taxid")]
        feats2glomby_feats <- feats2glomby_feats[!duplicated(feats2glomby_feats$LKT), ]
        colnames(feats2glomby_feats)[which(colnames(feats2glomby_feats) == glomby)] <- "Glomby_feats"
        colnames(feats2glomby_feats)[which(colnames(feats2glomby_feats) == "LKT")] <- "Feats"

        feats2glomby_feats <- fix_feats2glomby_feats(feats2glomby_feats = feats2glomby_feats, glomby = glomby)

        #Remove stale glomby information
        LKTs_to_glom_long[ , glomby] <- NULL

        #Replace with the correct info from feats2glomby_feats
        colnames(feats2glomby_feats)[which(colnames(feats2glomby_feats) == "Feats")] <- "LKT"
        LKTs_to_glom_long <- left_join(LKTs_to_glom_long, feats2glomby_feats, by = "LKT")
        #Rename relevant glomby column to make things more practical
        colnames(LKTs_to_glom_long)[which(colnames(LKTs_to_glom_long) == "Glomby_feats")] <- "Glomby"

        #Aggregate counts per row_id × Glomby
        Glomby_sum <- LKTs_to_glom_long %>% group_by(SampleAccession, Glomby) %>% summarise(NumBases = sum(NumBases, na.rm = TRUE), .groups = "drop")
        #Spread back to wide format (one row per original row)
        Glomby_wide <- Glomby_sum %>% pivot_wider(names_from = Glomby, values_from = NumBases, values_fill = 0)
        Glomby_wide <- as.data.frame(Glomby_wide)
        Glomby_wide <- left_join(SAsafeLookup, Glomby_wide, by = "SampleAccession")
        rownames(Glomby_wide) <- Glomby_wide$SampleAccession
        Glomby_wide$SampleAccession <- NULL
        allfeaturesbytaxa_interest <- Glomby_wide

        #Iron out taxonomy table if requested
        if (return_taxonomy_table == TRUE){
            featureorder <- colnames(Glomby_wide)[!(colnames(Glomby_wide) %in% c("Sample", "Accession"))]
            ftt <- as.data.frame(metadata(FuncExpObj)$tt)
            glom_ftt <- ftt

            rownames(glom_ftt) <- 1:nrow(glom_ftt)
            rownames(feats2glomby_feats) <- 1:nrow(feats2glomby_feats)

            #Mark current features as leaves
            colnames(feats2glomby_feats)[which(colnames(feats2glomby_feats) == "Feats")] <- "LKT"
            #Mark Glomby_feats as the current agglomerating level
            colnames(feats2glomby_feats)[which(colnames(feats2glomby_feats) == "Glomby_feats")] <- glomby

            glom_ftt[ , glomby] <- NULL
            glom_ftt$Taxid <- NULL
            glom_ftt <- left_join(glom_ftt, feats2glomby_feats, by = "LKT")
            #Eliminate unused rows
            glom_ftt <- glom_ftt[!is.na(glom_ftt[ , glomby]), ]

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
            taxtable <- glom_ftt
        }

    } else {

        if (return_taxonomy_table == TRUE){
            taxtable <- as.data.frame(metadata(FuncExpObj)$tt)
        }

    }

    #Normalize if applicable
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
        TaxonColumns <- colnames(taxsplit)[!(colnames(taxsplit) %in% c("Sample", "Accession", "NumBases"))]

        #Transform to PPM. Keep 2 significant decimal places because splitting a 1PPM gene several ways may lead to fractions.
        for (colm in TaxonColumns){
            taxsplit[ , colm] <- round(((taxsplit[ , colm] / taxsplit$NumBases) * 1000000), 2)
        }

        #Denoise
        TaxaMaxima <- sapply(TaxonColumns, function(x) {max(taxsplit[ , x])})
        TaxaToKeep <- names(which(TaxaMaxima > PPMthreshold))
        if (append_metatada){
            sample2metadata <- as.data.frame(colData(FuncExpObj))
            nonSamplecolms <- colnames(sample2metadata)[colnames(sample2metadata) != "Sample"]
            taxsplit <- left_join(taxsplit, sample2metadata, by = "Sample")
        } else {
            nonSamplecolms <- NULL
        }
        taxsplit <- taxsplit[ , c("Sample", "Accession", "NumBases", nonSamplecolms, TaxaToKeep)]
        taxsplit$NumBases <- NULL

    } else {

        taxsplit <- allfeaturesbytaxa_interest

    }

    if (return_taxonomy_table == TRUE){
        return(list(taxsplit = taxsplit, taxtable = taxtable))
    } else {
        return(taxsplit)
    }
}
