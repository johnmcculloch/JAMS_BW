#' contextualize_taxonomy(LKTdosesall = LKTdosesall, list.data = list.data, normalize_length = FALSE, dissimilarity_cutoff = 0.15, threads = 1, totmembytes = NULL)
#'
#' This is an internal function used exclusively within the make_SummarizedExperiments function to cluster taxonomic entities belonging to the same species taxid into functional clades. Do not attempt to use this out of this context.
#' @export

contextualize_taxonomy <- function(LKTdosesall = LKTdosesall, list.data = list.data, normalize_length = FALSE, dissimilarity_cutoff = 0.15, threads = 1, totmembytes = NULL){

    if (threads > 1){
        #Assure necessary packages are loaded
        library(foreach)
        library(doParallel)
    }

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

    #Helper: build the per-MAG gene data frame for a single working taxid, dropping any MAG
    #whose ConsolidatedGenomeBin has depth in the abundance table but NO annotated features in
    #the corresponding sample's featuredata (which would otherwise throw
    #"replacement has 1 row, data has 0" when columns are assigned into a zero-row data frame),
    #or whose features sum to zero total length (which would give NaN proportions).
    #Returns a list with the bound data frame (or NULL) and a character vector of skip reasons,
    #so that the parallel workers can hand reasons back to the master process for logging.
    build_curr_genes_df <- function(SampleEntities_df = NULL){
        skip_reasons <- character(0)
        curr_genes_df_list <- lapply(1:nrow(SampleEntities_df), function(rn) {
            curr_sample <- SampleEntities_df$Sample[rn]
            curr_CGB <- SampleEntities_df$ConsolidatedGenomeBin[rn]
            curr_MAG_Accession <- SampleEntities_df$MAG_Accession[rn]
            curr_featuredata <- list.data[[paste(curr_sample, "featuredata", sep = "_")]]

            #Defensive: the featuredata object itself may be missing.
            if (is.null(curr_featuredata)){
                skip_reasons <<- c(skip_reasons, paste0("Sample=", curr_sample, " | CGB=", curr_CGB, " | MAG=", curr_MAG_Accession, " -> featuredata object not found in list.data"))
                return(NULL)
            }

            curr_func_df <- subset(curr_featuredata, ConsolidatedGenomeBin == curr_CGB)[, c("Feature", "LengthDNA", "Product", "ConsolidatedGenomeBin")]

            #Skip bins that have depth in the abundance table but no annotated features in featuredata.
            if (nrow(curr_func_df) < 1){
                skip_reasons <<- c(skip_reasons, paste0("Sample=", curr_sample, " | CGB=", curr_CGB, " | MAG=", curr_MAG_Accession, " -> 0 features in featuredata (present in abundance table, absent in featuredata)"))
                return(NULL)
            }

            curr_func_df$LengthDNA <- as.numeric(curr_func_df$LengthDNA)
            colnames(curr_func_df)[which(colnames(curr_func_df) == "Product")] <- "Accession"

            #Guard against a zero (or NA) total length, which would give NaN proportions.
            total_len <- sum(curr_func_df$LengthDNA, na.rm = TRUE)
            if (is.na(total_len) || total_len <= 0){
                skip_reasons <<- c(skip_reasons, paste0("Sample=", curr_sample, " | CGB=", curr_CGB, " | MAG=", curr_MAG_Accession, " -> features present but total LengthDNA is 0 or NA"))
                return(NULL)
            }

            curr_func_df$ProportionLengthDNA <- curr_func_df$LengthDNA / total_len
            curr_func_df$Sample <- curr_sample
            curr_func_df$MAG_Accession <- curr_MAG_Accession

            return(curr_func_df)
        })

        #Drop bins with no usable feature data before binding.
        curr_genes_df_list <- curr_genes_df_list[!sapply(curr_genes_df_list, is.null)]
        if (length(curr_genes_df_list) < 1){
            curr_genes_df <- NULL
        } else {
            curr_genes_df <- do.call(rbind, curr_genes_df_list)
        }

        return(list(curr_genes_df = curr_genes_df, skip_reasons = skip_reasons))
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

    #############################################################################
    ## Self-audit: reconcile abundance-table CGBs against per-sample featuredata.
    ## This is cheap and runs regardless of thread count. It tells us, directly in
    ## the log, whether the mechanism behind the "replacement has 1 row, data has 0"
    ## crash (a CGB with sequencing depth but no annotated features) is present in
    ## this dataset. If everything reconciles, a single reassuring line is logged.
    #############################################################################
    audit_featuredata_reconciliation <- function(BinsDF = NULL, verbose = TRUE){
        n_checked <- 0
        offenders <- character(0)
        for (rn in 1:nrow(BinsDF)){
            curr_sample <- BinsDF$Sample[rn]
            curr_CGB <- BinsDF$ConsolidatedGenomeBin[rn]
            curr_featuredata <- list.data[[paste(curr_sample, "featuredata", sep = "_")]]
            n_checked <- n_checked + 1
            if (is.null(curr_featuredata)){
                offenders <- c(offenders, paste0("Sample=", curr_sample, " | CGB=", curr_CGB, " -> featuredata object missing"))
                next
            }
            n_feat <- sum(curr_featuredata$ConsolidatedGenomeBin == curr_CGB, na.rm = TRUE)
            if (n_feat == 0){
                offenders <- c(offenders, paste0("Sample=", curr_sample, " | CGB=", curr_CGB, " -> 0 features in featuredata"))
            }
        }

        if (verbose){
            if (length(offenders) == 0){
                flog.info(paste0("contextualize_taxonomy self-audit: all ", n_checked, " HQ/MHQ ConsolidatedGenomeBins reconcile with their sample featuredata (every bin with depth has annotated features). Taxonomic and functional domains are consistent."))
            } else {
                flog.warn(paste0("contextualize_taxonomy self-audit: ", length(offenders), " of ", n_checked, " HQ/MHQ ConsolidatedGenomeBins have sequencing depth in the abundance table but NO matching annotated features in the sample featuredata. These bins will be skipped during functional contextualization (this is the condition that previously caused the 'replacement has 1 row, data has 0' crash). Offending bins follow:"))
                #Cap the number reported so the log does not explode, but report how many were suppressed.
                max_report <- 50
                to_report <- offenders[1:min(length(offenders), max_report)]
                for (msg in to_report){
                    flog.warn(paste0("   ", msg))
                }
                if (length(offenders) > max_report){
                    flog.warn(paste0("   ... and ", (length(offenders) - max_report), " more not shown."))
                }
            }
        }

        return(invisible(NULL))
    }
    audit_featuredata_reconciliation(BinsDF = BinsDF, verbose = FALSE)

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

    #Collect working taxids which drop below 2 usable MAGs after the featuredata cull, so their
    #surviving bin can be routed through the singleton path rather than silently lost.
    demoted_to_singleton_MAGs <- character(0)

    if (length(WorkingTaxids_to_decon) > 0){

        # Calculate safe threads based on RAM
        safe_memory_threads <- calculate_safe_threads(
            requested_threads = threads, 
            totmembytes = totmembytes, 
            list.data = list.data, 
            BinsDF = BinsDF
        )

        # Don't use more threads than the number of tasks we actually have
        num_tasks <- length(WorkingTaxids_to_decon)
        optimal_threads <- min(safe_memory_threads, num_tasks)
        if (optimal_threads < threads) {
            flog.info(sprintf("Throttling threads from %d to %d to prevent memory exhaustion or overhead.", threads, optimal_threads))
        }

        # Determine whether deconvolution will proceed in single threaded or multi thread manner
        if (threads > 1) {

            flog.info(paste("Deconvoluting", num_tasks, "taxa using", optimal_threads, "parallel threads..."))
            #Setup cluster
            cl <- parallel::makeCluster(optimal_threads)
            doParallel::registerDoParallel(cl)

            #Each worker returns a list with the strain_df (or NULL), the skip reasons gathered
            #while building the gene data frame, and any bins demoted to singleton status. We
            #cannot reliably flog from inside workers, so we log everything from the master
            #process after the loop.
            worker_out <- foreach::foreach(WT = WorkingTaxids_to_decon, .packages = c('dplyr', 'tidyr', 'vegan')) %dopar% {
                SampleEntities_df <- subset(BinsDF, WorkingTaxid == WT)[, c("Sample", "ConsolidatedGenomeBin", "MAG_Accession")]
                built <- build_curr_genes_df(SampleEntities_df = SampleEntities_df)
                curr_genes_df <- built$curr_genes_df

                #If fewer than 2 MAGs survived the featuredata check, there is nothing to cluster.
                if (is.null(curr_genes_df) || length(unique(curr_genes_df$MAG_Accession)) < 2){
                    surviving_MAG <- if (!is.null(curr_genes_df)) unique(curr_genes_df$MAG_Accession) else character(0)
                    return(list(strain_df = NULL, skip_reasons = built$skip_reasons, demoted = surviving_MAG))
                }

                curr_strain_df <- cluster_strains(genes_df = curr_genes_df, normalize_length = normalize_length, cutoff = dissimilarity_cutoff)
                curr_strain_df$ContextualizedSpecies <- sapply(1:nrow(curr_strain_df), function(x) { 
                    rename_MAG_Accession(MAG_Accession = curr_strain_df[x, "MAG_Accession"], Cluster_Number = curr_strain_df[x, "Cluster_Number"], BinsDF = BinsDF) 
                })
                return(list(strain_df = curr_strain_df, skip_reasons = built$skip_reasons, demoted = character(0)))
            }
            
            # Stop cluster
            parallel::stopCluster(cl)

            #Harvest results, log skip reasons, and gather demoted singletons.
            for (wo in worker_out){
                if (length(wo$skip_reasons) > 0){
                    for (msg in wo$skip_reasons){
                        flog.warn(paste0("contextualize_taxonomy skipping bin during clustering: ", msg))
                    }
                }
                if (length(wo$demoted) > 0){
                    demoted_to_singleton_MAGs <- c(demoted_to_singleton_MAGs, wo$demoted)
                }
                if (!is.null(wo$strain_df)){
                    strain_df <- rbind(strain_df, wo$strain_df)
                }
            }

        } else {
            #Single threaded approach
            flog.info(paste("Deconvoluting", length(WorkingTaxids_to_decon), "taxa using a single thread..."))

            for (WT in WorkingTaxids_to_decon){
                SampleEntities_df <- subset(BinsDF, WorkingTaxid == WT)[, c("Sample", "ConsolidatedGenomeBin", "MAG_Accession")]
                built <- build_curr_genes_df(SampleEntities_df = SampleEntities_df)

                #Log skip reasons directly in the single-threaded case.
                if (length(built$skip_reasons) > 0){
                    for (msg in built$skip_reasons){
                        flog.warn(paste0("contextualize_taxonomy skipping bin during clustering: ", msg))
                    }
                }

                curr_genes_df <- built$curr_genes_df

                #If fewer than 2 MAGs survived the featuredata check, there is nothing to cluster.
                if (is.null(curr_genes_df) || length(unique(curr_genes_df$MAG_Accession)) < 2){
                    if (!is.null(curr_genes_df)){
                        demoted_to_singleton_MAGs <- c(demoted_to_singleton_MAGs, unique(curr_genes_df$MAG_Accession))
                    }
                    next
                }

                curr_strain_df <- cluster_strains(genes_df = curr_genes_df, normalize_length = normalize_length, cutoff = dissimilarity_cutoff)
                curr_strain_df$ContextualizedSpecies <- sapply(1:nrow(curr_strain_df), function(x) { 
                    rename_MAG_Accession(MAG_Accession = curr_strain_df[x, "MAG_Accession"], Cluster_Number = curr_strain_df[x, "Cluster_Number"], BinsDF = BinsDF) 
                })
                strain_df <- rbind(strain_df, curr_strain_df)
            }
        }
    }

    #Deal with singleton working taxids
    WorkingTaxids_singletons <- names(which(table(BinsDF$WorkingTaxid) == 1))
    singleton_MAGs <- BinsDF[which(BinsDF$WorkingTaxid %in% WorkingTaxids_singletons), "MAG_Accession"]

    #Add any MAGs demoted to singleton status (working taxid dropped below 2 usable MAGs after
    #the featuredata cull) so they are still contextualized rather than lost.
    if (length(demoted_to_singleton_MAGs) > 0){
        flog.warn(paste0("contextualize_taxonomy: ", length(demoted_to_singleton_MAGs), " bin(s) had their working taxid fall below 2 usable MAGs after the featuredata check and will be contextualized as singletons rather than clustered: ", paste0(demoted_to_singleton_MAGs, collapse = ", ")))
        singleton_MAGs <- unique(c(singleton_MAGs, demoted_to_singleton_MAGs))
    }

    if (length(singleton_MAGs) != 0){
        strain_df_singletons <- BinsDF[which(BinsDF$MAG_Accession %in% singleton_MAGs), "MAG_Accession", drop = FALSE]
        strain_df_singletons$Cluster_Number <- 1
        strain_df_singletons$ContextualizedSpecies <- sapply(1:nrow(strain_df_singletons), function (x) { rename_MAG_Accession(MAG_Accession = strain_df_singletons[x, "MAG_Accession"], Cluster_Number = strain_df_singletons[x, "Cluster_Number"], BinsDF = BinsDF) } )
        strain_df_singletons <- strain_df_singletons[ , c("Cluster_Number", "MAG_Accession", "ContextualizedSpecies")]
        strain_df <- rbind(strain_df, strain_df_singletons)
    }

    if (!is.null(strain_df) && nrow(strain_df) != 0){
        #Write changes to LKTdosesall
        #Replace the LKT value in LKTdosesall with the updated CGB
        LKTdosesall[rownames(strain_df), "LKT"] <- strain_df$ContextualizedSpecies
        rownames(LKTdosesall) <- 1:nrow(LKTdosesall)
    } else {
        flog.warn("contextualize_taxonomy: no entities were contextualized (strain_df is empty). Returning LKTdosesall unchanged.")
    }

    return(LKTdosesall)
}