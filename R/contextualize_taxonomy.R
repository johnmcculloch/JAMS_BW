#' contextualize_taxonomy(LKTdosesall = LKTdosesall, list.data = list.data, normalize_length = FALSE, dissimilarity_cutoff = 0.15, threads = 1, totmembytes = NULL)
#'
#' Internal function used exclusively within make_SummarizedExperiments to cluster taxonomic
#' entities belonging to the same species taxid into functional clades. Do not attempt to use
#' this out of this context.
#'
#' MEMORY REFACTOR (paths 1+2): the heavy per-sample featuredata is now pre-reduced ONCE on the
#' master process (build_compact_mag_profiles) into a small per-MAG x per-Product profile, and
#' the parallel clustering step operates only on tiny per-species slices of that profile. The
#' clustering helpers (cluster_strains, rename_MAG_Accession) are package top-level functions,
#' so the %dopar% closures no longer capture list.data. Previously, PSOCK workers each received
#' a full copy of list.data (tens of GB x N workers), which is what produced the ~1.3 TB peak
#' and the multi-hour, near-idle-CPU serialization stall.
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
    ## Exclude bins with NO taxonomic anchor from functional contextualization.
    ## A high-quality MetaBAT2 bin whose k-mer classification resolved to root (taxid 1),
    ## cellular organisms (taxid 131567) or unclassified (taxid 0) is named e.g.
    ## "MB2__12__Unclassified" and carries no d__/k__/p__/.../s__/is1__ tag. Several such
    ## bins from different samples share the same placeholder taxid, so they would otherwise be
    ## grouped into a single "working taxid" and clustered together by function, which is
    ## biologically meaningless (they are unrelated dark matter) and, worse, crashes
    ## rename_MAG_Accession downstream ("missing value where TRUE/FALSE needed", because no
    ## taxonomic tag can be found in the name). These bins keep their original LKT.
    #############################################################################
    has_taxonomic_tag <- function(x){
        split_elements <- unlist(strsplit(x, split = "_"))
        any(split_elements %in% c("d", "k", "p", "c", "o", "f", "g", "s", "is1"))
    }
    untagged <- (is.na(BinsDF$Taxid)) | (as.character(BinsDF$Taxid) %in% c("0", "1", "131567")) | (!sapply(BinsDF$ConsolidatedGenomeBin, has_taxonomic_tag))
    if (any(untagged)){
        flog.warn(paste0("contextualize_taxonomy: ", sum(untagged), " HQ/MHQ bin(s) have no taxonomic anchor (unclassified at every rank, e.g. root/cellular organisms) and will NOT be functionally contextualized. They keep their original name. Bins: ", paste0(BinsDF$MAG_Accession[untagged], collapse = ", ")))
        BinsDF <- BinsDF[which(!untagged), , drop = FALSE]
        if (nrow(BinsDF) == 0){
            return(LKTdosesall)
        }
    }

    #############################################################################
    ## Self-audit: reconcile abundance-table CGBs against per-sample featuredata.
    ## This is cheap and runs regardless of thread count, on the master. It tells us, directly
    ## in the log, whether a CGB with sequencing depth but no annotated features is present in
    ## this dataset (the condition that previously caused "replacement has 1 row, data has 0").
    ## If everything reconciles, a single reassuring line is logged.
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
                flog.warn(paste0("contextualize_taxonomy self-audit: ", length(offenders), " of ", n_checked, " HQ/MHQ ConsolidatedGenomeBins have sequencing depth in the abundance table but NO matching annotated features in the sample featuredata. These bins will be skipped during functional contextualization. Offending bins follow:"))
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
    audit_featuredata_reconciliation(BinsDF = BinsDF, verbose = TRUE)

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
    #Only cluster working taxids which have > 1 MAG in it. Otherwise, it defeats the purpose.
    WorkingTaxids_to_decon <- names(which(table(BinsDF$WorkingTaxid) > 1))

    #Collect working taxids whose usable-MAG count drops below 2 after the featuredata cull,
    #so their surviving bin is routed through the singleton path rather than silently lost.
    demoted_to_singleton_MAGs <- character(0)

    #Small, data-free lookup for renaming infraspecies MAGs to their species label. Tiny, safe
    #to ship to workers (replaces passing the whole BinsDF).
    mag2species <- setNames(BinsDF$Species, BinsDF$MAG_Accession)

    if (length(WorkingTaxids_to_decon) > 0){

        #########################################################################
        ## PATH 1: pre-reduce featuredata ONCE on the master into a compact profile.
        ## Only MAGs belonging to multi-MAG working taxids are profiled (singletons never need
        ## a gene profile). From here on, list.data is not read again, and in particular it is
        ## NEVER handed to a parallel worker.
        #########################################################################
        MAGs_to_profile <- BinsDF[which(BinsDF$WorkingTaxid %in% WorkingTaxids_to_decon), c("Sample", "ConsolidatedGenomeBin", "MAG_Accession"), drop = FALSE]

        flog.info(paste0("contextualize_taxonomy: pre-reducing featuredata for ", nrow(MAGs_to_profile), " MAG(s) across ", length(unique(MAGs_to_profile$Sample)), " sample(s) into a compact per-MAG functional profile (single master pass)."))

        profile_out <- build_compact_mag_profiles(MAGs_to_profile = MAGs_to_profile, list.data = list.data, normalize_length = normalize_length)
        compact <- profile_out$profiles

        #Log any skip reasons gathered during pre-reduction (master process).
        if (length(profile_out$skip_reasons) > 0){
            flog.warn(paste0("contextualize_taxonomy: ", length(profile_out$skip_reasons), " bin/sample(s) were skipped while building functional profiles (present in abundance tables, unusable in featuredata). First few:"))
            max_report <- 50
            to_report <- profile_out$skip_reasons[1:min(length(profile_out$skip_reasons), max_report)]
            for (msg in to_report){
                flog.warn(paste0("   ", msg))
            }
            if (length(profile_out$skip_reasons) > max_report){
                flog.warn(paste0("   ... and ", (length(profile_out$skip_reasons) - max_report), " more not shown."))
            }
        }

        if (!is.null(compact) && nrow(compact) > 0){

            #Attach working taxid and split into per-species slices.
            mag2wt <- setNames(BinsDF$WorkingTaxid, BinsDF$MAG_Accession)
            compact$WorkingTaxid <- mag2wt[compact$MAG_Accession]
            profile_by_wt <- split(compact, compact$WorkingTaxid)

            #Build a list of slices aligned to WorkingTaxids_to_decon (NULL where every MAG was culled).
            profile_slices <- lapply(WorkingTaxids_to_decon, function(wt) profile_by_wt[[as.character(wt)]])

            #Instrumentation: prove the exported footprint is small (compare to the old ~tens of GB).
            biggest_slice_bytes <- suppressWarnings(max(sapply(profile_slices, function(x) if (is.null(x)) 0 else as.numeric(object.size(x)))))
            flog.info(sprintf("contextualize_taxonomy: compact profile is %.1f MB total; largest single-species slice shipped to a worker is %.1f MB; MAG->Species lookup is %.2f MB.",
                as.numeric(object.size(compact)) / 1e6, biggest_slice_bytes / 1e6, as.numeric(object.size(mag2species)) / 1e6))

            #Decide thread count. We reuse calculate_safe_threads but pass the objects that are
            #ACTUALLY exported now (the compact profile and the tiny lookup), not list.data. This
            #is deliberately conservative (uses the full compact size as the per-worker estimate),
            #but because the compact profile is small it will normally permit the full thread count.
            #NOTE (future increment): make calculate_safe_threads profile-aware rather than reusing it like this.
            safe_memory_threads <- calculate_safe_threads(requested_threads = threads, totmembytes = totmembytes, list.data = compact, BinsDF = mag2species)

            num_tasks <- length(WorkingTaxids_to_decon)
            optimal_threads <- min(safe_memory_threads, num_tasks)
            if (optimal_threads < threads){
                flog.info(sprintf("Throttling threads from %d to %d (bounded by number of tasks and compact-profile memory, not by list.data anymore).", threads, optimal_threads))
            }

            #####################################################################
            ## PATH 2: cluster each species' slice. Workers call the top-level helpers and
            ## receive only their slice + the tiny mag2species lookup. list.data, LKTdosesall,
            ## BinsDF, the full compact table and profile_by_wt are explicitly NOT exported.
            #####################################################################
            if (threads > 1){

                flog.info(paste("Deconvoluting", num_tasks, "taxa using", optimal_threads, "parallel threads..."))
                cl <- parallel::makeCluster(optimal_threads)
                doParallel::registerDoParallel(cl)

                #n.b. mag2species, dissimilarity_cutoff and normalize_length are captured
                #automatically by foreach (they are referenced in the body), so they must NOT be
                #listed in .export as well, or foreach warns "already exporting variable(s)".
                worker_out <- foreach::foreach(slice = profile_slices, WT = WorkingTaxids_to_decon,
                    .packages = c("dplyr", "tidyr", "vegan"),
                    .export = c("cluster_strains", "rename_MAG_Accession"),
                    .noexport = c("list.data", "LKTdosesall", "BinsDF", "compact", "profile_by_wt", "profile_slices")) %dopar% {

                    #Valid MAGs for this taxon are exactly those present in its slice.
                    valid_in_wt <- if (is.null(slice)) character(0) else unique(slice$MAG_Accession)

                    if (length(valid_in_wt) < 2){
                        #Nothing to cluster; hand back any lone survivor to be treated as a singleton.
                        return(list(strain_df = NULL, demoted = valid_in_wt, WorkingTaxid = WT, error = NULL))
                    }

                    #Never let a single taxon kill the whole JAMSbeta run. On failure, report
                    #exactly which working taxid and MAGs were involved and fall back to treating
                    #the MAGs as singletons (they keep a species-level CSB name, FC1).
                    res <- tryCatch({
                        curr_strain_df <- cluster_strains(genes_df = slice, normalize_length = normalize_length, cutoff = dissimilarity_cutoff)
                        curr_strain_df$ContextualizedSpecies <- sapply(1:nrow(curr_strain_df), function(x) {
                            rename_MAG_Accession(MAG_Accession = curr_strain_df[x, "MAG_Accession"], Cluster_Number = curr_strain_df[x, "Cluster_Number"], mag2species = mag2species)
                        })
                        list(strain_df = curr_strain_df, demoted = character(0), WorkingTaxid = WT, error = NULL)
                    }, error = function(e){
                        list(strain_df = NULL, demoted = valid_in_wt, WorkingTaxid = WT, error = conditionMessage(e))
                    })
                    return(res)
                }

                parallel::stopCluster(cl)

                #Harvest results and gather demoted singletons.
                for (wo in worker_out){
                    if (!is.null(wo$error)){
                        flog.warn(paste0("contextualize_taxonomy: functional clustering FAILED for working taxid ", wo$WorkingTaxid, " with error \"", wo$error, "\". The ", length(wo$demoted), " MAG(s) in this taxon will be contextualized as singletons instead. MAGs: ", paste0(wo$demoted, collapse = ", ")))
                    }
                    if (length(wo$demoted) > 0){
                        demoted_to_singleton_MAGs <- c(demoted_to_singleton_MAGs, wo$demoted)
                    }
                    if (!is.null(wo$strain_df)){
                        strain_df <- rbind(strain_df, wo$strain_df)
                    }
                }

            } else {
                #Single-threaded path, same logic, same compact slices.
                flog.info(paste("Deconvoluting", num_tasks, "taxa using a single thread..."))

                for (ti in seq_along(WorkingTaxids_to_decon)){
                    slice <- profile_slices[[ti]]
                    valid_in_wt <- if (is.null(slice)) character(0) else unique(slice$MAG_Accession)

                    if (length(valid_in_wt) < 2){
                        demoted_to_singleton_MAGs <- c(demoted_to_singleton_MAGs, valid_in_wt)
                        next
                    }

                    curr_strain_df <- tryCatch({
                        csd <- cluster_strains(genes_df = slice, normalize_length = normalize_length, cutoff = dissimilarity_cutoff)
                        csd$ContextualizedSpecies <- sapply(1:nrow(csd), function(x) {
                            rename_MAG_Accession(MAG_Accession = csd[x, "MAG_Accession"], Cluster_Number = csd[x, "Cluster_Number"], mag2species = mag2species)
                        })
                        csd
                    }, error = function(e){
                        flog.warn(paste0("contextualize_taxonomy: functional clustering FAILED for working taxid ", WorkingTaxids_to_decon[ti], " with error \"", conditionMessage(e), "\". The ", length(valid_in_wt), " MAG(s) in this taxon will be contextualized as singletons instead. MAGs: ", paste0(valid_in_wt, collapse = ", ")))
                        NULL
                    })
                    if (is.null(curr_strain_df)){
                        demoted_to_singleton_MAGs <- c(demoted_to_singleton_MAGs, valid_in_wt)
                        next
                    }
                    strain_df <- rbind(strain_df, curr_strain_df)
                }
            }

        } else {
            flog.warn("contextualize_taxonomy: no usable functional profiles could be built for any multi-MAG taxon. No clustering will be performed; affected MAGs keep their original LKT.")
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
        strain_df_singletons$ContextualizedSpecies <- sapply(1:nrow(strain_df_singletons), function (x) { rename_MAG_Accession(MAG_Accession = strain_df_singletons[x, "MAG_Accession"], Cluster_Number = strain_df_singletons[x, "Cluster_Number"], mag2species = mag2species) } )
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