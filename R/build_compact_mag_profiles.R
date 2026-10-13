#' build_compact_mag_profiles(MAGs_to_profile = NULL, list.data = NULL, normalize_length = FALSE)
#'
#' PRE-REDUCTION PASS (the memory fix). Makes ONE streaming pass over the samples referenced
#' in MAGs_to_profile, reading each sample's featuredata exactly once on the master process,
#' and emits a single COMPACT long table of per-MAG, per-Product summed gene length. This is a
#' tiny fraction of the full featuredata (HQ/MHQ multi-MAG taxa only; Product + length only;
#' pre-aggregated), and critically it is completely decoupled from list.data, so the parallel
#' clustering step never needs to serialize list.data to workers.
#'
#' This replaces the old build_curr_genes_df closure, which reached into list.data via lexical
#' scope from inside the %dopar% block and therefore caused the whole cohort's featuredata to
#' be copied to every worker.
#'
#' @param MAGs_to_profile Data frame with columns Sample, ConsolidatedGenomeBin, MAG_Accession
#'   (one row per MAG that will actually be clustered).
#' @param list.data The JAMSbeta list.data object (read here, on the master, exactly once per sample).
#' @param normalize_length Logical. If TRUE, also compute ProportionLengthDNA (per-MAG fraction).
#' @return A list with: profiles (compact long data frame, or NULL), valid_MAGs (character
#'   vector of MAGs with usable features), and skip_reasons (character vector for logging).
#' @export

build_compact_mag_profiles <- function(MAGs_to_profile = NULL, list.data = NULL, normalize_length = FALSE){

    samples <- unique(MAGs_to_profile$Sample)
    profile_list <- vector("list", length(samples))
    skip_reasons <- character(0)

    for (i in seq_along(samples)){
        s <- samples[i]
        fd <- list.data[[paste(s, "featuredata", sep = "_")]]
        curr_map <- MAGs_to_profile[which(MAGs_to_profile$Sample == s), , drop = FALSE]

        #Defensive: the featuredata object itself may be missing.
        if (is.null(fd)){
            skip_reasons <- c(skip_reasons, paste0("Sample=", s, " -> featuredata object not found in list.data (", nrow(curr_map), " MAG(s) affected)"))
            next
        }

        wanted_cgbs <- curr_map$ConsolidatedGenomeBin
        fd_sub <- fd[which(fd$ConsolidatedGenomeBin %in% wanted_cgbs), c("Feature", "LengthDNA", "Product", "ConsolidatedGenomeBin"), drop = FALSE]

        #Skip whole sample if none of its wanted bins have annotated features.
        if (nrow(fd_sub) < 1){
            skip_reasons <- c(skip_reasons, paste0("Sample=", s, " -> none of its ", length(wanted_cgbs), " wanted MAG(s) have annotated features in featuredata"))
            next
        }
        fd_sub$LengthDNA <- as.numeric(fd_sub$LengthDNA)

        #Aggregate summed gene length per CGB per Product (Product is the functional accession used for clustering).
        agg <- fd_sub %>% group_by(ConsolidatedGenomeBin, Product) %>% summarise(LengthDNA = sum(LengthDNA, na.rm = TRUE), .groups = "drop")
        agg <- as.data.frame(agg)
        colnames(agg)[which(colnames(agg) == "Product")] <- "Accession"

        #Attach MAG_Accession (and Sample, for traceability).
        cgb2mag <- setNames(curr_map$MAG_Accession, curr_map$ConsolidatedGenomeBin)
        agg$MAG_Accession <- cgb2mag[agg$ConsolidatedGenomeBin]
        agg$ConsolidatedGenomeBin <- NULL
        agg$Sample <- s

        profile_list[[i]] <- agg
    }

    #Drop empty samples before binding.
    profile_list <- profile_list[!sapply(profile_list, is.null)]
    if (length(profile_list) < 1){
        return(list(profiles = NULL, valid_MAGs = character(0), skip_reasons = skip_reasons))
    }
    compact <- do.call(rbind, profile_list)
    rownames(compact) <- NULL

    #Determine per-MAG total length; a MAG with zero/NA total is unusable (matches the old
    #per-MAG "0 features / 0 length" skip, which would otherwise route it to singleton).
    mag_totals <- compact %>% group_by(MAG_Accession) %>% summarise(Total = sum(LengthDNA, na.rm = TRUE), .groups = "drop")
    mag_totals <- as.data.frame(mag_totals)
    valid_MAGs <- mag_totals$MAG_Accession[which(mag_totals$Total > 0)]
    zero_MAGs <- mag_totals$MAG_Accession[which(!(mag_totals$Total > 0))]
    if (length(zero_MAGs) > 0){
        skip_reasons <- c(skip_reasons, paste0("MAG=", zero_MAGs, " -> features present but total LengthDNA is 0 or NA"))
    }

    #Also flag MAGs that were requested but produced no profile rows at all (bin absent from
    #featuredata even though the sample's featuredata existed) - informational only.
    missing_MAGs <- setdiff(unique(MAGs_to_profile$MAG_Accession), unique(compact$MAG_Accession))
    missing_MAGs <- setdiff(missing_MAGs, zero_MAGs)
    if (length(missing_MAGs) > 0){
        skip_reasons <- c(skip_reasons, paste0("MAG=", missing_MAGs, " -> no annotated features found for this bin in its sample featuredata"))
    }

    compact <- compact[which(compact$MAG_Accession %in% valid_MAGs), , drop = FALSE]

    if (normalize_length){
        totmap <- setNames(mag_totals$Total, mag_totals$MAG_Accession)
        compact$ProportionLengthDNA <- compact$LengthDNA / totmap[compact$MAG_Accession]
    }

    return(list(profiles = compact, valid_MAGs = valid_MAGs, skip_reasons = skip_reasons))
}
