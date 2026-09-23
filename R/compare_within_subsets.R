#' compare_within_subsets(ExpObj = NULL, compareby = NULL, subsetby = NULL, ...)
#'
#' Orchestrates the three core JAMS comparison plots (ordination, relative abundance
#' heatmap, and alpha diversity) for a SINGLE analysis space, organised subset-by-subset.
#' For each subset tier produced by multiple_subsetting_sample_selector, the three plots
#' are emitted in the order ordination -> heatmap -> alpha diversity, answering, in turn:
#' "is the overall community different?", "which features drive/differ?", and "is richness
#' different?". This keeps one subset's complete story together before moving to the next.
#'
#' This function does NOT open a graphics device. Call it between your own pdf(...) and
#' dev.off(), because the heatmap (plot_relabund_heatmap) draws directly to the active
#' device and cannot be deferred into a list. Ordination and alpha diversity ggplot objects
#' are printed on the fly so that all output lands in your open device in subset order.
#'
#' The comparator delegates all object handling (CSB restriction, agglomeration, filtration)
#' to the individual plotting functions, passing samplesToKeep per subset with subsetby = NULL.
#' It never pre-vets the SummarizedExperiment itself, so ordering invariants (e.g. CSB
#' restriction must precede agglomeration) remain owned by the plotting functions.
#'
#' @param ExpObj A single JAMS-style SummarizedExperiment object (i.e. expvec[[analysis]]).
#' @param compareby String; the metadata variable defining the groups being compared.
#' @param subsetby String or vector (up to 3) of metadata variables for tiered subsetting. If NULL, a single "no subset" pass is made over all samples.
#' @param glomby String; taxonomic level to agglomerate to (taxonomic spaces only). Default NULL.
#' @param only_allow_CSBs Logical; restrict to Consolidated Species Bins (taxonomic ConsolidatedGenomeBin space only). Passed through to each plotting function, which apply it before any agglomeration. Default FALSE.
#' @param do_ordination,do_heatmap,do_alpha Logical switches to include/exclude each of the three plot types. All default TRUE.
#' @param ordination_highlight_subset_in_context Logical; if set to TRUE, ordination will be done with all samples within the SummarizedExperiment object, but only the subset will be highlighted. PERMANOVA and centroid calculations will be performed on the highlighted samples ONLY.
#' @param ordination_args,heatmap_args,alpha_args Named lists of extra arguments passed to plot_Ordination, plot_relabund_heatmap and plot_alpha_diversity respectively, overriding the comparator's shared defaults for that function only.
#' @param split_heatmap_columns_by_compareby Logical; if set to TRUE, heatmap columns will be split into the same samples split by compareby
#' @param applyfilters,featcutoff,GenomeCompletenessCutoff,PPM_normalize_to_bases_sequenced,cdict,class_to_ignore Shared arguments forwarded to all three functions unless overridden per-function via the *_args lists.
#' @param return_stats Logical; if TRUE, invisibly returns a named list collecting the returnstats output of the heatmap and alpha functions per subset. Default TRUE.
#'
#' @export

compare_within_subsets <- function(ExpObj = NULL, compareby = NULL, subsetby = NULL,
                                   glomby = NULL, only_allow_CSBs = FALSE,
                                   do_ordination = TRUE, ordination_highlight_subset_in_context = FALSE, do_heatmap = TRUE, split_heatmap_columns_by_compareby = TRUE, do_alpha = TRUE,
                                   ordination_args = list(), heatmap_args = list(), alpha_args = list(),
                                   applyfilters = "light", featcutoff = NULL,
                                   GenomeCompletenessCutoff = NULL,
                                   PPM_normalize_to_bases_sequenced = TRUE,
                                   cdict = NULL, class_to_ignore = "N_A",
                                   return_stats = TRUE){

    #Basic input checks
    if (as.character(class(ExpObj)[1]) != "SummarizedExperiment"){
        stop("ExpObj must be a single SummarizedExperiment object, e.g. expvec[[\"ConsolidatedGenomeBin\"]].")
    }
    if (is.null(compareby)){
        stop("You must supply a compareby variable.")
    }

    analysis <- metadata(ExpObj)$analysis
    flog.info(paste("Comparator running on analysis space:", analysis))

    #Warn on likely-misapplied taxonomic-only options
    taxonomic_spaces <- c("LKT", "Contig_LKT", "ConsolidatedGenomeBin", "MB2bin", "16S")
    if (only_allow_CSBs && !(analysis %in% taxonomic_spaces)){
        flog.warn(paste0("only_allow_CSBs = TRUE was passed, but analysis space \"", analysis, "\" is not taxonomic. The plotting functions will ignore it."))
    }
    if (!is.null(glomby) && !(analysis %in% taxonomic_spaces)){
        flog.warn(paste0("glomby = \"", glomby, "\" was passed, but analysis space \"", analysis, "\" is not taxonomic. Agglomeration may be ignored or error; consider glomby = NULL for functional spaces."))
    }

    #Resolve subset iteration order and names ONLY. We do not pre-filter the object;
    #each plotting function will vet itself when handed samplesToKeep per subset.
    if (!is.null(subsetby)){
        subset_list <- multiple_subsetting_sample_selector(SEobj = ExpObj, phenotable = NULL, subsetby = subsetby, compareby = compareby, cats_to_ignore = class_to_ignore)
        subset_df <- subset_list$Subsets_stats
        subset_df <- subset_df[which(subset_df$Subset_Tier_Level != 0), , drop = FALSE]
        #Drop subsets too small to compare at all.
        if (any(subset_df$Num_samples_in_subset < 2)){
            LowSampSubsets <- subset_df[which(subset_df$Num_samples_in_subset < 2), "Subset_Tier_Class_Name"]
            flog.warn(paste("Subsets", paste0(LowSampSubsets, collapse = ", "), "contain fewer than 2 samples and will be skipped entirely."))
            subset_df <- subset_df[which(subset_df$Num_samples_in_subset >= 2), , drop = FALSE]
        }
        if (nrow(subset_df) < 1){
            flog.warn("No surviving subsets. Falling back to a single no-subset pass over all samples.")
            subset_names <- "no_sub"
        } else {
            subset_names <- subset_df$Subset_Tier_Class_Name
        }
    } else {
        subset_list <- NULL
        subset_names <- "no_sub"
    }

    #Helper to resolve which samples belong to the current subset.
    get_subset_samples <- function(subname){
        if (subname == "no_sub" || is.null(subset_list)){
            return(rownames(colData(ExpObj)))
        } else {
            return(subset_list[[subname]])
        }
    }

    #Helper to merge shared defaults with per-function override lists, without
    #letting an override silently duplicate an argument (later wins).
    build_args <- function(shared, overrides){
        for (nm in names(overrides)){
            shared[[nm]] <- overrides[[nm]]
        }
        return(shared)
    }

    stats_collected <- list()

    for (subname in subset_names){

        STK <- get_subset_samples(subname)
        STK <- STK[STK %in% rownames(colData(ExpObj))]

        if (length(STK) < 2){
            flog.warn(paste("Subset", subname, "has fewer than 2 samples after resolution; skipping."))
            next
        }

        #A per-subset banner page, so the PDF is navigable subset-by-subset.
        subtit <- if (subname == "no_sub") paste("All samples |", analysis) else paste(subname, "|", analysis)
        tryCatch({
            plot.new()
            grid.table(c(paste("SUBSET:", subtit),
                         paste("compareby:", compareby),
                         paste("n samples:", length(STK))),
                       rows = NULL, cols = NULL,
                       theme = ttheme_default(base_size = 14))
        }, error = function(e){ flog.warn(paste("Could not draw banner for", subname, ":", conditionMessage(e))) })

        flog.info(paste("==== Comparator: subset", subname, "----", length(STK), "samples ===="))

        ##############################
        ## 1. Ordination (is it different overall?)
        ##############################
        if (do_ordination){
            if (ordination_highlight_subset_in_context){
                ordSTK <- NULL
                ordSTH <- STK
            } else {
                ordSTK <- STK
                ordSTH <- NULL
            }

            ord_shared <- list(ExpObj = ExpObj, samplesToKeep = ordSTK, samplesToHighlight = ordSTH, subsetby = NULL,
                               glomby = glomby, only_allow_CSBs = only_allow_CSBs,
                               compareby = compareby, colourby = compareby,
                               applyfilters = applyfilters, featcutoff = featcutoff,
                               GenomeCompletenessCutoff = GenomeCompletenessCutoff,
                               PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced,
                               cdict = cdict, class_to_ignore = class_to_ignore,
                               addtit = paste("SUBSET:", subname))
            ord_call <- build_args(ord_shared, ordination_args)
            tryCatch({
                ord_plots <- do.call(plot_Ordination, ord_call)
                #plot_Ordination returns a list of ggplots; print each into the open device.
                if (!is.null(ord_plots)){
                    for (gg in ord_plots){
                        print(gg)
                    }
                }
            }, error = function(e){
                flog.warn(paste("Ordination failed for subset", subname, ":", conditionMessage(e)))
            })
        }

        ##############################
        ## 2. Heatmap (which features drive/differ?)
        ##############################
        if (do_heatmap){
            if (split_heatmap_columns_by_compareby){
                splitcolsby <- compareby
            } else {
                splitcolsby <- NULL
            }

            hm_shared <- list(ExpObj = ExpObj, samplesToKeep = STK, subsetby = NULL,
                             glomby = glomby, only_allow_CSBs = only_allow_CSBs,
                             hmtype = "comparative", compareby = compareby, splitcolsby = splitcolsby,
                             applyfilters = applyfilters, featcutoff = featcutoff,
                             GenomeCompletenessCutoff = GenomeCompletenessCutoff,
                             PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced,
                             cdict = cdict, class_to_ignore = class_to_ignore,
                             addtit = paste("SUBSET:", subname),
                             returnstats = return_stats)
            hm_call <- build_args(hm_shared, heatmap_args)
            tryCatch({
                #plot_relabund_heatmap draws directly to the device; capture only stats.
                hm_stats <- plot_relabund_heatmap_safe(hm_call, return_stats = return_stats)
                if (return_stats && !is.null(hm_stats)){
                    stats_collected[[paste0("HM_", subname)]] <- hm_stats
                }
            }, error = function(e){
                flog.warn(paste("Heatmap failed for subset", subname, ":", conditionMessage(e)))
            })
        }

        ##############################
        ## 3. Alpha diversity (is richness different?)
        ##############################
        if (do_alpha){
            al_shared <- list(ExpObj = ExpObj, samplesToKeep = STK, subsetby = NULL,
                             glomby = glomby, only_allow_CSBs = only_allow_CSBs,
                             compareby = compareby, colourby = compareby, fillby = compareby,
                             applyfilters = applyfilters, featcutoff = featcutoff,
                             GenomeCompletenessCutoff = GenomeCompletenessCutoff,
                             PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced,
                             cdict = cdict, class_to_ignore = class_to_ignore,
                             addtit = paste("SUBSET:", subname),
                             returnstats = return_stats)
            al_call <- build_args(al_shared, alpha_args)
            tryCatch({
                al_out <- do.call(plot_alpha_diversity, al_call)
                #With returnstats = TRUE, gvec contains ggplots AND appended stats data frames.
                #ggplots print; data frames are collected. Distinguish by class.
                if (!is.null(al_out)){
                    for (nm in names(al_out)){
                        item <- al_out[[nm]]
                        if (inherits(item, "ggplot")){
                            print(item)
                        } else if (return_stats && is.data.frame(item)){
                            stats_collected[[paste0("ALPHA_", subname, "_", nm)]] <- item
                        }
                    }
                }
            }, error = function(e){
                flog.warn(paste("Alpha diversity failed for subset", subname, ":", conditionMessage(e)))
            })
        }

    } #End per-subset loop

    flog.info("Comparator complete.")
    if (return_stats){
        invisible(stats_collected)
    }
}


#' Internal helper: run plot_relabund_heatmap, which draws to the device as a side effect
#' and (optionally) returns its stats list. Isolated so the returnstats plumbing is explicit.
#' @keywords internal
plot_relabund_heatmap_safe <- function(arglist, return_stats = TRUE){
    out <- do.call(plot_relabund_heatmap, arglist)
    #When returnstats = TRUE, plot_relabund_heatmap returns svec (a list of stats data frames)
    #AND has already drawn to the device. When FALSE, it returns invisibly after drawing.
    if (return_stats){
        return(out)
    } else {
        return(NULL)
    }
}