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
#' ## Feature-set consistency across subsets (IMPORTANT)
#' The individual plotting functions, when called standalone with a subsetby, first establish
#' a GLOBAL set of surviving features across ALL samples (using the user's filtration
#' parameters), and then hold that set fixed while subsetting, passing c(0,0) filters within
#' each subset. This prevents the surviving-feature count from dwindling (and differing) from
#' subset to subset, since featcutoff/GenomeCompletenessCutoff are prevalence-based and the
#' sample count changes per subset.
#'
#' Because this comparator subsets OUTSIDE the plotting functions (handing each a per-subset
#' samplesToKeep with subsetby = NULL), that built-in global-feature machinery never engages.
#' To preserve the same guarantee, this function now computes a single global feature set
#' (global_FTK) ONCE, up front, across the retained samples, using the user's filtration
#' parameters. It then passes featuresToKeep = global_FTK to every plotting call with all
#' filtration arguments set to non-filtering values (applyfilters = NULL, featcutoff = c(0,0),
#' GenomeCompletenessCutoff = c(0,0)). Every panel in every subset therefore considers exactly
#' the same features, and comparisons across subsets are honest.
#'
#' The global feature set is computed on an object that has ALREADY been vetted in the same
#' order the plotting functions use (ExpObjVetting applies samplesToKeep/featuresToKeep, then
#' only_allow_CSBs, then glomby). This matters because CSB restriction must precede
#' agglomeration and agglomeration rewrites feature names: computing global_FTK in the final
#' (post-CSB, post-glom) namespace guarantees the names match what each plotting function sees.
#'
#' @param ExpObj A single JAMS-style SummarizedExperiment object (i.e. expvec[[analysis]]).
#' @param compareby String; the metadata variable defining the groups being compared.
#' @param subsetby String or vector (up to 3) of metadata variables for tiered subsetting. If NULL, a single "no subset" pass is made over all samples.
#' @param samplesToKeep Optional vector of sample names to retain BEFORE any global filtration or subsetting. Use this to drop samples you know you do not want in the comparison at all (e.g. outliers). If NULL (default), all samples are kept. Applied first, ahead of everything else.
#' @param featuresToKeep Optional vector of feature names to retain BEFORE any global filtration or subsetting. Use this to restrict to (or drop) features you care about up front (e.g. remove noise, or focus on a shortlist). If NULL (default), all features are considered. Applied first, ahead of the global filtration that derives global_FTK. For taxonomic objects that will be agglomerated via glomby, these must be names in the ORIGINAL (pre-glom) namespace, since culling happens before agglomeration; see ExpObjVetting.
#' @param glomby String; taxonomic level to agglomerate to (taxonomic spaces only). Default NULL.
#' @param only_allow_CSBs Logical; restrict to Consolidated Species Bins (taxonomic ConsolidatedGenomeBin space only). Applied before any agglomeration. Default FALSE.
#' @param do_ordination,do_heatmap,do_alpha Logical switches to include/exclude each of the three plot types. All default TRUE.
#' @param ordination_highlight_subset_in_context Logical; if set to TRUE, ordination will be done with all samples within the SummarizedExperiment object, but only the subset will be highlighted. PERMANOVA and centroid calculations will be performed on the highlighted samples ONLY.
#' @param ordination_args,heatmap_args,alpha_args Named lists of extra arguments passed to plot_Ordination, plot_relabund_heatmap and plot_alpha_diversity respectively, overriding the comparator's shared defaults for that function only.
#' @param split_heatmap_columns_by_compareby Logical; if set to TRUE, heatmap columns will be split into the same samples split by compareby
#' @param applyfilters,featcutoff,GenomeCompletenessCutoff,PPM_normalize_to_bases_sequenced,cdict,class_to_ignore Shared arguments. The filtration arguments (applyfilters, featcutoff, GenomeCompletenessCutoff) are used ONCE here to derive global_FTK and are NOT forwarded to the plotting functions (which receive non-filtering values instead). PPM_normalize_to_bases_sequenced, cdict and class_to_ignore are forwarded to all three functions unless overridden per-function via the *_args lists.
#' @param return_stats Logical; if TRUE, invisibly returns a named list collecting the returnstats output of the heatmap and alpha functions per subset. Default TRUE.
#'
#' @export

compare_within_subsets <- function(ExpObj = NULL, compareby = NULL, subsetby = NULL,
                                   samplesToKeep = NULL, featuresToKeep = NULL,
                                   glomby = NULL, only_allow_CSBs = FALSE,
                                   do_ordination = TRUE, ordination_highlight_subset_in_context = FALSE, do_heatmap = TRUE, split_heatmap_columns_by_compareby = TRUE, do_alpha = TRUE,
                                   ordination_args = list(), heatmap_args = list(), alpha_args = list(),
                                   applyfilters = NULL, featcutoff = NULL,
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
    #n.b. CSBs only exist in ConsolidatedGenomeBin space
    if (only_allow_CSBs && !(analysis %in% "ConsolidatedGenomeBin")){
        flog.warn(paste0("only_allow_CSBs = TRUE was passed, but analysis space \"", analysis, "\" is not taxonomic. The plotting functions will ignore it."))
    }
    if (!is.null(glomby) && !(analysis %in% taxonomic_spaces)){
        flog.warn(paste0("glomby = \"", glomby, "\" was passed, but analysis space \"", analysis, "\" is not taxonomic. Agglomeration may be ignored or error; consider glomby = NULL for functional spaces."))
    }

    #############################################################################
    ## STEP 0: Up-front user culling (samplesToKeep / featuresToKeep).
    ## These are applied FIRST, before any global filtration or subsetting, so a
    ## user can knowingly drop outlier samples or noise/uninteresting features
    ## before anything else happens. We do this via ExpObjVetting so that the
    ## ordering invariants (samples/features cull -> CSB restriction -> glom) are
    ## owned by the same vetting routine the plotting functions use, keeping the
    ## feature namespace consistent.
    ##
    ## NOTE on glomby here: we deliberately DO NOT agglomerate at this step. We
    ## only resolve the retained sample set and the global feature set. Each
    ## plotting function will perform its own CSB restriction and agglomeration
    ## (we forward only_allow_CSBs and glomby to them). However, global_FTK must
    ## be expressed in the SAME namespace the plotting functions will filter in,
    ## which is the post-CSB, post-glom namespace. So to DERIVE global_FTK we vet
    ## a throwaway copy WITH glomby applied, read off its surviving feature names,
    ## and hand those back as featuresToKeep. filter_experiment/ExpObjVetting in
    ## each plotting function accept post-glom feature names in featuresToKeep
    ## (see their documentation), so this lines up correctly.
    #############################################################################

    #Resolve the retained sample universe up front (used for both global_FTK and subsetting).
    all_samples <- rownames(colData(ExpObj))
    if (!is.null(samplesToKeep)){
        samplesToKeep <- unique(samplesToKeep)
        retained_samples <- all_samples[all_samples %in% samplesToKeep]
        dropped_n <- length(all_samples) - length(retained_samples)
        if (dropped_n > 0){
            flog.info(paste0("samplesToKeep: retaining ", length(retained_samples), " of ", length(all_samples), " samples (dropping ", dropped_n, ") before any filtration or subsetting."))
        }
        if (length(retained_samples) < 2){
            stop("Fewer than 2 samples remain after applying samplesToKeep. Nothing to compare.")
        }
    } else {
        retained_samples <- all_samples
    }

    #Report any explicitly requested features we will try to honour.
    if (!is.null(featuresToKeep)){
        featuresToKeep <- unique(featuresToKeep)
        present_ftk <- featuresToKeep[featuresToKeep %in% rownames(ExpObj)]
        if (length(present_ftk) < length(featuresToKeep)){
            flog.warn(paste0("featuresToKeep: ", (length(featuresToKeep) - length(present_ftk)), " of ", length(featuresToKeep), " requested feature(s) were not found in the (pre-glom) object and will be ignored."))
        }
        if (length(present_ftk) < 1){
            flog.warn("featuresToKeep: none of the requested features were found. Proceeding as if featuresToKeep were NULL.")
            featuresToKeep <- NULL
        } else {
            featuresToKeep <- present_ftk
        }
    }

    #############################################################################
    ## STEP 1: Derive the single GLOBAL feature set (global_FTK).
    ## Computed ONCE, across the retained samples, in the final (post-CSB,
    ## post-glom) namespace, using the user's filtration parameters. This mirrors
    ## the global_FTK logic inside plot_relabund_heatmap / plot_Ordination, which
    ## we are otherwise bypassing by subsetting outside those functions.
    #############################################################################

    #Translate applyfilters/featcutoff/GenomeCompletenessCutoff into a concrete preset once.
    presetlist <- declare_filtering_presets(analysis = analysis, applyfilters = applyfilters, featcutoff = featcutoff, GenomeCompletenessCutoff = GenomeCompletenessCutoff)

    #Vet a copy in the same order the plotting functions use: user cull -> CSB -> glom.
    #This yields an object whose feature names are in the final namespace.
    vetted_for_FTK <- ExpObjVetting(ExpObj = ExpObj, samplesToKeep = retained_samples, featuresToKeep = featuresToKeep, glomby = glomby, only_allow_CSBs = only_allow_CSBs, variables_to_fix = NULL, class_to_ignore = NULL)

    #Decide whether any filtration is actually requested. If not, the global feature
    #set is simply every feature in the vetted object (still a fixed, shared set).
    any_filtration <- any(c(all(presetlist$featcutoff != c(0, 0)),
                            (!is.null(presetlist$GenomeCompletenessCutoff) && all(presetlist$GenomeCompletenessCutoff != c(0, 0))),
                            !is.null(applyfilters)))

    if (any_filtration){
        #Apply the user's filters ONCE across all retained samples to get the global set.
        global_obj <- filter_experiment(SEobj = vetted_for_FTK,
                                        featcutoff = presetlist$featcutoff,
                                        samplesToKeep = NULL,
                                        featuresToKeep = NULL,
                                        only_allow_CSBs = FALSE,
                                        GenomeCompletenessCutoff = presetlist$GenomeCompletenessCutoff,
                                        normalization = "relabund",
                                        PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced,
                                        flush_out_empty_samples = FALSE,
                                        give_info = TRUE)
        global_FTK <- rownames(global_obj)
        global_obj <- NULL
        flog.info(paste0("Established a global feature set of ", length(global_FTK), " feature(s) across all ", length(retained_samples), " retained samples. This set is held FIXED across every subset and every plot type."))
    } else {
        global_FTK <- rownames(vetted_for_FTK)
        flog.info(paste0("No filtration requested; the global feature set is all ", length(global_FTK), " feature(s) in the (vetted) object. Held fixed across every subset and plot type."))
    }
    vetted_for_FTK <- NULL

    if (length(global_FTK) < 2){
        stop("Fewer than 2 features survive the global filtration. Relax your filtration parameters (applyfilters/featcutoff/GenomeCompletenessCutoff) or check featuresToKeep.")
    }

    #From here on, the plotting functions must NOT re-filter. They receive the fixed
    #global feature set and non-filtering parameters.
    FIXED_applyfilters <- NULL
    FIXED_featcutoff <- c(0, 0)
    FIXED_GenomeCompletenessCutoff <- c(0, 0)

    #############################################################################
    ## STEP 2: Resolve subset iteration order and names, scoped to retained_samples.
    ## We still do not pre-filter features here; each plotting function receives
    ## featuresToKeep = global_FTK and a per-subset samplesToKeep.
    #############################################################################
    if (!is.null(subsetby)){
        #Subset only within the retained sample universe. We hand multiple_subsetting
        #a phenotable restricted to retained_samples so dropped samples never reappear.
        retained_pheno <- as.data.frame(colData(ExpObj))[retained_samples, , drop = FALSE]
        subset_list <- multiple_subsetting_sample_selector(SEobj = ExpObj, phenotable = retained_pheno, subsetby = subsetby, compareby = compareby, cats_to_ignore = class_to_ignore)
        if (is.null(subset_list)){
            flog.warn("Subsetting could not be resolved (check that subsetby variables are discrete). Falling back to a single no-subset pass over the retained samples.")
            subset_list <- NULL
            subset_names <- "no_sub"
        } else {
            subset_df <- subset_list$Subsets_stats
            subset_df <- subset_df[which(subset_df$Subset_Tier_Level != 0), , drop = FALSE]
            #Drop subsets too small to compare at all.
            if (any(subset_df$Num_samples_in_subset < 2)){
                LowSampSubsets <- subset_df[which(subset_df$Num_samples_in_subset < 2), "Subset_Tier_Class_Name"]
                flog.warn(paste("Subsets", paste0(LowSampSubsets, collapse = ", "), "contain fewer than 2 samples and will be skipped entirely."))
                subset_df <- subset_df[which(subset_df$Num_samples_in_subset >= 2), , drop = FALSE]
            }
            if (nrow(subset_df) < 1){
                flog.warn("No surviving subsets. Falling back to a single no-subset pass over the retained samples.")
                subset_names <- "no_sub"
            } else {
                subset_names <- subset_df$Subset_Tier_Class_Name
            }
        }
    } else {
        subset_list <- NULL
        subset_names <- "no_sub"
    }

    #Helper to resolve which samples belong to the current subset (always a subset of retained_samples).
    get_subset_samples <- function(subname){
        if (subname == "no_sub" || is.null(subset_list)){
            return(retained_samples)
        } else {
            smp <- subset_list[[subname]]
            return(smp[smp %in% retained_samples])
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
        subtit <- if (subname == "no_sub") paste("All retained samples |", analysis) else paste(subname, "|", analysis)
        tryCatch({
            plot.new()
            grid.table(c(paste("SUBSET:", subtit),
                         paste("compareby:", compareby),
                         paste("n samples:", length(STK)),
                         paste("n features (global, fixed):", length(global_FTK))),
                       rows = NULL, cols = NULL,
                       theme = ttheme_default(base_size = 14))
        }, error = function(e){ flog.warn(paste("Could not draw banner for", subname, ":", conditionMessage(e))) })

        flog.info(paste("==== Comparator: subset", subname, "----", length(STK), "samples ===="))

        ##############################
        ## 1. Ordination (is it different overall?)
        ##############################
        if (do_ordination){
            if (ordination_highlight_subset_in_context){
                #Context mode: plot against the full retained backdrop, highlight the subset.
                ordSTK <- retained_samples
                ordSTH <- STK
            } else {
                ordSTK <- STK
                ordSTH <- NULL
            }

            ord_shared <- list(ExpObj = ExpObj, samplesToKeep = ordSTK, samplesToHighlight = ordSTH, subsetby = NULL,
                               glomby = glomby, only_allow_CSBs = only_allow_CSBs,
                               featuresToKeep = global_FTK,
                               compareby = compareby, colourby = compareby,
                               applyfilters = FIXED_applyfilters, featcutoff = FIXED_featcutoff,
                               GenomeCompletenessCutoff = FIXED_GenomeCompletenessCutoff,
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
                             featuresToKeep = global_FTK,
                             hmtype = "comparative", compareby = compareby, splitcolsby = splitcolsby,
                             applyfilters = FIXED_applyfilters, featcutoff = FIXED_featcutoff,
                             GenomeCompletenessCutoff = FIXED_GenomeCompletenessCutoff,
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
                             featuresToKeep = global_FTK,
                             compareby = compareby, colourby = compareby, fillby = NULL,
                             applyfilters = FIXED_applyfilters, featcutoff = FIXED_featcutoff,
                             GenomeCompletenessCutoff = FIXED_GenomeCompletenessCutoff,
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