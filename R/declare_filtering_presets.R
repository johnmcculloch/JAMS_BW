#' declare_filtering_presets(SEobj = NULL, analysis = NULL, is16S = NULL, glomby = NULL, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, PPM_normalize_to_bases_sequenced = FALSE, max_ranks_to_filter = 5, max_ranks_for_completeness = 3, fit_if_missing = TRUE, verbose = TRUE, ...)
#'
#' Resolves the JAMS filtration presets into the concrete cutoffs used by several functions. This
#' is the single place in JAMS2 where applyfilters is turned into numbers; everything downstream
#' receives the numbers, never the string.
#'
#' A filtration preset decides one thing only: whether a feature survives into the analysis at
#' all. It therefore sets featcutoff and, for taxonomic spaces, GenomeCompletenessCutoff, and
#' nothing else. Effect-size thresholds such as minl2fc, maxl2fc and minabscorrcoeff are NOT
#' filtration criteria. They decide how a surviving feature is judged or drawn, which is an
#' entirely separate question, and they are consequently set by their own arguments in the
#' functions which offer them, each with its own sensible default. Prior to JAMS 2.3.2 the
#' presets did set minl2fc and minabscorrcoeff, which coupled "is this feature real" to "is
#' this feature interesting" and made a plain applyfilters string silently change the
#' statistical verdict of a plot. That coupling has been removed.
#'
#' WHAT applyfilters CAN BE, as of JAMS 2.3.4. Every value is now FITTED to the object in hand.
#' The hard-coded c(50, 5) / c(250, 15) / c(2000, 15) constants are gone: they were chosen by eye
#' years ago on other people's data and referred to nothing about the cohort they were applied
#' to, which was the whole complaint. The four tiers are one ladder, walked further each step.
#'
#'   "denoise"    The gentlest fitted thresholds: the cell that removes the most noise per unit
#'                of real community lost. Gets rid of the thousands of spurious low-abundance
#'                LKT fragments and nothing else. "auto" is a synonym.
#'   "light"      Denoise plus 5 further percentage points of median sequencing depth spent on a
#'                cleaner feature set. At a native taxonomic rank this buys a genome-completeness
#'                floor, which denoise alone does not impose. THE DEFAULT of every plotting
#'                function.
#'   "moderate"   As light, but 10 points.
#'   "stringent"  As light, but 20 points.
#'   "presence"   The denoise tier with the prevalence rule forced to a single sample. The
#'                feature must still clear a fitted abundance floor and, for taxa, a completeness
#'                floor, but it need only do so once. This is the engraftment case: a complete
#'                genome seen in one subject is the finding, and every percentage-based
#'                prevalence rule destroys it.
#'   NULL         No filtration whatsoever. An explicit, supported choice, warned about once per
#'                call because an unfiltered metagenomic feature universe is overwhelmingly
#'                low-abundance noise and a reader who is not a bioinformatician will take those
#'                rows as findings.
#'
#' WHY THE TIERS ARE A DEPTH BUDGET. Spending 5, 10 or 20 points of the median sample's retained
#' sequencing depth is as arbitrary a choice as 250 PPM was, but it is arbitrary in a unit the
#' user can reason about, it is computed against their own data, and the fit says what it bought.
#' The tiers are nested, so a stringent plot can never contain a feature a light plot dropped.
#'
#' PRECEDENCE. An explicitly supplied featcutoff or GenomeCompletenessCutoff always wins,
#' whatever applyfilters says and even when applyfilters is NULL. This is how a historical
#' analysis is reproduced now that the legacy constants are gone: read the numbers off the
#' subtitle of the old figure and pass them directly.
#'
#' @param SEobj The JAMS2 SummarizedExperiment object being filtered. Required for "auto" and
#'   "presence"; optional otherwise, in which case analysis and is16S are derived from it.
#' @param analysis Legacy alias for naming the analysis space when no SEobj is passed. Ignored
#'   when SEobj is supplied and carries metadata()$analysis.
#' @param is16S Whether this is a 16S object. Derived from SEobj when NULL.
#' @param glomby The taxonomic level the caller is about to agglomerate to, if any. Thresholds are
#'   rank-specific, so this selects which stamped entry is read and what is fitted on the fly.
#' @param applyfilters One of the strings above, or NULL.
#' @param featcutoff Explicit c(PPM, percent of samples). Overrides the preset.
#' @param GenomeCompletenessCutoff Explicit c(percent complete, percent of samples). Overrides.
#' @param PPM_normalize_to_bases_sequenced Which PPM denominator this call will filter against.
#'   Forwarded to the fit, and compared against the denominator a stamped fit used, because a
#'   PPM cutoff is not transferable between the two.
#' @param max_ranks_to_filter How many taxonomic ranks, counting up from the bottom, may be
#'   filtered at all. Default 5: LKT (or the native terminal rank), Species, Genus, Family,
#'   Order. At Class and above "auto" resolves to no filtration and logs why.
#' @param max_ranks_for_completeness How many ranks, counting up from the bottom, may use the
#'   genome completeness axis. Default 3: LKT, Species and Genus. See resolve_filtration_rank.
#' @param fit_if_missing If TRUE (default), "auto" fits on the fly when the object carries no
#'   stamped entry for the rank in question. If FALSE, a missing entry means no filtration.
#' @param verbose Log the resolution and its justification.
#' @param ... Passed to fit_filtering_thresholds when a fit happens. Stale effect-size arguments
#'   are named in a warning rather than silently swallowed.
#' @export

declare_filtering_presets <- function(SEobj = NULL, analysis = NULL, is16S = NULL, glomby = NULL, applyfilters = NULL, featcutoff = NULL, GenomeCompletenessCutoff = NULL, PPM_normalize_to_bases_sequenced = FALSE, max_ranks_to_filter = 5, max_ranks_for_completeness = 3, fit_if_missing = TRUE, verbose = TRUE, ...){

    have_SEobj <- (!is.null(SEobj)) && (as.character(class(SEobj)[1]) == "SummarizedExperiment")

    #Derive what the object can tell us, so that the eight call sites do not each have to.
    if (have_SEobj){
        if (!is.null(metadata(SEobj)$analysis)){
            analysis <- metadata(SEobj)$analysis
        }
        if (is.null(is16S)){
            is16S <- ("asv2lkt" %in% names(metadata(SEobj)))
        }
    }
    if (is.null(is16S)){
        is16S <- FALSE
    }

    taxonomic_spaces <- c("LKT", "Contig_LKT", "ConsolidatedGenomeBin", "MB2bin", "16S")
    #Resolved once. A bare `analysis %in% taxonomic_spaces` is logical(0) when analysis is NULL,
    #and `if (logical(0))` is a hard error, so every branch below would abort on a NULL analysis.
    is_taxonomic <- ((!is.null(analysis)) && (length(analysis) == 1) && (analysis %in% taxonomic_spaces))

    #Named, rather than silently swallowed, so that a stale caller is visible in the log
    #instead of quietly having its effect-size threshold ignored. Arguments destined for
    #fit_filtering_thresholds are passed through rather than reported.
    dots <- list(...)
    deadargs <- names(dots)
    if (is.null(deadargs)){
        deadargs <- character(0)
    }
    deadargs <- deadargs[nzchar(deadargs)]
    effectargs <- intersect(deadargs, c("minl2fc", "maxl2fc", "minabscorrcoeff"))
    if (length(effectargs) > 0){
        flog.warn(paste0("declare_filtering_presets no longer handles effect-size thresholds, so the following argument(s) were ignored: ", paste0(effectargs, collapse = ", "), ". Filtration presets now govern feature survival only. Set these on the plotting function itself, where they have their own defaults (minl2fc = 1, minabscorrcoeff = 0.4)."))
        dots <- dots[setdiff(names(dots), effectargs)]
    }
    #Anything left in ... is destined for fit_filtering_thresholds. Arguments it does not accept
    #are reported and dropped rather than forwarded, because do.call would turn a typo into an
    #"unused argument" error from a function the user never called by name.
    if (length(dots) > 0){
        fitformals <- names(formals(fit_filtering_thresholds))
        otherargs <- setdiff(names(dots), fitformals)
        if (length(otherargs) > 0){
            flog.warn(paste0("declare_filtering_presets does not recognize the following argument(s), which were ignored: ", paste0(otherargs, collapse = ", "), "."))
            dots <- dots[intersect(names(dots), fitformals)]
        }
    }

    if ((!is_taxonomic) && (!(is.null(GenomeCompletenessCutoff)))){
        flog.warn("Genome completeness only makes sense for taxa. Please choose a taxonomic (non functional) analysis.")
        GenomeCompletenessCutoff <- NULL
    }

    presetlist <- list()
    presetlist$applyfilters_resolved <- applyfilters
    presetlist$fit <- NULL

    #Validate the string up front. An unrecognised preset used to fall through every branch and
    #silently become no filtration at all, which is the worst of the available behaviours: the
    #user believes they filtered and they did not.
    validpresets <- c("denoise", "auto", "light", "moderate", "stringent", "presence")
    if (!is.null(applyfilters)){
        if ((length(applyfilters) != 1) || (!is.character(applyfilters)) || (is.na(applyfilters)) || (!(applyfilters %in% validpresets))){
            flog.warn(paste0("applyfilters must be one of ", paste0(paste0('"', validpresets, '"'), collapse = ", "), ", or NULL for no filtration. The supplied value (", paste0(as.character(applyfilters), collapse = ", "), ") is not recognised and was ignored, so NO filtration will be applied."))
            applyfilters <- NULL
            presetlist$applyfilters_resolved <- NULL
        }
    }

    if (is.null(applyfilters)){
        #Deliberate, supported, and worth a line in the log. Silence here is how an unfiltered
        #feature universe gets mistaken for a biological result. It is only worth saying when the
        #caller supplied no explicit cutoffs either: filter_experiment is called dozens of times
        #per figure with concrete featcutoff values and a NULL preset, and announcing the obvious
        #each time would bury the messages that matter.
        if (verbose && is.null(featcutoff) && is.null(GenomeCompletenessCutoff)){
            flog.info("applyfilters is NULL and no explicit cutoffs were given, so no filtration will be applied. Every feature, including low-abundance and low-completeness ones, will be reported. If that was not intended, use applyfilters = \"light\" (the default of the plotting functions) or \"denoise\" for the gentlest fitted thresholds.")
        }
    } else {

        #Every preset is now a tier of the same fitted ladder. "presence" is the denoise tier with
        #the prevalence rule forced to a single sample; "auto" is a synonym for "denoise".
        wanted_tier <- if (identical(applyfilters, "presence")) "denoise" else if (identical(applyfilters, "auto")) "denoise" else applyfilters

        if (!have_SEobj){
            flog.warn(paste0("applyfilters = \"", applyfilters, "\" fits thresholds to the data in the SummarizedExperiment object, but no object was passed to declare_filtering_presets, so NO filtration will be applied. This is a bug in the calling function; please report it."))
        } else {
            #One resolver decides the rank, its nativity and whether it is deep enough to filter,
            #so that declare_filtering_presets, stamp_filtering_thresholds and
            #fit_filtering_thresholds cannot disagree. In particular it reads the rank off
            #rowData rather than off metadata$analysis, which agglomerate_features leaves naming
            #the native rank even on a glommed object, and it treats a terminal glomby ("LKT" on
            #a ConsolidatedGenomeBin object) as the native rank, because agglomerate_features
            #returns such a request unchanged.
            rankinfo <- resolve_filtration_rank(SEobj = SEobj, glomby = glomby, max_ranks_to_filter = max_ranks_to_filter, max_ranks_for_completeness = max_ranks_for_completeness)
            rankkey <- rankinfo$rank
            stamped <- metadata(SEobj)$taxon_filtration_thresholds
            thisfit <- NULL

            #Too high up the ladder to filter at all. Said once, here, so that the user is told
            #why their Phylum heatmap is unfiltered rather than left to wonder.
            if (!rankinfo$filterable){
                if (verbose){
                    flog.info(rankinfo$reason)
                }
                rankkey <- NULL
            }

            #The cache is keyed [[rank]][[tier]]. Objects stamped by JAMS 2.3.3 used a flat
            #[[rank]] holding a single fit, which is detectable because a fit has a featcutoff
            #and a tier list does not. Such an entry is the denoise tier under the old name, so
            #it is honoured for that tier and ignored for the others rather than silently
            #serving denoise thresholds to someone who asked for stringent.
            rankentry <- if (!is.null(rankkey)) stamped[[rankkey]] else NULL
            legacy_flat <- (is.list(rankentry) && ("featcutoff" %in% names(rankentry)))
            if (legacy_flat){
                if (identical(wanted_tier, "denoise")){
                    flog.warn(paste0("The filtration thresholds stamped on this object for rank ", rankkey, " are in the pre-tier format written by JAMS 2.3.3. They are being read as the \"denoise\" tier. Run stamp_filtering_thresholds(overwrite = TRUE) to record the full light/moderate/stringent ladder."))
                } else {
                    flog.warn(paste0("The filtration thresholds stamped on this object for rank ", rankkey, " are in the pre-tier format written by JAMS 2.3.3 and describe the \"denoise\" tier only, so tier \"", wanted_tier, "\" was refitted on the fly. Run stamp_filtering_thresholds(overwrite = TRUE) to record the full ladder and stop paying for this on every call."))
                }
                rankentry <- list(denoise = rankentry)
            }

            have_stamped <- ((!is.null(rankkey)) && is.list(rankentry) && (wanted_tier %in% names(rankentry)) && is.list(rankentry[[wanted_tier]]))

            #Does a stamped entry actually describe the object in hand? The sample set is one
            #question and the FEATURE set is a different, sharper one. agglomerate_features copies
            #metadata() wholesale and leaves metadata$analysis naming the native rank, so a user
            #who agglomerates by hand and then plots with glomby = NULL presents an object whose
            #rankkey still says "ConsolidatedGenomeBin" while its rows are genera. Reading the
            #stamped native thresholds there is exactly the cross-rank inheritance this design
            #forbids, and the sample names are identical so nothing else would catch it. The row
            #count is recorded at fit time for this purpose.
            if (have_stamped){
                fitted_nfeat <- rankentry[[wanted_tier]]$provenance$nfeat_object
                if ((!is.null(fitted_nfeat)) && is.finite(fitted_nfeat) && (fitted_nfeat != nrow(SEobj))){
                    flog.warn(paste0("Filtration thresholds stamped for rank ", rankkey, " were fitted on an object with ", fitted_nfeat, " features, but the object in hand has ", nrow(SEobj), ". They therefore do not describe this feature universe - the usual cause is an object agglomerated before the plotting call, which keeps the native metadata$analysis - so they were ignored and the thresholds refitted. Pass glomby to the plotting function rather than pre-agglomerating, and the correct rank will be used."))
                    have_stamped <- FALSE
                }
            }

            #A PPM threshold is only meaningful against the denominator it was fitted on, and the
            #two denominators differ: TotalBasesSequenced counts everything sequenced,
            #TotalBasesSequencedinAnalysis only the bases assigned to features in this analysis.
            #The objects are stamped with FALSE (the titration default, matching the plotting
            #functions), while filter_experiment, compare_within_subsets and export_expvec_to_XL
            #default to TRUE. The flags are therefore compared - but only after checking that the
            #two denominators actually differ on this object, because for a taxonomic space
            #make_SummarizedExperiments sets them equal, which makes the flag irrelevant and a
            #warning on every call pure noise.
            if (have_stamped){
                fitted_ppmnorm <- rankentry[[wanted_tier]]$provenance$PPM_normalize_to_bases_sequenced
                denominators_differ <- (!isTRUE(all.equal(metadata(SEobj)$TotalBasesSequenced, metadata(SEobj)$TotalBasesSequencedinAnalysis)))
                if ((!is.null(fitted_ppmnorm)) && (!identical(isTRUE(fitted_ppmnorm), isTRUE(PPM_normalize_to_bases_sequenced))) && denominators_differ){
                    flog.warn(paste0("Filtration thresholds stamped for rank ", rankkey, " were fitted with PPM_normalize_to_bases_sequenced = ", isTRUE(fitted_ppmnorm), " but this call uses ", isTRUE(PPM_normalize_to_bases_sequenced), ", and the two PPM denominators differ on this object. A PPM cutoff is not transferable between them, so the thresholds were refitted."))
                    have_stamped <- FALSE
                }
            }

            #"presence" changes the prevalence rule, so a stamped entry fitted at the default
            #prevalence is not an answer to the question being asked and is deliberately not read.
            if ((!identical(applyfilters, "presence")) && have_stamped){
                candidate <- rankentry[[wanted_tier]]
                #The thresholds were fitted on a particular set of samples; if the object in hand
                #is not that set (subsetby, samplesToKeep, a merge that slipped through) the
                #numbers may be stale. Say so and use them anyway rather than refitting, because
                #refitting per subset would make subsets incomparable with each other, which is a
                #subtler and worse failure than a slightly stale cutoff.
                fitted_samples <- candidate$provenance$samplenames
                nowsamples <- colnames(SEobj)
                if ((!is.null(fitted_samples)) && (!setequal(fitted_samples, nowsamples))){
                    flog.warn(paste0("Filtration thresholds stamped on this ", analysis, " object at rank ", rankkey, " were fitted on ", length(fitted_samples), " samples, but the object now in hand has ", length(nowsamples), " (", length(intersect(fitted_samples, nowsamples)), " in common). The stamped thresholds are being used as they stand, so that every subset of this cohort is filtered identically and remains comparable. Re-run stamp_filtering_thresholds(overwrite = TRUE) if you want them refitted to the current sample set."))
                }
                thisfit <- candidate
                if (verbose){
                    flog.info(paste0("Using filtration thresholds stamped on the object for rank ", rankkey, ", tier \"", wanted_tier, "\": ", candidate$justification))
                }
            } else if (fit_if_missing && (!is.null(rankkey))){
                #Duplicated names would make do.call abort with "formal argument matched by
                #multiple actual arguments", so the ones this function sets itself are stripped
                #from the forwarded dots first.
                #glomby is deliberately NOT forwarded. By the time any plotting function calls
                #this, ExpObjVetting has already agglomerated the object it hands over, so asking
                #the fitter to agglomerate again would either double-glom or fail. The rank is
                #communicated as a label, and whether the completeness axis is usable comes from
                #the resolver's structural reading of the object rather than from glomby, which
                #is blind to an object the user agglomerated before calling.
                selfset <- c("SEobj", "glomby", "rankname", "fit_completeness", "max_ranks_to_filter", "max_ranks_for_completeness", "tier", "PPM_normalize_to_bases_sequenced", "verbose")
                fitargs <- c(list(SEobj = SEobj, glomby = NULL, rankname = rankkey, fit_completeness = rankinfo$completeness_fittable, max_ranks_to_filter = max_ranks_to_filter, max_ranks_for_completeness = max_ranks_for_completeness, tier = wanted_tier, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, verbose = verbose), dots[setdiff(names(dots), selfset)])
                if (identical(applyfilters, "presence")){
                    #One sample is the whole point of this preset; a caller-supplied value in ...
                    #would defeat it, so it is set last and wins.
                    fitargs$prevalence_min_samples <- 1
                }
                thisfit <- try(do.call(fit_filtering_thresholds, fitargs), silent = TRUE)
                if (inherits(thisfit, "try-error")){
                    flog.warn(paste0("Could not fit filtration thresholds for rank ", rankkey, ", so NO filtration will be applied. Consider passing featcutoff explicitly. The error was: ", as.character(thisfit)))
                    thisfit <- NULL
                } else if (verbose && (!have_stamped) && (!identical(applyfilters, "presence"))){
                    flog.info(paste0("This object carries no usable stamped filtration thresholds for rank ", rankkey, ", tier \"", wanted_tier, "\", so they were fitted on the fly. Run expvec <- stamp_filtering_thresholds(expvec) once to record them on the objects and skip this step in future."))
                }
            } else if (!is.null(rankkey)) {
                flog.warn(paste0("No filtration thresholds are stamped on this object for rank ", rankkey, " and fit_if_missing is FALSE, so NO filtration will be applied."))
            }

            if (!is.null(thisfit)){
                presetlist$fit <- thisfit
                presetlist$featcutoff <- thisfit$featcutoff
                presetlist$GenomeCompletenessCutoff <- thisfit$GenomeCompletenessCutoff
                presetlist$fitmsg <- thisfit$justification
            }
        }
    }

    #Replace with any values explicitly set by the user. This block runs after EVERY branch
    #above, including the NULL one, which is what makes the precedence rule unconditional: an
    #explicit featcutoff or GenomeCompletenessCutoff always wins, whatever applyfilters said and
    #even when applyfilters is NULL. That is the documented way to reproduce a historical
    #analysis, now that the hard-coded light/moderate/stringent constants are gone: read the
    #numbers off the subtitle of the old figure and pass them here.
    argstoset <- c("featcutoff", "GenomeCompletenessCutoff")[!unlist(lapply(list(featcutoff, GenomeCompletenessCutoff), is.null))]

    if (length(argstoset) > 0){
        for (ats in argstoset){
            presetlist[[ats]] <- get(ats)
        }
    }

    #Generate a filtration message. This string is what reaches plot subtitles, so it has to
    #describe the filter that will actually run. Two traps. First, a resolved featcutoff of
    #c(0, 0) means NO filtration - filter_experiment skips the block entirely - so announcing
    #"Feature must be >=0 PPM in at least 0% of samples" would caption a plot with a filter that
    #was never applied. That case now yields a NULL message, as it did before "auto" existed and
    #featcutoff could arrive here as NULL. Second, an absolute-sample-count prevalence is carried
    #as a fractional percentage (the half-sample margin described in fit_filtering_thresholds),
    #so "in at least 0.5% of samples" is arithmetically right and humanly useless; where the fit
    #recorded a sample count, the message says so.
    presetlist$filtermsg <- NULL
    prevalence_phrase <- function(pct){
        k <- presetlist$fit$rule$prevalence_min_samples
        if ((!is.null(k)) && is.finite(k)){
            return(paste0("at least ", k, " sample(s)"))
        }

        return(paste0("at least ", signif(pct, 4), "% of samples"))
    }

    if (is.null(presetlist$featcutoff)){
        presetlist$featcutoff <- c(0, 0)
    }

    if (!all(presetlist$featcutoff == c(0, 0))){
        presetlist$thresholdPPM <- presetlist$featcutoff[1]
        presetlist$sampcutoffpct <- min(presetlist$featcutoff[2], 100)
        #">=" and not ">": filter_experiment keeps a feature whose passing proportion is >= this,
        #and its per-sample test is also >=. The message used to say otherwise.
        presetlist$filtermsg <- paste0("Feature must be >=", presetlist$thresholdPPM, " PPM in ", prevalence_phrase(presetlist$sampcutoffpct))
    }

    if (!(is.null(presetlist$GenomeCompletenessCutoff))){
        if ((!is16S) && (!all(presetlist$GenomeCompletenessCutoff == c(0, 0)))){
            presetlist$thresholdGenomeCompleteness <- presetlist$GenomeCompletenessCutoff[1]
            presetlist$sampcutoffpctGenomeCompleteness <- min(presetlist$GenomeCompletenessCutoff[2], 100)
            gcmsg <- paste0("Taxon genome completeness must be >=", presetlist$thresholdGenomeCompleteness, "% in ", prevalence_phrase(presetlist$sampcutoffpctGenomeCompleteness))
            presetlist$filtermsg <- if (is.null(presetlist$filtermsg)) gcmsg else paste(presetlist$filtermsg, gcmsg, sep = "\n")
        }
    }

    return(presetlist)
}
