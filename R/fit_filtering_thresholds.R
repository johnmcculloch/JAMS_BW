#' resolve_filtration_rank(SEobj = NULL, glomby = NULL, max_ranks_to_filter = 5, max_ranks_for_completeness = 3)
#'
#' Works out which taxonomic rank a filtration threshold would actually apply to, and whether
#' filtering at that rank is worth doing at all. Used by fit_filtering_thresholds,
#' stamp_filtering_thresholds and declare_filtering_presets so that all three agree.
#'
#' WHY THE RANK CANNOT BE READ OFF metadata$analysis OR glomby. agglomerate_features copies
#' metadata() wholesale, so an object agglomerated to Genus still reports
#' metadata$analysis == "ConsolidatedGenomeBin". A caller that trusts that name will label a
#' genus-level fit as native, and - worse - will treat the summed GenomeCompleteness assay as a
#' percentage. What agglomerate_features DOES change is rowData: it keeps the taxonomic columns
#' down to glomby and drops everything below, so the deepest taxonomic column present is a
#' reliable, structural statement of what the rows actually are. A native object still has LKT;
#' a genus-level one stops at Genus.
#'
#' WHY ONLY THE FIVE LOWEST RANKS ARE FILTERED. Filtration exists to remove spurious
#' low-abundance assignments, and that is a problem of the deep ranks, where a cohort carries
#' thousands of features of which most are noise. By Class and above there are only tens of
#' features, every one of them a real clade, and a threshold chosen by any data-driven criterion
#' discards whole clades rather than noise: run against the MR50 benchmark, an unrestricted fit
#' kept 1 Domain out of 4 and 2 Kingdoms out of 17, which is arithmetically defensible and
#' biologically absurd. Above Order the honest answer is to show everything.
#'
#' @param SEobj A JAMS2 SummarizedExperiment object.
#' @param glomby The rank the caller is about to agglomerate to, or NULL. A terminal level
#'   (LKT, ConsolidatedGenomeBin, Contig_LKT, MB2bin, 16S) is treated as NULL, because
#'   agglomerate_features returns such a request unchanged.
#' @param max_ranks_to_filter How many ranks, counting up from the bottom, may be filtered.
#'   Default 5, i.e. LKT (or the native terminal rank), Species, Genus, Family, Order.
#' @param max_ranks_for_completeness How many ranks, counting up from the bottom, may use the
#'   genome completeness axis. Default 3, i.e. LKT, Species and Genus.
#' @return A list with rank (the label to key the cache by), slot (the position on the taxonomic
#'   ladder, with the native terminal rank reported as "LKT"), is_native, is_taxonomic,
#'   filterable, reason (populated when filterable is FALSE), completeness_fittable, and
#'   completeness_reason (populated when completeness_fittable is FALSE).
#' @export

resolve_filtration_rank <- function(SEobj = NULL, glomby = NULL, max_ranks_to_filter = 5, max_ranks_for_completeness = 3){

    #Lowest first. This is the order filtration cares about, and the first max_ranks_to_filter of
    #it is the filterable set.
    ranks_low_to_high <- c("LKT", "Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Domain", "Gram")
    #Levels agglomerate_features refuses to agglomerate by, returning the object unchanged.
    terminal_levels <- c("LKT", "Contig_LKT", "ConsolidatedGenomeBin", "MB2bin", "16S", "IS1", "is1_")
    taxonomic_spaces <- c("LKT", "Contig_LKT", "ConsolidatedGenomeBin", "MB2bin", "16S")

    analysis <- metadata(SEobj)$analysis
    is_taxonomic <- ((!is.null(analysis)) && (length(analysis) == 1) && (analysis %in% taxonomic_spaces))

    out <- list(rank = analysis, slot = NA_character_, is_native = TRUE, is_taxonomic = is_taxonomic, filterable = TRUE, reason = NULL, completeness_fittable = is_taxonomic, completeness_reason = NULL)

    if (!is_taxonomic){
        #Functional spaces have one rank and no ladder. Whether they are filtered at all is
        #decided by allow_functional in fit_filtering_thresholds, not here.
        out$completeness_fittable <- FALSE

        return(out)
    }

    #A terminal glomby is a no-op in agglomerate_features, so it must not be treated as a
    #separate rank. Stamping one would duplicate the native fit under a second key, and looking
    #one up would miss the native entry.
    if ((!is.null(glomby)) && (length(glomby) == 1) && (!(glomby %in% terminal_levels))){
        out$rank <- glomby
        out$slot <- glomby
        out$is_native <- FALSE
    } else {
        #No glomby (or a terminal one): read the rank off the object itself.
        present <- ranks_low_to_high[ranks_low_to_high %in% colnames(rowData(SEobj))]
        deepest <- if (length(present) > 0) present[1] else NA_character_
        if (is.na(deepest) || identical(deepest, "LKT")){
            #Still carries its leaves: this is the native object.
            out$rank <- analysis
            out$slot <- "LKT"
            out$is_native <- TRUE
        } else {
            #The leaf columns are gone, so somebody agglomerated this object before handing it
            #over. Believe the rowData, not metadata$analysis.
            out$rank <- deepest
            out$slot <- deepest
            out$is_native <- FALSE
        }
    }

    filterable_slots <- ranks_low_to_high[1:max(1, min(max_ranks_to_filter, length(ranks_low_to_high)))]
    if (!is.na(out$slot) && (!(out$slot %in% filterable_slots))){
        out$filterable <- FALSE
        out$reason <- paste0("Filtration is not applied at ", out$rank, ". Only the ", max_ranks_to_filter, " lowest taxonomic ranks (", paste0(filterable_slots, collapse = ", "), ") are filtered, because above those there are few features and every one of them is a real clade, so any threshold discards clades rather than noise. Everything is shown at this rank.")
    }

    #WHERE GENOME COMPLETENESS MAY BE USED. agglomerate_features aggregates the assay by SUMMING
    #it across child taxa, so above the native rank it is not strictly a percentage. That is a
    #reason for care, not for abandoning the axis, and the asymmetry is what matters:
    #
    #  Summation is MONOTONE. The sum is never smaller than the largest child, so an aggregated
    #  cutoff is strictly MORE PERMISSIVE than the same cutoff applied to the best-recovered
    #  child. A species whose CSB was 80% recovered cannot be removed by a 30% floor at Species.
    #  The filter can therefore fail to remove something, but it cannot remove something that had
    #  a lot of genome behind it, which is the error that would actually matter.
    #
    #  The converse risk is fragments adding up to a spurious "complete" genome, and that risk
    #  scales with FAN-OUT - how many children a node has. At Species the fan-out is small and
    #  the children are species-coherent by construction, because contextualize_taxonomy builds
    #  CSBs within a single species-level NCBI taxid; the pieces being summed are genuinely
    #  pieces of one organism, so the sum is a better estimate of recovery than any single child.
    #  At Genus the fan-out is still modest and summation is what recovers the case JAMS exists
    #  to catch: a novel species absent from the Kraken2 database, unbinned by MetaBAT2, whose
    #  contigs scatter across several pseudo-species of a real genus and only add up to a
    #  recognisable genome once agglomerated.
    #
    #  From Family upward the fan-out grows without bound and summation stops meaning anything,
    #  so the axis is stood down there. max_ranks_for_completeness is where that line sits.
    #
    #NOTE that enabling the axis does not force its use: the fit only moves along completeness if
    #doing so removes features cheaply in depth, so on a cohort where aggregated completeness is
    #uninformative the criterion simply declines to use it.
    completeness_slots <- ranks_low_to_high[1:max(1, min(max_ranks_for_completeness, length(ranks_low_to_high)))]
    if (!is.na(out$slot) && (!(out$slot %in% completeness_slots))){
        out$completeness_fittable <- FALSE
        out$completeness_reason <- paste0("Genome completeness was not fitted at rank ", out$rank, " because agglomerate_features sums that assay across child taxa, and only the ", max_ranks_for_completeness, " lowest ranks (", paste0(completeness_slots, collapse = ", "), ") have a small enough fan-out for the sum to still describe one organism's genome. The fit used the PPM and prevalence axes only.")
    }
    if (!out$filterable){
        out$completeness_fittable <- FALSE
    }

    return(out)
}


#' fit_filtering_thresholds(SEobj = NULL, glomby = NULL, rankname = NULL, fit_completeness = NULL, tier = "denoise", tier_keep_fraction = c(denoise = 1, light = 0.5, moderate = 0.25, stringent = 0.10), tier_depth_rails = c(denoise = 80, light = 80, moderate = 75, stringent = 65), min_noise_fraction_for_PPM_axis = 0.05, only_allow_CSBs = FALSE, samplesToKeep = NULL, prevalence_pct = NULL, prevalence_min_samples = 2, min_depth_median = 80, min_depth_min = NULL, min_depth_lowq = NULL, depth_quantile_probs = 0.1, divergence_tolerance = 2, max_ranks_to_filter = 5, max_ranks_for_completeness = 3, min_features_to_filter = 50, allow_functional = FALSE, PPM_thresholds = NULL, completeness_thresholds = NULL, n_ppm = 30, n_comp = 30, max_PPM_threshold = 10000, ignoreunclassified = TRUE, PPM_normalize_to_bases_sequenced = TRUE, class_to_ignore = "N_A", precomputed_titration = NULL, return_plots = FALSE, verbose = TRUE)
#'
#' Fits filtration thresholds to the data actually present in a JAMS2 SummarizedExperiment object,
#' replacing the arbitrary "light" / "moderate" / "stringent" presets with a choice that can be
#' justified from the object itself.
#'
#' WHY THIS EXISTS. Until JAMS 2.3.2 the only way to decide a feature cutoff was to pick one of
#' three hard-coded preset strings whose numbers were chosen by eye, years ago, on other people's
#' data. Those numbers were never wrong so much as unaccountable: nothing about c(250, 15) refers
#' to the cohort it is applied to. At the same time, running with no filtration at all is not a
#' neutral act either, because the unfiltered feature universe of a deeply sequenced metagenome is
#' overwhelmingly composed of spurious low-abundance assignments, and a user who is not a
#' bioinformatician will read those as findings. titrate_feature_counts was written to let a human
#' see, for their own data, where tightening the thresholds stops removing noise and starts
#' removing community. This function encodes that same read as an algorithm so that it can be the
#' default rather than an expert manoeuvre.
#'
#' THE SELECTION RULE, in three steps, over the grid returned by titrate_feature_counts:
#'
#' 1. SAFETY RAIL. A cell is admissible only if the surviving features still carry at least
#'    min_depth_median percent of the median sample's sequencing depth. This is a rail, not the
#'    selector: it exists to exclude catastrophe, and the choosing is done in step 2 by a
#'    criterion that already penalises depth loss. An earlier version made the rail the selector
#'    (median 95%, plus a floor on the WORST sample at 85%) and the MR50 benchmark showed why
#'    that fails. Cohorts routinely contain blanks, controls and poorly sequenced samples holding
#'    almost no real community; such a sample legitimately loses most of its depth to any
#'    threshold at all, and on MR50 the worst sample fell to 74% at the very first non-trivial
#'    cell. A worst-sample floor therefore let two or three blanks veto the whole scheme,
#'    returning "no filtration" for an 8625-LKT object from one mouse room on one diet.
#'    min_depth_min is consequently NULL by default: the worst sample is reported in every
#'    justification but constrains nothing unless asked for. min_depth_lowq is the middle road, a
#'    floor on a low quantile of per-sample retention rather than on the single worst sample, so
#'    a handful of blanks cannot veto while a systematic massacre still can. The
#'    (0 PPM, 0 completeness) cell always retains 100% of everything, so the admissible set is
#'    never empty and this function never fails to return an answer.
#'
#' 2. SELECTION. Among admissible cells, maximise divergence = pct_depth_median - pct_retained,
#'    the quantity plotted by the count_vs_depth_divergence panel. Read plainly it is "how many
#'    more percentage points of sequencing depth did I keep than of features", i.e. how much
#'    noise was axed per unit of real community lost. Its maximum is normally interior rather
#'    than at the tight corner, because past a point each further notch starts costing depth.
#'
#' 3. BACK-OFF. Among admissible cells within divergence_tolerance percentage points of that
#'    maximum, take the LOOSEST rather than the best-scoring one. This is the one-standard-error
#'    rule familiar from cross-validated model selection: when several candidates are
#'    statistically indistinguishable, prefer the gentler one. It stops the fit from sliding into
#'    the tight corner on cohorts where depth retention holds up all the way across the grid.
#'
#' 4. TIER. Steps 1-3 produce the "denoise" cell: the one that removes the most noise per unit
#'    of real community lost. The tiers above it walk the same surface further, each aiming to
#'    RETAIN A TARGET FRACTION of the denoise cell's feature count - a half, a quarter, a tenth
#'    by default. "Light halves the feature list" is a sentence a biologist can act on. Tiers are
#'    walked cumulatively, each constrained by the one below it, so a stringent plot can never
#'    contain a feature a light plot dropped. Each tier has its own depth rail
#'    (80/80/75/65): denoise and light are what a user gets without asking and keep the strict
#'    rail, while reaching moderate or stringent is a deliberate act.
#'
#'    This replaced a depth BUDGET ("spend N points of median retained depth"), which read
#'    better but does not survive a quantile-spaced ladder: rungs chosen to remove comparable
#'    numbers of features remove wildly uneven amounts of depth, so one rung routinely cost more
#'    than an entire tier budget and the tiers collapsed onto each other.
#'
#'    Worked example, the MR50 benchmark at its native rank: denoise keeps 409 features of 8625
#'    at 7% completeness, light 190 at 49%, moderate 95 at 92%. Restricted to Consolidated
#'    Species Bins, where the abundance axis stands itself down, it is 87 / 64 / 44 / 33 bins at
#'    90 / 96 / 98 / 99% completeness.
#'
#' ONLY THE FIVE LOWEST RANKS ARE FILTERED, and only ranks with at least min_features_to_filter
#' features. See resolve_filtration_rank for the reasoning: filtration removes spurious
#' low-abundance assignments, which is a problem of the deep ranks, and at Class and above every
#' feature is a real clade, so any data-driven threshold prunes clades rather than noise.
#' Unrestricted, the fit kept 1 Domain out of 4 on MR50.
#'
#' WHAT IT DELIBERATELY DOES NOT DO. It does not optimise prevalence. Prevalence is the parameter
#' with the largest effect on power and the least visibility, and silently raising it is how a
#' legitimately rare feature - a complete genome seen in one subject, which is the whole point of
#' an engraftment study - disappears without anyone being told. Prevalence is therefore a declared
#' axis, held fixed across the grid, and the fit works within it. It also never sees the group
#' labels. That is not an oversight but the condition that makes fitting thresholds on the same
#' data that is later tested legitimate: a filter chosen from marginal abundance and depth alone
#' is independent of the test statistic under the null, which is the independent-filtering
#' argument of Bourgon, Gentleman and Huber (PNAS 2010). Never fit within compareby strata.
#'
#' GENOME COMPLETENESS IS FITTED AT LKT, SPECIES AND GENUS. agglomerate_features aggregates the
#' GenomeCompleteness assay by SUMMING it across child taxa, exactly as it does BaseCounts, which
#' is a reason for care but not for abandoning the axis above the native rank. The asymmetry is
#' what licenses it: summation is monotone, so the sum is never smaller than the largest child,
#' and an aggregated cutoff is therefore strictly MORE PERMISSIVE than the same cutoff applied to
#' the best-recovered child. A species whose CSB was 80% recovered cannot be removed by a 30%
#' floor at Species. The filter can fail to remove something; it cannot remove something with a
#' well-recovered genome behind it, and that asymmetry is the opposite of the one on the PPM
#' axis, where aggregation actively manufactures artefacts.
#'
#' The converse risk - sparse fragments summing to a figure none of them earns - scales with
#' FAN-OUT, which is why the axis is capped at the three lowest ranks. At Species the fan-out is
#' small and the children are species-coherent by construction: contextualize_taxonomy builds
#' CSBs within a single species-level NCBI taxid, so the pieces being summed really are pieces of
#' one organism and the sum is a better estimate of recovery than any single child. At Genus the
#' fan-out is still modest and summation recovers the case JAMS exists to catch - a novel species
#' absent from the Kraken2 database and unbinned by MetaBAT2, whose contigs scatter across
#' several pseudo-species of a real genus and only add up to a recognisable genome once
#' agglomerated. From Family upward the fan-out is unbounded and the axis is stood down.
#'
#' Enabling the axis does not force its use: the fit only moves along completeness when doing so
#' removes features cheaply in depth, so where aggregated completeness is uninformative the
#' criterion simply declines to use it.
#'
#' Whether a rank is native is decided STRUCTURALLY, from the deepest taxonomic column in rowData
#' (see resolve_filtration_rank), and explicitly NOT from whether completeness exceeds 100%. That
#' test was tried and was wrong: make_SummarizedExperiments sums Completeness across the entities
#' consolidated into one ConsolidatedGenomeBin row, so the native MR50 assay reaches 188.5% with
#' nothing agglomerated, and evaluate_LKTs scores non-prokaryotic contigs by
#' percent-of-expected-genome-size rather than with CheckM2, which also overshoots before
#' reformat_completeness_to_CheckM_style caps it. Treating >100 as evidence of agglomeration cost
#' MR50 its most discriminating axis.
#'
#' A NOTE ON WHY INHERITING THRESHOLDS ACROSS RANKS IS UNSAFE. PPM aggregates by summation, so a
#' cutoff that a feature clears at LKT is cleared trivially at Genus. That is not conservative. A
#' genus assembled from two hundred spurious 3 PPM LKT assignments arrives at 600 PPM and passes a
#' floor of 50 that every one of its constituents failed. Aggregation does not merely weaken an
#' inherited filter, it manufactures the exact artefact the filter exists to prevent, and it does
#' so most aggressively for the ragbag taxa with the noisiest children. Fit each rank separately;
#' it costs a fraction of a second.
#'
#' @param SEobj A JAMS2 SummarizedExperiment object.
#' @param glomby Taxonomic level to agglomerate to before fitting. NULL (default) fits the object
#'   as handed over. Thresholds are rank-specific and must not be carried between ranks.
#' @param rankname Label for the rank being fitted, overriding what would be inferred from glomby.
#'   Needed because declare_filtering_presets receives an object that ExpObjVetting has already
#'   agglomerated: glomby must stay NULL there to avoid agglomerating twice, yet the rank being
#'   fitted is the glommed one.
#' @param fit_completeness Force the genome completeness axis on (TRUE) or off (FALSE). NULL
#'   (default) decides from resolve_filtration_rank, i.e. from whether the object's rows are
#'   still its native leaves. A caller which agglomerated before calling should pass FALSE.
#' @param tier One of "denoise", "light", "moderate", "stringent" ("auto" is a synonym for
#'   "denoise"). See step 4 above.
#' @param tier_keep_fraction Named numeric vector giving, per tier, the share of the denoise
#'   cell's surviving feature count that tier aims to retain. Default c(denoise = 1,
#'   light = 0.5, moderate = 0.25, stringent = 0.10), i.e. light halves the feature list.
#' @param tier_depth_rails Named numeric vector giving, per tier, the floor on the median
#'   sample's retained depth. Default c(denoise = 80, light = 80, moderate = 75,
#'   stringent = 65). The upper tiers are allowed to spend further because reaching them is a
#'   deliberate act, whereas denoise and light are what a user gets without asking. 65 rather
#'   than something lower because on a multi-cohort human stool object a floor of 55 let
#'   "stringent" discard 41% of the median sample's sequencing depth, and it separated the
#'   tiers no better than 65 did.
#' @param min_noise_fraction_for_PPM_axis The abundance axis is stood down when fewer than this
#'   fraction of features carry no genome evidence at all, i.e. when the universe has no noise
#'   population to remove. Default 0.05. Completeness is never stood down this way: it is a
#'   quality statement, abundance is not.
#' @param only_allow_CSBs Restrict to Consolidated Species Bins before fitting.
#' @param only_allow_CSBs Restrict to Consolidated Species Bins before fitting.
#' @param samplesToKeep Optional vector of samples to restrict the fit to. Recorded in the
#'   provenance so that a later mismatch is detectable.
#' @param prevalence_pct The prevalence rule held fixed across the grid, as a percentage of
#'   samples. Default 5. Not optimised; see above.
#' @param prevalence_min_samples Optional absolute sample count, overriding prevalence_pct. Set to
#'   1 for the "present anywhere, but only if the genome is real" rule used in engraftment work.
#' @param min_depth_median Safety rail on the median sample's retained depth, in percent.
#'   Default 80. Deliberately low: this excludes catastrophe, it does not choose the cell.
#' @param min_depth_min Optional floor on the WORST-hit sample's retained depth, in percent.
#'   Default NULL, i.e. reported but not enforced, because blanks and poorly sequenced controls
#'   legitimately lose most of their depth to any threshold and must not be able to veto a
#'   cohort's filtration. Set a number to enforce it.
#' @param min_depth_lowq Optional floor on a low QUANTILE of per-sample retained depth, in
#'   percent. Default NULL. The robust middle road between the median and the worst sample: a
#'   few blanks cannot veto, but depth loss across many samples still can.
#' @param depth_quantile_probs Which quantile pct_depth_lowq reports. Default 0.1.
#' @param max_ranks_to_filter How many taxonomic ranks, counting up from the bottom, may be
#'   filtered at all. Default 5: LKT (or the native terminal rank), Species, Genus, Family,
#'   Order. Above that the fit refuses and nothing is filtered.
#' @param max_ranks_for_completeness How many ranks, counting up from the bottom, may use the
#'   genome completeness axis. Default 3: LKT, Species, Genus. See resolve_filtration_rank for
#'   why summed completeness is a usable - and conservative - filter at those ranks and not
#'   above them.
#' @param min_features_to_filter Ranks with fewer starting features than this are not filtered.
#'   Default 50.
#' @param divergence_tolerance Back-off width in percentage points of divergence. Default 2.
#' @param allow_functional If FALSE (default), a non-taxonomic analysis space returns a
#'   no-filtration result with status "not_applicable". Functional spaces have no genome
#'   completeness axis and, more importantly, a rare ECNumber or Product is frequently a real gene
#'   that is genuinely rare rather than an artefact, so an abundance floor is a poor proxy for
#'   annotation error. Set TRUE to fit the PPM axis there anyway.
#' @param PPM_thresholds,completeness_thresholds,n_ppm,n_comp,max_PPM_threshold Grid controls,
#'   passed to titrate_feature_counts. max_PPM_threshold (default 10000 PPM, i.e. 1% relative
#'   abundance) caps the auto-built abundance ladder, so the fit cannot impose an abundance
#'   floor so high that it is selecting dominant taxa rather than removing noise.
#' @param ignoreunclassified Exclude the unclassified bin when counting. Default TRUE, because
#'   dark matter should not drive the choice of threshold. Recorded in the provenance, since
#'   filter_experiment does NOT drop that row and the delivered feature count will therefore
#'   exceed the fitted one by however many unclassified rows survive.
#' @param PPM_normalize_to_bases_sequenced Passed to titrate_feature_counts.
#' @param class_to_ignore Metadata classes to drop (signature parity).
#' @param precomputed_titration A titration already computed for THIS rank and prevalence rule,
#'   so that stamp_filtering_thresholds can derive all four tiers from one grid instead of
#'   titrating four times. Internal convenience; supplying a grid from a different rank or
#'   prevalence would silently produce wrong thresholds.
#' @param return_plots If TRUE, returns the four titration panels with the chosen cell outlined,
#'   for auditing or for pasting into a report.
#' @param verbose Emit the decision and its justification via flog.info.
#' @return A list of class "JAMSfilterfit" carrying featcutoff, GenomeCompletenessCutoff, the
#'   chosen grid row, the full grid, the rule parameters, a human-readable justification, and a
#'   provenance record. Designed to be consumed by declare_filtering_presets and stored on the
#'   object by stamp_filtering_thresholds.
#' @export

fit_filtering_thresholds <- function(SEobj = NULL, glomby = NULL, rankname = NULL, fit_completeness = NULL, tier = "denoise", tier_keep_fraction = c(denoise = 1, light = 0.5, moderate = 0.25, stringent = 0.10), tier_depth_rails = c(denoise = 80, light = 80, moderate = 75, stringent = 65), min_noise_fraction_for_PPM_axis = 0.05, only_allow_CSBs = FALSE, samplesToKeep = NULL, prevalence_pct = NULL, prevalence_min_samples = 2, min_depth_median = 80, min_depth_min = NULL, min_depth_lowq = NULL, depth_quantile_probs = 0.1, divergence_tolerance = 2, max_ranks_to_filter = 5, max_ranks_for_completeness = 3, min_features_to_filter = 50, allow_functional = FALSE, PPM_thresholds = NULL, completeness_thresholds = NULL, n_ppm = 30, n_comp = 30, max_PPM_threshold = 10000, ignoreunclassified = TRUE, PPM_normalize_to_bases_sequenced = TRUE, class_to_ignore = "N_A", precomputed_titration = NULL, return_plots = FALSE, verbose = TRUE){

    assert_JAMS2_SEobj(SEobj = SEobj, objname = "SEobj", calling_function = "fit_filtering_thresholds")

    #"auto" is retained as a synonym for the gentlest tier, because it is the obvious thing for a
    #user to type and it is what this argument was called before the tiers existed.
    if (identical(tier, "auto")){
        tier <- "denoise"
    }
    tier_names <- names(tier_keep_fraction)
    if ((length(tier) != 1) || (!is.character(tier)) || (!(tier %in% tier_names))){
        stop(paste0("tier must be one of ", paste0(paste0('"', tier_names, '"'), collapse = ", "), " (or \"auto\", a synonym for \"denoise\")."))
    }
    if ((!is.numeric(tier_keep_fraction)) || any(!is.finite(tier_keep_fraction)) || any(tier_keep_fraction <= 0) || any(tier_keep_fraction > 1)){
        stop("tier_keep_fraction must be a named numeric vector of fractions in (0, 1], giving the share of the denoise feature count each tier aims to retain.")
    }
    if (!all(tier_names %in% names(tier_depth_rails))){
        stop("tier_depth_rails must name every tier in tier_keep_fraction.")
    }
    tier_index <- match(tier, tier_names)
    #Retained for the record and for the message strings.
    tier_budget <- as.numeric(tier_keep_fraction[[tier]])
    analysis <- metadata(SEobj)$analysis

    #Everything about which rank this is, whether it is the native one, and whether it is deep
    #enough to be worth filtering, comes from one place. See resolve_filtration_rank for why the
    #answer is read off rowData rather than off metadata$analysis or glomby.
    rankinfo <- resolve_filtration_rank(SEobj = SEobj, glomby = glomby, max_ranks_to_filter = max_ranks_to_filter, max_ranks_for_completeness = max_ranks_for_completeness)
    #rankname is a label the caller may override: declare_filtering_presets receives an object
    #that ExpObjVetting has already agglomerated, so glomby must stay NULL here (agglomerating a
    #second time would be wrong) while the rank being fitted is nonetheless the glommed one.
    rankname <- if (!is.null(rankname)) rankname else rankinfo$rank
    is_taxonomic <- rankinfo$is_taxonomic
    #fit_completeness is the caller's override; otherwise the resolver decides, from the rank's
    #position on the ladder rather than from nativity. See its header for why Species and Genus
    #are included and Family upward is not.
    if (is.null(fit_completeness)){
        fit_completeness <- rankinfo$completeness_fittable
    }

    #A no-filtration result, used for the functional opt-out and as the shape every caller can
    #rely on receiving. featcutoff c(0, 0) is the documented "do nothing" value in
    #filter_experiment, so a consumer which simply uses these numbers does the right thing
    #without having to test the status field first.
    null_fit <- function(status = "not_applicable", justification = ""){
        list(featcutoff = c(0, 0), GenomeCompletenessCutoff = NULL, tier = tier, status = status, justification = justification, chosen = NULL, grid = NULL, plots = NULL,
             rule = list(tier = tier, tier_keep_fraction = tier_budget, prevalence_pct = prevalence_pct, prevalence_min_samples = prevalence_min_samples, min_depth_median = min_depth_median, min_depth_min = min_depth_min, min_depth_lowq = min_depth_lowq, divergence_tolerance = divergence_tolerance, max_ranks_to_filter = max_ranks_to_filter, max_ranks_for_completeness = max_ranks_for_completeness, min_features_to_filter = min_features_to_filter),
             provenance = list(analysis = analysis, rank = rankname, slot = rankinfo$slot, glomby = glomby, only_allow_CSBs = only_allow_CSBs, samplenames = colnames(SEobj), nsamp = ncol(SEobj), nfeat_total = NA_integer_, nfeat_object = nrow(SEobj), ignoreunclassified = ignoreunclassified, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, completeness_fitted = FALSE, JAMSversion = as.character(metadata(SEobj)$version)[1], fitted_on = as.character(Sys.time())))
    }

    if ((!is_taxonomic) && (!allow_functional)){
        msg <- paste0("fit_filtering_thresholds: ", analysis, " is a functional analysis space, so no thresholds were fitted and no filtration will be applied. Functional spaces carry no genome completeness axis, and a rare feature there is frequently a genuinely rare gene rather than an artefact, so an abundance floor is a poor proxy for annotation error. Pass allow_functional = TRUE to fit the PPM axis anyway, or set featcutoff explicitly.")
        if (verbose){
            flog.info(msg)
        }

        return(structure(null_fit(status = "not_applicable", justification = msg), class = c("JAMSfilterfit", "list")))
    }

    #Too high up the taxonomic ladder to be worth filtering. Refused before the titration rather
    #than after it, because the grid is not the problem - the rank is.
    if (!rankinfo$filterable){
        if (verbose){
            flog.info(paste0("fit_filtering_thresholds: ", rankinfo$reason))
        }

        return(structure(null_fit(status = "not_applicable", justification = rankinfo$reason), class = c("JAMSfilterfit", "list")))
    }

    #Run the titration. return_plots is driven by what the caller asked for, since building the
    #four ggplots is the only part of this that is not free.
    #precomputed_titration lets stamp_filtering_thresholds titrate a rank ONCE and derive all
    #four tiers from that single grid: the grid depends on the rank and the prevalence rule, not
    #on the tier, so re-running it per tier would quadruple the work (and the agglomerations)
    #for identical numbers. Supplying a titration computed for a different rank or prevalence
    #would silently produce wrong thresholds, so it is an internal convenience, not an invitation.
    if (!is.null(precomputed_titration)){
        #Shape-checked rather than trusted: titrate_feature_counts returns NULL when it was run
        #with plot_on_the_fly = TRUE, and a silent NULL here would surface several steps later
        #as a "feature universe is empty" refusal blaming ignoreunclassified.
        ok_shape <- is.data.frame(precomputed_titration) || (is.list(precomputed_titration) && is.data.frame(precomputed_titration$grid))
        if (!ok_shape){
            stop("precomputed_titration must be the value returned by titrate_feature_counts with return_plots = FALSE (a data frame) or TRUE (a list with a $grid data frame). A titration run with plot_on_the_fly = TRUE returns NULL and cannot be reused.")
        }
        titration <- precomputed_titration
    } else {
        titration <- titrate_feature_counts(ExpObj = SEobj, glomby = glomby, only_allow_CSBs = only_allow_CSBs, samplesToKeep = samplesToKeep, PPM_thresholds = PPM_thresholds, completeness_thresholds = completeness_thresholds, n_ppm = n_ppm, n_comp = n_comp, max_PPM_threshold = max_PPM_threshold, prevalence_pct = prevalence_pct, prevalence_min_samples = prevalence_min_samples, depth_quantile_probs = depth_quantile_probs, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, ignoreunclassified = ignoreunclassified, class_to_ignore = class_to_ignore, return_plots = return_plots, plot_on_the_fly = FALSE, verbose = FALSE)
    }

    #A titration carries its grid directly when plots were not asked for, and inside $grid when
    #they were. Tested on shape rather than on the return_plots flag, because a precomputed
    #titration may have been built either way.
    if (is.data.frame(titration)){
        grid <- titration
        plots <- NULL
    } else {
        grid <- titration$grid
        plots <- titration$plots
    }

    universe <- attr(grid, "titration_universe")
    if (is.null(universe)){
        universe <- list(nfeat_total = NA_integer_, nfeat_object = nrow(SEobj), nsamp = ncol(SEobj), samplenames = colnames(SEobj), prev_desc = paste0(">= ", prevalence_pct, "% of samples"), have_GenomeCompleteness = FALSE, min_samples = NA_real_)
    }

    #Degenerate universe guard. If the starting feature set is empty - everything dropped by the
    #ignoreunclassified filter, say - every pct_retained is 100 * 0/0 = NaN, and the selection
    #below then takes max() of all-NA (returning -Inf with a warning), finds no contenders, pulls
    #row 1 of a zero-row frame (a phantom all-NA row) and dies at `if (NA)`. Caught here instead,
    #where the reason can be stated.
    degenerate <- ((!is.finite(universe$nfeat_total)) || (universe$nfeat_total < 1) || all(!is.finite(grid$pct_retained)) || all(!is.finite(grid$pct_depth_median)))
    if (degenerate){
        msg <- paste0("fit_filtering_thresholds: the feature universe for ", rankname, " is empty or carries no sequencing depth, so there is nothing to titrate and no filtration will be applied. Check whether ignoreunclassified removed every row, or whether these samples have any bases assigned to this analysis at all.")
        flog.warn(msg)

        return(structure(null_fit(status = "refused", justification = msg), class = c("JAMSfilterfit", "list")))
    }

    #Too few features at this rank for any threshold to be denoising rather than pruning. The
    #rank cap above handles the usual case; this catches a shallow cohort whose Family or Order
    #space is as sparse as a Class space normally is.
    if (is.finite(universe$nfeat_total) && (universe$nfeat_total < min_features_to_filter)){
        msg <- paste0("fit_filtering_thresholds: ", rankname, " has only ", universe$nfeat_total, " features, fewer than min_features_to_filter (", min_features_to_filter, "), so no filtration will be applied. At this few features a data-driven threshold removes real taxa rather than noise.")
        if (verbose){
            flog.info(msg)
        }

        return(structure(null_fit(status = "not_applicable", justification = msg), class = c("JAMSfilterfit", "list")))
    }

    #Is the completeness axis interpretable at all? Two reasons it might not be: the object may
    #not carry the assay, and the rank may be high enough up the ladder that the summed assay no
    #longer describes one organism's genome. resolve_filtration_rank owns that second decision
    #and carries the reasoning; this function only honours it, or the caller's override.
    #
    #WHAT IS DELIBERATELY *NOT* USED AS EVIDENCE: a completeness above 100%. That was tried and is
    #wrong. JAMS produces values above 100 at the native rank by two legitimate routes. First,
    #make_SummarizedExperiments sums Completeness across the entities consolidated into a single
    #ConsolidatedGenomeBin row, so a row assembled from a MetaBAT2 bin plus taxonomy-only contigs
    #can exceed 100 with nothing agglomerated at all - on the MR50 benchmark the native assay
    #reaches 188.5%. Second, evaluate_LKTs scores non-prokaryotic contigs with a
    #"Percent contig sum over expected genome size" model instead of CheckM2, and although
    #reformat_completeness_to_CheckM_style caps that at 100 and tips the excess into
    #Contamination, it only runs when such taxa are present and older objects predate it.
    #Treating >100 as "this has been agglomerated" cost MR50 its completeness axis, which is the
    #most discriminating axis it has.
    completeness_available <- isTRUE(universe$have_GenomeCompleteness)
    completeness_reason <- NULL
    #Recorded so that the justification can say what a cutoff above the native rank means.
    completeness_is_summed <- (completeness_available && isTRUE(fit_completeness) && (!isTRUE(rankinfo$is_native)))

    if (completeness_available && (!isTRUE(fit_completeness))){
        completeness_available <- FALSE
        completeness_reason <- if (!is.null(rankinfo$completeness_reason)) rankinfo$completeness_reason else paste0("Genome completeness was not fitted at rank ", rankname, " because the caller stood that axis down. The fit used the PPM and prevalence axes only.")
    }

    if (completeness_available){
        observed_max_comp <- suppressWarnings(max(grid$Completeness, na.rm = TRUE))
        #A flat completeness axis (a single level) is itself the tell that there was nothing to
        #titrate: the assay is absent, constant, or all zero.
        if (!is.finite(observed_max_comp) || (length(unique(grid$Completeness)) < 2)){
            completeness_available <- FALSE
            completeness_reason <- paste0("Genome completeness was not fitted at rank ", rankname, " because the completeness ladder collapsed to a single level, which means the assay is absent, constant or zero throughout.")
        }
    }
    if ((!completeness_available) && (!is.null(completeness_reason)) && verbose){
        flog.info(completeness_reason)
    }

    #If completeness is not in play, restrict the search to the completeness == 0 row so that the
    #fit cannot choose a cell it is not allowed to express as a cutoff.
    searchgrid <- grid
    if (!completeness_available){
        searchgrid <- searchgrid[searchgrid$Completeness == min(searchgrid$Completeness), , drop = FALSE]
    }
    #Computed on the whole search space rather than only on the admissible subset, because the
    #tier escalation below walks cells the depth rail excludes.
    searchgrid$divergence <- searchgrid$pct_depth_median - searchgrid$pct_retained

    #STAND THE ABUNDANCE AXIS DOWN ON A UNIVERSE THAT CONTAINS NO NOISE.
    #Divergence has no notion that a feature is REAL. It sees only depth, so on a universe where
    #every feature is a genuine organism it will happily spend its budget removing rare ones and
    #report a fine score for doing so. The discriminator is in the data and was going unused: the
    #fraction of features carrying NO genome evidence at all. Where that fraction is large there
    #is a noise population and abundance filtering is denoising; where it is near zero - a
    #Consolidated Species Bin universe, where MetaBAT2 binning has already done the denoising -
    #there is nothing to denoise, and an abundance floor is just a choice to keep the dominant
    #members of the community. On MR50 restricted to CSBs the figure is exactly 0%, and the
    #abundance axis was nonetheless being driven to a 1% relative-abundance floor.
    #Completeness is never stood down this way: it is a quality statement, abundance is not.
    ppm_axis_reason <- NULL
    if (completeness_available && is.finite(universe$n_no_evidence) && (universe$nfeat_total > 0)){
        noise_fraction <- universe$n_no_evidence / universe$nfeat_total
        if (noise_fraction < min_noise_fraction_for_PPM_axis){
            searchgrid <- searchgrid[searchgrid$PPM == min(searchgrid$PPM), , drop = FALSE]
            ppm_axis_reason <- paste0("The abundance axis was stood down at ", rankname, " because only ", sprintf("%.1f%%", 100 * noise_fraction), " of its features carry no genome evidence at all, i.e. this universe has no noise population to remove - MetaBAT2 binning or an upstream restriction has already done that job. Filtering it further on abundance would select the dominant members of the community rather than remove noise, so the fit used genome completeness alone.")
            if (verbose){
                flog.info(ppm_axis_reason)
            }
        }
    }

    #STEP 1: feasibility. This is a SAFETY RAIL, not the selector. The original version made it
    #the selector - median >= 95% and worst sample >= 85% - and the MR50 benchmark showed why
    #that is wrong twice over. The worst-sample floor is the more serious error: a cohort
    #normally contains blanks, controls and poorly sequenced samples which hold almost no real
    #community, and such a sample legitimately loses most of its depth to any threshold at all.
    #On MR50 the very first non-trivial cell (1 PPM) took the worst sample to 74.3%, so a floor
    #of 85 vetoed the entire surface and the fit returned "no filtration" on an object with 8625
    #LKTs from one mouse room on one diet. Two or three blanks should not be able to veto a
    #cohort's filtration scheme. The worst-sample floor is therefore OFF by default (reported,
    #never enforced, unless the user asks for it), and min_depth_median has come down to a level
    #at which it only excludes catastrophe, leaving the divergence criterion - which already
    #penalises depth loss - to do the actual choosing. min_depth_lowq is the middle road: a floor
    #on a low quantile rather than on the single worst sample, so a handful of blanks cannot
    #veto, but a systematic massacre across many samples still can.
    depth_ok <- (!is.na(searchgrid$pct_depth_median)) & (searchgrid$pct_depth_median >= min_depth_median)
    if (!is.null(min_depth_min)){
        depth_ok <- depth_ok & (is.na(searchgrid$pct_depth_min) | (searchgrid$pct_depth_min >= min_depth_min))
    }
    if ((!is.null(min_depth_lowq)) && ("pct_depth_lowq" %in% colnames(searchgrid))){
        depth_ok <- depth_ok & (is.na(searchgrid$pct_depth_lowq) | (searchgrid$pct_depth_lowq >= min_depth_lowq))
    }
    feasible <- searchgrid[which(depth_ok), , drop = FALSE]

    #The (0, 0) cell retains everything and is always in the grid, so this should be unreachable.
    #It is guarded anyway because an empty feasible set would otherwise surface as an NA subscript
    #several lines later, which is a far worse way to find out.
    if (nrow(feasible) < 1){
        msg <- paste0("fit_filtering_thresholds: no grid cell met the depth rail (median >= ", min_depth_median, "%", (if (!is.null(min_depth_min)) paste0(", worst sample >= ", min_depth_min, "%") else ""), (if (!is.null(min_depth_lowq)) paste0(", low quantile >= ", min_depth_lowq, "%") else ""), ") for ", rankname, ", not even the unfiltered cell. This should be impossible and suggests the depth accounting is degenerate (all-zero samples?). No filtration will be applied.")
        flog.warn(msg)

        return(structure(null_fit(status = "refused", justification = msg), class = c("JAMSfilterfit", "list")))
    }

    #STEP 2: maximise divergence among admissible cells.
    best_divergence <- max(feasible$divergence, na.rm = TRUE)

    #STEP 3: back off to the loosest cell that is within tolerance of the best.
    ppm_levels <- sort(unique(grid$PPM))
    comp_levels <- sort(unique(grid$Completeness))
    loosest_of <- function(df){
        df$looseness <- match(df$PPM, ppm_levels) + match(df$Completeness, comp_levels)
        df <- df[order(df$looseness, df$PPM, df$Completeness), , drop = FALSE]

        return(df[1, , drop = FALSE])
    }
    contenders <- feasible[which(feasible$divergence >= (best_divergence - divergence_tolerance)), , drop = FALSE]
    #"Loosest" = lowest combined position on the two ladders, tie-broken towards the lower PPM.
    denoise_row <- loosest_of(contenders)
    chosen <- denoise_row
    escalation_note <- NULL

    #STEP 4: TIER ESCALATION. "denoise" stops at the cell chosen above - the one that removes the
    #most noise per unit of real community lost. The tiers above it are not different criteria
    #but the same surface walked further, each aiming to RETAIN A TARGET FRACTION of the denoise
    #cell's surviving feature count: a half, a quarter, a tenth. "Light halves the feature list"
    #is a sentence a biologist can act on.
    #
    #WHY THIS IS NOT A DEPTH BUDGET ANY MORE. The tiers used to spend a fixed number of
    #percentage points of median retained depth, which read better but does not survive a
    #quantile-spaced ladder. Quantile rungs are chosen so that each removes a comparable number
    #of FEATURES, and the depth those features carry is wildly uneven, so a single rung routinely
    #costs more depth than an entire tier budget. The escalation then had nowhere to land and the
    #tiers collapsed onto each other - on MR50, denoise, light and moderate all returned the same
    #cell at several ranks. Targeting features instead puts the tier and the ladder in the same
    #currency. The depth rail is what stops this being reckless, and it is the rail, not the
    #tier, that carries the "do not gut the cohort" guarantee.
    #
    #THE LADDER IS WALKED CUMULATIVELY, each tier nesting against the cell chosen for the tier
    #below it rather than against the denoise cell. Nesting every tier against denoise
    #independently does NOT give a monotone ladder, which was an earlier version of this and was
    #wrong: "loosest" is the sum of two ladder indices, so a tighter target can move the winner
    #sideways across the (PPM, completeness) plane instead of outward, and a feature dropped by
    #light could then reappear under moderate. Constraining each step by its predecessor makes "a
    #stringent plot can never contain a feature a light plot dropped" true by construction.
    #
    #THE RAIL RELAXES FOR THE UPPER TIERS. min_depth_median protects the default, which is now
    #"light" in every plotting function, so a user reaches denoise and light without asking for
    #them. Asking for moderate or stringent is a deliberate act, and on a cohort whose depth is
    #concentrated in a few features a fixed rail would otherwise pin all three upper tiers to the
    #same cell. tier_depth_rails carries one rail per tier.
    if (tier_index > 1){
        ladder_tiers <- tier_names[2:tier_index]
        prev_row <- denoise_row
        n0 <- as.numeric(denoise_row$n_features[1])
        for (thistier in ladder_tiers){
            thisrail <- as.numeric(tier_depth_rails[[thistier]])
            if (!is.finite(thisrail)){
                thisrail <- min_depth_median
            }
            target_n <- n0 * as.numeric(tier_keep_fraction[[thistier]])
            #Cells at least as strict as the tier below, strictly stricter on at least one axis,
            #and still clearing this tier's rail.
            nested <- searchgrid[which((searchgrid$PPM >= prev_row$PPM[1]) & (searchgrid$Completeness >= prev_row$Completeness[1]) & ((searchgrid$PPM > prev_row$PPM[1]) | (searchgrid$Completeness > prev_row$Completeness[1])) & (!is.na(searchgrid$pct_depth_median)) & (searchgrid$pct_depth_median >= thisrail)), , drop = FALSE]
            if (nrow(nested) < 1){
                break
            }
            #The cell whose surviving count sits closest to the target, loosest among ties.
            gap <- abs(nested$n_features - target_n)
            prev_row <- loosest_of(nested[which(gap <= (min(gap) + 1e-9)), , drop = FALSE])
        }
        chosen <- prev_row
        #The escalation can legitimately fail to move, when nothing stricter clears the rail.
        #That must be said rather than left for the justification to claim a step that never
        #happened.
        if ((chosen$PPM[1] == denoise_row$PPM[1]) && (chosen$Completeness[1] == denoise_row$Completeness[1])){
            escalation_note <- paste0("Tier \"", tier, "\" aims to retain about ", round(100 * tier_budget), "% of the denoise feature count, but no cell stricter than the one denoising chose still clears that tier's depth rail, so the denoise thresholds were used unchanged. This rank cannot be filtered harder without discarding more depth than the rail allows.")
            if (verbose){
                flog.info(escalation_note)
            }
        }
    }

    chosen_ppm <- as.numeric(chosen$PPM[1])
    chosen_comp <- as.numeric(chosen$Completeness[1])

    #The prevalence rule, expressed as the percentage that filter_experiment expects. This is
    #READ BACK from the titration rather than recomputed here, because the two must describe the
    #same rule to the sample and a second derivation is a second opportunity to disagree. The
    #titration places the percentage halfway between the counts it has to separate (counts are
    #integers, so ((k - 0.5) / nsamp) * 100 admits exactly k or more and cannot be flipped by a
    #floating-point round trip); see its header.
    nsamp <- universe$nsamp
    prevalence_as_pct <- universe$prevalence_as_pct
    if (is.null(prevalence_as_pct) || (!is.finite(prevalence_as_pct))){
        k <- max(1, min(if (!is.null(prevalence_min_samples)) prevalence_min_samples else 2, nsamp))
        prevalence_as_pct <- ((k - 0.5) / nsamp) * 100
    }

    unfiltered <- ((chosen_ppm <= 0) && (chosen_comp <= 0))

    if (unfiltered){
        featcutoff <- c(0, 0)
        GenomeCompletenessCutoff <- NULL
    } else {
        featcutoff <- c(chosen_ppm, prevalence_as_pct)
        if (completeness_available && (chosen_comp > 0)){
            GenomeCompletenessCutoff <- c(chosen_comp, prevalence_as_pct)
        } else {
            GenomeCompletenessCutoff <- NULL
        }
    }

    #Diagnostics for the justification string. The decline rate is reported rather than enforced:
    #the depth rail plus the back-off already do the constraining, and a second tunable
    #interacting with the first is how a selection rule becomes unexplainable. It is reported
    #without alarm, because on a log-spaced PPM ladder every step is large and a double-digit
    #rate at a perfectly good cell is the norm rather than a warning sign.
    chosen_rate <- suppressWarnings(as.numeric(chosen$decline_rate[1]))
    rate_desc <- if (is.na(chosen_rate)) "the unfiltered corner, where a local rate of loss is undefined" else paste0("a local feature-loss rate of ", sprintf("%.0f%%", 100 * chosen_rate), " per grid step (", (if (chosen_rate <= 0.1) "a plateau" else "a sloping part of the surface"), ")")
    depth_desc <- paste0(chosen$pct_depth_median[1], "% of the median sample's sequencing depth", (if (is.finite(chosen$pct_depth_min[1])) paste0(" and ", chosen$pct_depth_min[1], "% of the worst sample's (not a constraint unless min_depth_min is set)") else ""))

    if (unfiltered){
        justification <- paste0("fit_filtering_thresholds chose NO filtration for ", rankname, " at tier \"", tier, "\". Every cell that removed any features fell below the depth rail (median >= ", min_depth_median, "%", (if (!is.null(min_depth_min)) paste0(", worst sample >= ", min_depth_min, "%") else ""), (if (!is.null(min_depth_lowq)) paste0(", low quantile >= ", min_depth_lowq, "%") else ""), "), which means this cohort's features cannot be thinned without discarding real community. Treat the unfiltered feature count as including low-abundance noise. If a few blanks or poorly sequenced samples are responsible, look at the depth_retained panel and consider lowering min_depth_median or dropping those samples with samplesToKeep.")
        flog.warn(justification)
    } else {
        tier_desc <- if ((tier_index > 1) && is.null(escalation_note)) paste0(" Tier \"", tier, "\" then walked on from the denoise cell (", denoise_row$PPM[1], " PPM", (if (denoise_row$Completeness[1] > 0) paste0(" / ", denoise_row$Completeness[1], "% completeness") else ""), ", ", denoise_row$n_features[1], " features) aiming to retain about ", round(100 * tier_budget), "% of that count.") else ""
        justification <- paste0("fit_filtering_thresholds chose ", chosen_ppm, " PPM", (if (!is.null(GenomeCompletenessCutoff)) paste0(" and ", chosen_comp, "% genome completeness") else ""), ", each in ", universe$prev_desc, ", for ", rankname, " at tier \"", tier, "\". That cell keeps ", chosen$n_features[1], " of ", universe$nfeat_total, " features (", chosen$pct_retained[1], "%) while retaining ", depth_desc, ", a divergence of ", sprintf("%+.0f", chosen$divergence[1]), " points, at ", rate_desc, ".", tier_desc, " The denoise cell was the loosest within ", divergence_tolerance, " points of the best divergence (", sprintf("%.0f", best_divergence), ") among ", nrow(feasible), " cells clearing the depth rail.")
        if (verbose){
            flog.info(justification)
        }
    }

    if (!is.null(escalation_note)){
        justification <- paste(justification, escalation_note, sep = " ")
    }

    if (!is.null(completeness_reason)){
        justification <- paste(justification, completeness_reason, sep = " ")
    }
    if (completeness_is_summed && (!is.null(GenomeCompletenessCutoff))){
        justification <- paste(justification, paste0("At ", rankname, " the GenomeCompleteness assay is the SUM of its child taxa, so this cutoff reads as \"at least ", chosen_comp, "% of a genome recovered in total across the children of this ", rankname, "\". Because summing can only raise a value, the cutoff is strictly more permissive than the same number applied to the best-recovered child: nothing with a well-recovered genome behind it can be removed by it. What it can do is let many sparse fragments add up to a figure none of them earns, which is why the axis is stood down from Family upward."), sep = " ")
    }
    if (isTRUE(ignoreunclassified)){
        justification <- paste(justification, "Counts exclude the unclassified bin, which filter_experiment does not remove, so the delivered feature count will exceed the fitted one by however many unclassified rows survive.", sep = " ")
    }

    #Outline the chosen cell on each returned panel, so that a fit pasted into a report shows the
    #decision in the context of the surface it was taken from rather than as a bare number.
    if (return_plots && (!is.null(plots))){
        chosen_marker <- data.frame(PPM_f = factor(chosen_ppm, levels = ppm_levels), Comp_f = factor(chosen_comp, levels = comp_levels))
        plots <- lapply(plots, function(p){ p + geom_tile(data = chosen_marker, mapping = aes(x = PPM_f, y = Comp_f), fill = NA, colour = "black", linewidth = 0.9, inherit.aes = FALSE) })
    }

    fitobj <- list(featcutoff = featcutoff,
                   GenomeCompletenessCutoff = GenomeCompletenessCutoff,
                   tier = tier,
                   status = if (unfiltered) "unfiltered" else "fitted",
                   justification = justification,
                   chosen = chosen,
                   denoise_cell = denoise_row,
                   grid = grid,
                   plots = plots,
                   rule = list(tier = tier, tier_keep_fraction = tier_budget, tier_depth_rails = tier_depth_rails, min_noise_fraction_for_PPM_axis = min_noise_fraction_for_PPM_axis, prevalence_pct = prevalence_pct, prevalence_min_samples = prevalence_min_samples, prevalence_as_pct = prevalence_as_pct, min_depth_median = min_depth_median, min_depth_min = min_depth_min, min_depth_lowq = min_depth_lowq, depth_quantile_probs = depth_quantile_probs, divergence_tolerance = divergence_tolerance, max_ranks_to_filter = max_ranks_to_filter, max_ranks_for_completeness = max_ranks_for_completeness, min_features_to_filter = min_features_to_filter, best_divergence = best_divergence, n_feasible = nrow(feasible)),
                   provenance = list(analysis = analysis, rank = rankname, slot = rankinfo$slot, glomby = glomby, only_allow_CSBs = only_allow_CSBs, samplenames = universe$samplenames, nsamp = universe$nsamp, nfeat_total = universe$nfeat_total, nfeat_object = universe$nfeat_object, ignoreunclassified = ignoreunclassified, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, completeness_fitted = (!is.null(GenomeCompletenessCutoff)), completeness_is_summed = completeness_is_summed, JAMSversion = as.character(metadata(SEobj)$version)[1], fitted_on = as.character(Sys.time())))

    return(structure(fitobj, class = c("JAMSfilterfit", "list")))
}


#' print.JAMSfilterfit(x, ...)
#'
#' Prints a fitted filtration threshold object as the decision plus the evidence for it. The grid
#' is shown deliberately: a threshold the user cannot interrogate is indistinguishable from a
#' magic number, which is the problem this machinery exists to solve.
#'
#' @param x A JAMSfilterfit object from fit_filtering_thresholds.
#' @param show_grid Print the titration grid beneath the decision. Default TRUE.
#' @param ... Ignored.
#' @export

print.JAMSfilterfit <- function(x, show_grid = TRUE, ...){
    cat("JAMS fitted filtration thresholds\n")
    cat(paste0("  Rank:              ", x$provenance$rank, " (", x$provenance$analysis, ")\n"))
    cat(paste0("  Tier:              ", (if (!is.null(x$tier)) x$tier else "denoise"), "\n"))
    cat(paste0("  Status:            ", x$status, "\n"))
    cat(paste0("  featcutoff:        c(", paste0(signif(x$featcutoff, 6), collapse = ", "), ")\n"))
    if (!is.null(x$GenomeCompletenessCutoff)){
        cat(paste0("  GenomeCompleteness: c(", paste0(signif(x$GenomeCompletenessCutoff, 6), collapse = ", "), ")\n"))
    } else {
        cat("  GenomeCompleteness: not fitted\n")
    }
    cat(paste0("  Fitted on:         ", x$provenance$nsamp, " samples, ", x$provenance$nfeat_total, " starting features, ", x$provenance$fitted_on, "\n"))
    cat("\n")
    cat(strwrap(x$justification, width = 100, prefix = "  "), sep = "\n")
    if (show_grid && (!is.null(x$grid))){
        cat("\n  Grid it chose from:\n")
        gridcols <- intersect(c("PPM", "Completeness", "n_features", "pct_retained", "pct_depth_median", "pct_depth_lowq", "pct_depth_min", "decline_rate"), colnames(x$grid))
        print(x$grid[ , gridcols])
    } else if (show_grid){
        #stamp_filtering_thresholds keeps one copy of the grid per rank, on the first tier
        #stamped, because every tier of a rank is chosen from the identical surface and a
        #hundred-row data frame per tier per rank is real memory on an object that is saved and
        #reloaded.
        cat("\n  (The titration grid is shared across tiers and is stored on the first tier of this rank.)\n")
    }

    return(invisible(x))
}


#' stamp_filtering_thresholds(obj = NULL, ranks = NULL, tiers = c("denoise", "light", "moderate", "stringent"), max_ranks_to_filter = 5, max_ranks_for_completeness = 3, overwrite = FALSE, verbose = TRUE, ...)
#'
#' Fits filtration thresholds and writes them onto the object, at
#' metadata(SEobj)$taxon_filtration_thresholds, as a list keyed FIRST by rank and then by tier:
#' metadata(SEobj)$taxon_filtration_thresholds$Genus$moderate. Accepts either a single JAMS2
#' SummarizedExperiment object or an expvec (a named list of them) and returns the same shape,
#' so the idiomatic call is expvec <- stamp_filtering_thresholds(expvec).
#'
#' Each rank is TITRATED ONCE and all four tiers are derived from that single grid, because the
#' grid depends on the rank and the prevalence rule and not on the tier. The grid itself is kept
#' on the first tier of each rank only; storing a hundred-row data frame four times per rank per
#' analysis space is real memory on an object that gets saved and reloaded.
#'
#' WHY THE SLOT IS A CACHE AND NOT A REQUIREMENT. Objects built before JAMS 2.3.3 do not carry it,
#' and neither does an object the user built by hand. declare_filtering_presets therefore treats
#' an absent entry as an instruction to fit on the fly and say so in the log, so nothing breaks if
#' this function is never run. What running it buys is (a) not repeating the fit on every
#' plotting call and, more importantly, (b) reproducibility: the thresholds a figure used are
#' recorded on the object that produced it, rather than recomputed each time from whatever subset
#' of samples happened to be in play.
#'
#' AGGLOMERATION CARRIES THE SLOT. agglomerate_features copies metadata() wholesale onto the
#' glommed object, so a single stamped native object serves every glomby derived from it, provided
#' the rank in question was stamped. Ranks are therefore fitted and stored separately, never
#' inherited: see the header of fit_filtering_thresholds for why an LKT-level cutoff applied at
#' Genus is not conservative but actively harmful.
#'
#' MERGING DROPS THE SLOT. merge_JAMS_SEobj removes it, because thresholds fitted on cohort A say
#' nothing about cohort A plus B. Re-stamp after merging.
#'
#' @param obj A JAMS2 SummarizedExperiment object, or a named list of them (an expvec).
#' @param ranks Character vector of ranks to fit. NULL (default) fits the object's own rank plus
#'   every agglomerable taxonomic level that is both present in rowData and low enough to be
#'   filtered at all - by default Species, Genus, Family and Order. Class and above are skipped,
#'   as are terminal levels (fitting "LKT" on a ConsolidatedGenomeBin object just duplicates the
#'   native fit, since agglomerate_features returns terminal requests unchanged). "native" fits
#'   only the object's own rank.
#' @param tiers Which tiers to fit per rank. Default all four. All are derived from one
#'   titration per rank, so fitting all four costs barely more than fitting one.
#' @param max_ranks_to_filter Passed to resolve_filtration_rank and to each fit; how many ranks
#'   up from the bottom may be filtered. Default 5.
#' @param overwrite If FALSE (default), ranks already present on the object are left alone.
#' @param verbose Log each fit.
#' @param ... Passed to fit_filtering_thresholds, e.g. prevalence_pct, min_depth_median,
#'   min_depth_min, allow_functional.
#' @return The object, or list of objects, with metadata()$taxon_filtration_thresholds populated.
#' @export

stamp_filtering_thresholds <- function(obj = NULL, ranks = NULL, tiers = c("denoise", "light", "moderate", "stringent"), max_ranks_to_filter = 5, max_ranks_for_completeness = 3, overwrite = FALSE, verbose = TRUE, ...){

    #Dispatch over an expvec without requiring the caller to know which they are holding.
    if (as.character(class(obj)[1]) != "SummarizedExperiment"){
        if (is.list(obj) && (length(obj) > 0) && all(sapply(obj, function(x){ as.character(class(x)[1]) == "SummarizedExperiment" }))){
            #seq_along rather than names(): an unnamed expvec would make lapply(names(obj), ...)
            #iterate over NULL and silently return an empty list.
            stamped <- lapply(seq_along(obj), function(x){ stamp_filtering_thresholds(obj = obj[[x]], ranks = ranks, tiers = tiers, max_ranks_to_filter = max_ranks_to_filter, max_ranks_for_completeness = max_ranks_for_completeness, overwrite = overwrite, verbose = verbose, ...) })
            names(stamped) <- names(obj)

            return(stamped)
        } else {
            stop("stamp_filtering_thresholds requires a JAMS2 SummarizedExperiment object or a named list of them (an expvec).")
        }
    }

    assert_JAMS2_SEobj(SEobj = obj, objname = "obj", calling_function = "stamp_filtering_thresholds")

    analysis <- metadata(obj)$analysis
    nativeinfo <- resolve_filtration_rank(SEobj = obj, glomby = NULL, max_ranks_to_filter = max_ranks_to_filter, max_ranks_for_completeness = max_ranks_for_completeness)

    #Which ranks to fit. Three rules.
    #  (a) The object's own rank is always first, keyed by nativeinfo$rank, which is what
    #      declare_filtering_presets looks up when no glomby is given.
    #  (b) Only ranks that are actually filtered are fitted - the five lowest by default. Fitting
    #      Class or Phylum costs an agglomeration each and produces a refusal, so it is pure
    #      waste, and before the rank cap existed it produced nonsense (1 Domain kept out of 4).
    #  (c) Terminal levels are never fitted as separate ranks. agglomerate_features returns a
    #      request to glom by LKT unchanged, so stamping "LKT" on a ConsolidatedGenomeBin object
    #      duplicates the native fit under a second key and doubles the work for nothing.
    terminal_levels <- c("LKT", "Contig_LKT", "ConsolidatedGenomeBin", "MB2bin", "16S", "IS1", "is1_")
    if (is.null(ranks)){
        wanted <- nativeinfo$rank
        if (nativeinfo$is_taxonomic){
            agglomerable <- c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Domain")
            present <- agglomerable[agglomerable %in% colnames(rowData(obj))]
            #Keep only those the rank cap will actually allow, and never a terminal level.
            keepable <- present[sapply(present, function(x){ isTRUE(resolve_filtration_rank(SEobj = obj, glomby = x, max_ranks_to_filter = max_ranks_to_filter, max_ranks_for_completeness = max_ranks_for_completeness)$filterable) })]
            wanted <- unique(c(nativeinfo$rank, setdiff(keepable, terminal_levels)))
        }
    } else if (identical(ranks, "native")){
        wanted <- nativeinfo$rank
    } else {
        wanted <- unique(ranks)
    }

    existing <- metadata(obj)$taxon_filtration_thresholds
    if (is.null(existing)){
        existing <- list()
    }

    #Arguments this function sets itself must never also arrive through ..., or do.call aborts
    #with "formal argument matched by multiple actual arguments" - and because each fit is
    #wrapped in try(), that abort would surface as a warning and an unstamped object rather than
    #as an error, which is the worse outcome. stamp_filtering_thresholds(expvec, glomby = "Genus")
    #or a tier/tiers typo would otherwise silently stamp nothing at all.
    dots <- list(...)
    if (is.null(names(dots))){
        names(dots) <- rep("", length(dots))
    }
    selfset <- c("ExpObj", "SEobj", "obj", "glomby", "tier", "rankname", "max_ranks_to_filter", "precomputed_titration", "return_plots", "plot_on_the_fly", "verbose")
    clashing <- intersect(names(dots), selfset)
    if (length(clashing) > 0){
        flog.warn(paste0("stamp_filtering_thresholds sets the following argument(s) itself, so the value(s) supplied were ignored: ", paste0(clashing, collapse = ", "), ". Use the ranks and tiers arguments to control which ranks and tiers are fitted."))
        dots <- dots[setdiff(names(dots), selfset)]
    }
    titargs <- intersect(names(dots), names(formals(titrate_feature_counts)))

    for (rk in wanted){
        if ((rk %in% names(existing)) && (!overwrite) && (!is.list(existing[[rk]]) || (!("featcutoff" %in% names(existing[[rk]]))))){
            if (verbose){
                flog.info(paste0("stamp_filtering_thresholds: ", analysis, " already carries thresholds for rank ", rk, "; leaving them alone. Pass overwrite = TRUE to refit."))
            }
            next
        }
        if ((rk %in% names(existing)) && is.list(existing[[rk]]) && ("featcutoff" %in% names(existing[[rk]]))){
            #A pre-tier entry written by JAMS 2.3.3: a bare fit where a list of tiers belongs.
            #Refitted regardless of overwrite, because otherwise such an object could never be
            #upgraded by a plain stamp_filtering_thresholds(expvec) and every plotting call
            #would keep warning about it.
            flog.info(paste0("stamp_filtering_thresholds: rank ", rk, " of ", analysis, " carries pre-tier thresholds from JAMS 2.3.3; refitting it as the full tier ladder."))
        }
        #The object's own rank is fitted as it stands; anything else is a genuine agglomeration.
        glomby_arg <- if (identical(rk, nativeinfo$rank)) NULL else rk

        #Ask whether this rank will be fitted at all BEFORE paying for a titration. The fit
        #refuses functional spaces and ranks above the cap, and a functional space can carry
        #tens of thousands of rows, so titrating first and refusing afterwards would spend an
        #ExpObjVetting, a filter_experiment and a hundred cell counts to produce a refusal.
        thisrank <- resolve_filtration_rank(SEobj = obj, glomby = glomby_arg, max_ranks_to_filter = max_ranks_to_filter, max_ranks_for_completeness = max_ranks_for_completeness)
        allow_functional_here <- isTRUE(dots[["allow_functional"]])
        if (((!thisrank$is_taxonomic) && (!allow_functional_here)) || (!thisrank$filterable)){
            if (verbose){
                flog.info(paste0("stamp_filtering_thresholds: not fitting rank ", rk, " of ", analysis, ". ", (if (!thisrank$is_taxonomic) "It is a functional analysis space; pass allow_functional = TRUE to fit the PPM axis there." else thisrank$reason)))
            }
            next
        }

        #Titrate ONCE per rank and derive every tier from that single grid. The grid depends on
        #the rank and the prevalence rule, not on the tier, so titrating per tier would repeat
        #the agglomeration and the hundred cell counts four times over for identical numbers.
        thistitration <- try(do.call(titrate_feature_counts, c(list(ExpObj = obj, glomby = glomby_arg, return_plots = FALSE, plot_on_the_fly = FALSE, verbose = FALSE), dots[titargs])), silent = TRUE)
        if (inherits(thistitration, "try-error")){
            flog.warn(paste0("stamp_filtering_thresholds: could not titrate rank ", rk, " of ", analysis, ", so that rank was left unstamped and will be fitted on demand instead. The error was: ", as.character(thistitration)))
            next
        }

        tierfits <- list()
        for (tr in tiers){
            thisfit <- try(do.call(fit_filtering_thresholds, c(list(SEobj = obj, glomby = glomby_arg, tier = tr, max_ranks_to_filter = max_ranks_to_filter, max_ranks_for_completeness = max_ranks_for_completeness, precomputed_titration = thistitration, verbose = verbose), dots)), silent = TRUE)
            #A rank or tier can legitimately fail to fit. That must not cost the user the other
            #ranks or tiers, and it must not be silent either.
            if (inherits(thisfit, "try-error")){
                flog.warn(paste0("stamp_filtering_thresholds: could not fit tier \"", tr, "\" for rank ", rk, " of ", analysis, ", so that tier was left unstamped and will be fitted on demand instead. The error was: ", as.character(thisfit)))
                next
            }
            #The grid is identical across tiers, so it is kept once on the first tier only. A
            #hundred-row data frame per tier per rank per analysis space is real memory on an
            #object that gets saved to disk and loaded again.
            if (length(tierfits) > 0){
                thisfit$grid <- NULL
            }
            tierfits[[tr]] <- thisfit
        }

        if (length(tierfits) > 0){
            existing[[rk]] <- tierfits
        }
    }

    metadata(obj)$taxon_filtration_thresholds <- existing

    return(obj)
}
