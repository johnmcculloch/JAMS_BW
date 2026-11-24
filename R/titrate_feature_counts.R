#' titrate_feature_counts(ExpObj = NULL, glomby = NULL, only_allow_CSBs = FALSE, samplesToKeep = NULL, PPM_thresholds = NULL, completeness_thresholds = NULL, n_ppm = 30, n_comp = 30, max_PPM_threshold = 10000, prevalence_pct = NULL, prevalence_min_samples = 2, depth_quantile_probs = 0.1, PPM_normalize_to_bases_sequenced = TRUE, ignoreunclassified = TRUE, class_to_ignore = "N_A", return_plots = TRUE, plot_on_the_fly = FALSE, verbose = TRUE)
#'
#' A FAST, correlation-free companion to titrate_correlation_stability. For a taxonomic
#' SummarizedExperiment object, it counts how many features survive each combination of a
#' relative-abundance (PPM) threshold and a genome-completeness threshold, on a grid, at a fixed
#' prevalence. No correlations, no permutations: this is just repeated feature counting, so it runs
#' in a fraction of a second on a laptop even for thousands of features.
#'
#' The goal is to let a user see, for THEIR dataset, where the surviving-feature count stops moving
#' as thresholds tighten. Because raw counts decline smoothly (there is rarely a sharp elbow), the
#' function also reports the LOCAL RATE of decline (the marginal feature loss per grid step), which
#' is where turning points such as the recurring ~10% genome-completeness floor actually show up.
#'
#' Crucially, it ALSO reports how much SEQUENCING DEPTH the surviving features represent, per cell.
#' Feature-count loss and depth loss can diverge enormously: axing 90% of features might cost only
#' 5% of depth (you removed low-abundance noise) or 50% of depth (you removed real community),
#' depending on each sample's evenness. Depth is summarised as the MEDIAN retained fraction across
#' samples with the WORST-hit sample (min) alongside, and a divergence panel shows depth% minus
#' feature% directly, so a threshold that looks aggressive by count but cheap by depth is visible.
#'
#' Returned plots: feature_count, decline_rate, depth_retained, count_vs_depth_divergence. The
#' returned grid data frame carries n_features, pct_retained, pct_depth_median, pct_depth_lowq
#' and pct_depth_min.
#'
#' @param ExpObj A JAMS-style taxonomic SummarizedExperiment object (e.g. ConsolidatedGenomeBin).
#' @param glomby Taxonomic level to agglomerate to (e.g. "Species"). NULL for the native level.
#' @param only_allow_CSBs Restrict to Consolidated Species Bins (MetaBAT2 bins) only.
#' @param samplesToKeep Optional vector of samples to restrict to.
#' @param PPM_thresholds Optional explicit vector of PPM thresholds. If NULL, a sensible ladder of
#'   length n_ppm is built automatically from the data (see n_ppm).
#' @param completeness_thresholds Optional explicit vector of completeness thresholds. If NULL, an
#'   evenly spaced ladder of length n_comp from 0 to a data-driven maximum is built.
#' @param n_ppm Number of PPM thresholds to auto-generate when PPM_thresholds is NULL. The ladder
#'   starts at 0 and follows a roughly log spacing (0, then quantile-like steps up to a high value)
#'   so the informative low-PPM region is well sampled. Default 10.
#' @param n_comp Number of completeness thresholds to auto-generate when completeness_thresholds is
#'   NULL. Evenly spaced from 0. Default 10. Where the assay exceeds 100% - which it legitimately
#'   does, since both make_SummarizedExperiments and agglomerate_features sum it - roughly half
#'   as many further rungs are appended between 100 and the observed maximum, so that a universe
#'   whose features mostly sit above 100 still has something to titrate against.
#' @param max_PPM_threshold Ceiling on the auto-built PPM ladder, in PPM. Default 10000, i.e. 1%
#'   relative abundance. Beyond that an abundance floor selects the dominant members of a
#'   community rather than removing noise, which is a biological choice the user should make
#'   explicitly via featcutoff or PPM_thresholds. Only applies to the auto-built ladder; an
#'   explicit PPM_thresholds vector is used as given. Set to Inf for the old uncapped behaviour.
#' @param prevalence_pct Fixed prevalence (the "in at least X percent of samples" companion applied
#'   to BOTH the PPM and completeness rules). Held constant across the grid. Default 5. Set the
#'   partner argument prevalence_min_samples to override this with an absolute sample count.
#' @param prevalence_min_samples Optional absolute sample-count prevalence. If given (non-NULL), a
#'   feature must meet the threshold in at least this many samples, OVERRIDING prevalence_pct. This
#'   exists because "5% of samples" means very different things at n=10 vs n=300; a count is often
#'   the more honest criterion, and setting it to 1 recovers the "present in at least one sample"
#'   behaviour (useful for detecting a single-sample bug while still censoring empty/low-quality
#'   bins via the completeness axis, e.g. completeness >= 20 in >= 1 sample).
#' @param depth_quantile_probs Quantile of per-sample retained depth reported as pct_depth_lowq,
#'   alongside the median and the outright minimum. Default 0.1. The minimum is hostage to a
#'   single blank or failed library; this says whether the depth loss is one sample or many.
#' @param PPM_normalize_to_bases_sequenced Passed to filter_experiment for PPM computation.
#' @param ignoreunclassified Drop the unclassified/none feature before counting.
#' @param class_to_ignore Metadata classes to drop (signature parity; selection is by samplesToKeep).
#' @param return_plots If TRUE, returns list(grid = <df>, plots = <list>). If FALSE, returns the
#'   data frame invisibly.
#' @param plot_on_the_fly If TRUE, prints plots and returns NULL.
#' @param verbose Emit brief progress via flog.info.
#' @export

titrate_feature_counts <- function(ExpObj = NULL, glomby = NULL, only_allow_CSBs = FALSE, samplesToKeep = NULL, PPM_thresholds = NULL, completeness_thresholds = NULL, n_ppm = 30, n_comp = 30, max_PPM_threshold = 10000, prevalence_pct = NULL, prevalence_min_samples = 2, depth_quantile_probs = 0.1, PPM_normalize_to_bases_sequenced = TRUE, ignoreunclassified = TRUE, class_to_ignore = "N_A", return_plots = TRUE, plot_on_the_fly = FALSE, verbose = TRUE){

    require(ggplot2)

    #Vet + agglomerate + CSB-restrict up front, once.
    obj <- ExpObjVetting(ExpObj = ExpObj, samplesToKeep = samplesToKeep, featuresToKeep = NULL, only_allow_CSBs = only_allow_CSBs, glomby = glomby, class_to_ignore = class_to_ignore)

    analysis <- metadata(obj)$analysis
    analysisname <- if (!is.null(glomby)) glomby else analysis
    facet_label <- paste0(analysisname, if (only_allow_CSBs) " | CSB-only" else " | all-contigs")

    #Build the PPM matrix (for the abundance rule and % feature retention) AND the raw BaseCounts
    #matrix (sequencing depth, for the % depth retention). Both come from the same vetted object.
    baseobj <- filter_experiment(SEobj = obj, featcutoff = c(0, 0), samplesToKeep = NULL, featuresToKeep = NULL, only_allow_CSBs = FALSE, normalization = "relabund", PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced, GenomeCompletenessCutoff = c(0, 0), give_info = FALSE)
    PPMmat <- as.matrix(assays(baseobj)$PPM)
    BaseMat <- as.matrix(assays(baseobj)$BaseCounts)

    have_GC <- "GenomeCompleteness" %in% names(assays(baseobj))
    GCmat <- if (have_GC) as.matrix(assays(baseobj)$GenomeCompleteness) else NULL

    #Captured before the unclassified drop below, so that it reflects the object as it stands.
    nfeat_total_before_drop <- nrow(PPMmat)

    #Drop unclassified/none rows from ALL matrices consistently.
    if (ignoreunclassified){
        dunno <- c(paste(analysis, "none", sep = "_"), "LKT__d__Unclassified", "LKT__Unclassified", paste0(analysisname, "__Unclassified"))
        keep <- which(!(rownames(PPMmat) %in% dunno) & rownames(PPMmat) != "" & !is.na(rownames(PPMmat)))
        PPMmat <- PPMmat[keep, , drop = FALSE]
        BaseMat <- BaseMat[rownames(PPMmat), , drop = FALSE]
        if (have_GC) GCmat <- GCmat[rownames(PPMmat), , drop = FALSE]
    }

    #DROP PHANTOM FEATURES: rows with zero counts in EVERY sample of this object.
    #
    #A feature that was never observed here is not a feature of this object, and including it
    #corrupts the two things the fit is built on. It pads nfeat_total, so pct_retained - and
    #therefore divergence, which is depth% minus feature% - is compressed on the feature side,
    #biasing the whole surface towards looser cells. And it poisons the noise-fraction test that
    #decides whether the abundance axis is meaningful, because a phantom row has no genome
    #evidence and is counted as noise.
    #
    #This is not hypothetical and it is not rare. JAMS objects acquire phantom rows whenever
    #feature spaces are unioned or samples are dropped: merge_JAMS_SEobj unions the features of
    #both cohorts, and subsetting samples afterwards leaves behind every feature that belonged
    #only to the samples that went. On a 178-sample multi-cohort human stool object, 32124 of
    #46691 rows - 68.8% - were phantoms. Restricted to Consolidated Species Bins the figure was
    #11305 of 13651, and the consequence was concrete: the noise-fraction test read 82.8%
    #instead of the true 0.0%, left the abundance axis switched on, and drove a universe of
    #binned genomes to a ~1% relative-abundance floor. With the phantoms gone the same object
    #filters purely on genome completeness, exactly as the MR50 benchmark does. The MR50 cohort
    #has no phantom rows at all, which is why this went unseen until a merged object was tried.
    n_phantom <- 0L
    if (nrow(BaseMat) > 0){
        real_rows <- which(rowSums(BaseMat, na.rm = TRUE) > 0)
        n_phantom <- nrow(BaseMat) - length(real_rows)
        if ((n_phantom > 0) && (length(real_rows) > 0)){
            PPMmat <- PPMmat[real_rows, , drop = FALSE]
            BaseMat <- BaseMat[real_rows, , drop = FALSE]
            if (have_GC) GCmat <- GCmat[real_rows, , drop = FALSE]
        }
    }

    nsamp <- ncol(PPMmat)
    nfeat_total <- nrow(PPMmat)
    #Row count of the vetted (and, where asked, agglomerated) object BEFORE the unclassified drop
    #above. This is the number of rows the object actually has at this rank, and it is recorded so
    #that fit_filtering_thresholds can fingerprint the feature universe a fit belongs to. Without
    #it, an object the user agglomerated by hand before plotting still reports its native
    #metadata$analysis, so a stamped native-rank entry would be found and applied at the wrong
    #rank with nothing to detect the mismatch.
    nfeat_object <- nfeat_total_before_drop

    #IMPORTANT for honest depth accounting: % depth retained is measured against the TOTAL depth of
    #the (vetted) object BEFORE the ignoreunclassified drop would understate the denominator. We use
    #the per-sample total of the kept feature universe as the denominator, i.e. the same feature set
    #the grid draws from, so retention fractions are internally consistent with the count panel.
    sample_total_depth <- colSums(BaseMat)
    #Guard against any all-zero sample column (would divide by zero).
    sample_total_depth[sample_total_depth == 0] <- NA_real_

    #Resolve the prevalence rule: absolute count takes precedence over percentage.
    #PARITY WITH filter_experiment IS LOAD-BEARING HERE. fit_filtering_thresholds chooses a grid
    #cell on the counts produced below and then hands the corresponding thresholds to
    #filter_experiment, so if the two functions disagree about who survives, the fitted record
    #promises a feature count that the filter does not deliver. filter_experiment keeps a feature
    #when (passing samples / ncol) >= (sampcutoffpct / 100). The test below is written in exactly
    #that form - a PROPORTION against a proportion - rather than as a count against
    #(pct/100)*nsamp, because the two are algebraically identical but not numerically identical:
    #at nsamp 225 and pct 28, (28/100)*225 is 63.00000000000001 so a count of 63 fails, while
    #63/225 >= 0.28 passes. The default prevalence of 5% happens to be clean at every cohort size
    #tested, which is exactly why this would have gone unnoticed until someone changed it.
    #The one place the two used to part company outright is a prevalence of zero: filter_experiment
    #converts that to "at least one sample" (and warns), whereas a bare threshold of 0 here would
    #pass every feature including those seen in no sample at all. Floored to the half-sample
    #margin so the two agree.
    #RESOLVE THE PREVALENCE RULE, AND THEN NEVER MENTION IT AGAIN.
    #
    #The filter this titration explores is "value >= X in >= k samples". For integer k that is
    #exactly equivalent to "the feature's k-th LARGEST value across samples is >= X", and the
    #rest of this function uses that equivalent form. Two reasons. It is far cheaper - one sort
    #per axis instead of a full pass per grid cell - and, more importantly, it makes the
    #prevalence rule visible. The k-th largest value is the single number the threshold is
    #compared against, so a feature whose k-th largest is ZERO is removed by PREVALENCE, at every
    #rung of the ladder, no matter what the threshold is. That population is counted and reported
    #below, because it was invisible before and it turned out to be enormous: on the MR50
    #benchmark restricted to Consolidated Species Bins, 291 of 384 bins were being discarded by
    #the prevalence rule alone, every one of them >= 70% complete wherever it appeared, while the
    #log said they had failed a GENOME COMPLETENESS threshold of 11%.
    #
    #Prevalence now defaults to an absolute COUNT rather than a percentage. Whether a detection
    #is credible does not scale with cohort size - two independent observations are two
    #independent observations whether the study has 20 samples or 2000 - whereas "5% of samples"
    #silently means one sample at n <= 20, three at n = 59 and fifteen at n = 300. JAMS also has
    #something most tools lack, a direct quality signal in genome completeness, so it does not
    #need prevalence as a proxy for credibility. prevalence_pct is still honoured when it is
    #given explicitly and prevalence_min_samples is NULL.
    if (!is.null(prevalence_min_samples)){
        min_samples <- max(1, min(round(prevalence_min_samples), nsamp))
        prev_desc <- paste0(">= ", min_samples, " sample(s)")
    } else if (!is.null(prevalence_pct)){
        #Smallest integer count satisfying count/nsamp >= pct/100, floored at one sample so that
        #a prevalence of zero cannot admit a feature seen in no sample at all (filter_experiment
        #makes the same reinterpretation, with a warning).
        min_samples <- max(1, ceiling((min(prevalence_pct, 100) / 100) * nsamp - 1e-9))
        prev_desc <- paste0(">= ", prevalence_pct, "% of samples (= ", min_samples, " of ", nsamp, ")")
    } else {
        min_samples <- max(1, min(2, nsamp))
        prev_desc <- paste0(">= ", min_samples, " sample(s)")
    }
    #The percentage that filter_experiment must be given so that it readmits exactly this count.
    #It keeps a feature when (passing / ncol) >= (pct / 100); counts are integers, so putting the
    #percentage HALFWAY between min_samples - 1 and min_samples admits exactly min_samples or
    #more and cannot be flipped by floating-point round-trip.
    prevalence_as_pct <- ((min_samples - 0.5) / nsamp) * 100

    #The per-feature statistic each axis is thresholded against.
    kth_largest <- function(M){
        if (min_samples >= ncol(M)){
            return(apply(M, 1, min))
        }

        return(apply(M, 1, function(x){ sort(x, decreasing = TRUE, na.last = TRUE)[min_samples] }))
    }
    PPMk <- kth_largest(PPMmat)
    PPMk[!is.finite(PPMk)] <- 0
    if (have_GC){
        GCk <- kth_largest(GCmat)
        GCk[!is.finite(GCk)] <- 0
        GCmax <- suppressWarnings(apply(GCmat, 1, max, na.rm = TRUE))
        GCmax[!is.finite(GCmax)] <- 0
    } else {
        GCk <- rep(0, nfeat_total)
        GCmax <- rep(0, nfeat_total)
    }
    #Two populations worth naming. "No genome evidence" is a feature whose completeness is zero
    #in every sample: unambiguous noise, and the thing denoising exists to remove. "Lost to
    #prevalence" is a feature that DOES carry genome evidence but in fewer than min_samples
    #samples, so no threshold on any axis can save it.
    n_no_evidence <- sum(GCmax <= 0)
    n_lost_to_prevalence <- sum((GCmax > 0) & (GCk <= 0))

    #Build the ladders from the QUANTILES of those same per-feature statistics, so that each rung
    #sits where features actually cross it. Even spacing was badly wrong in both directions: at
    #the native rank the non-zero completeness values pile up between 2% and 8%, so the single
    #rung at 11% was above essentially all of them, while on a Consolidated Species Bin universe
    #the minimum completeness is 70% and the six rungs from 11 to 67 all returned the identical
    #answer. Worse, spacing the completeness axis evenly and the PPM axis logarithmically made
    #the two incomparable, so which axis the fit chose depended partly on how finely each had
    #happened to be sampled - capping the PPM ladder once re-spaced it, tripled its rungs in the
    #band that mattered, and flipped a fit from completeness-based to abundance-based on its own.
    quantile_ladder <- function(stat, n, cap = NULL){
        pos <- stat[is.finite(stat) & (stat > 0)]
        if (length(pos) < 1){
            return(0)
        }
        rungs <- unique(round(as.numeric(quantile(pos, probs = seq(0, 1, length.out = max(2, n)), names = FALSE))))
        rungs <- sort(unique(c(0, rungs)))
        if (!is.null(cap) && is.finite(cap)){
            rungs <- rungs[rungs <= cap]
        }
        if (length(rungs) < 1){
            rungs <- 0
        }

        return(rungs)
    }

    if (is.null(PPM_thresholds)){
        PPM_thresholds <- quantile_ladder(PPMk, n_ppm, cap = max_PPM_threshold)
    }

    if (is.null(completeness_thresholds)){
        completeness_thresholds <- if (have_GC) quantile_ladder(GCk, n_comp) else 0
    }

    if (!have_GC && any(completeness_thresholds > 0)){
        if (verbose) flog.warn("No GenomeCompleteness assay present; completeness thresholds > 0 collapsed to 0.")
        completeness_thresholds <- 0
    }

    #Walk the grid one COMPLETENESS rung at a time, and inside each rung use a running
    #cumulative sum rather than re-subsetting the count matrix for every abundance rung.
    #
    #Why it is exact. Within the features that clear a completeness rung, order them by their
    #abundance statistic from largest to smallest. The survivors of abundance rung p are then
    #precisely the first m of them, where m is how many have PPMk >= p, so their depth is row m
    #of the cumulative sum. No approximation, just the same arithmetic done once instead of once
    #per cell. It turns the cost from (rungs x rungs x features x samples) into
    #(rungs x features x samples), which is what makes a 30 x 30 grid affordable - and 30 is not
    #a luxury: at 20 rungs the moderate and stringent tiers land on the same cell at several
    #ranks of the MR50 benchmark, because the targets they are chasing fall between rungs.
    ppm_levels_sorted <- sort(unique(PPM_thresholds))
    comp_levels_sorted <- sort(unique(completeness_thresholds))
    grid <- expand.grid(PPM = ppm_levels_sorted, Completeness = comp_levels_sorted, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    grid$n_features <- 0L
    grid$pct_depth_median <- 0
    grid$pct_depth_lowq <- NA_real_
    grid$pct_depth_min <- NA_real_

    for (c_thr in comp_levels_sorted){
        sel <- if (have_GC && (c_thr > 0)) which(GCk >= c_thr) else seq_len(nfeat_total)
        if (length(sel) > 0){
            ord <- sel[order(PPMk[sel], decreasing = TRUE)]
            pk_sorted <- PPMk[ord]
            csum <- apply(BaseMat[ord, , drop = FALSE], 2, cumsum)
            #apply() drops to a vector when there is a single surviving feature.
            if (is.null(dim(csum))){
                csum <- matrix(csum, nrow = 1, dimnames = list(NULL, names(csum)))
            }
        } else {
            ord <- integer(0)
            pk_sorted <- numeric(0)
            csum <- NULL
        }
        for (p_thr in ppm_levels_sorted){
            m <- if (length(ord) < 1) 0L else if (p_thr <= 0) length(ord) else sum(pk_sorted >= p_thr)
            here <- which((grid$PPM == p_thr) & (grid$Completeness == c_thr))
            grid$n_features[here] <- m
            if (m > 0){
                depth_frac <- csum[m, ] / sample_total_depth
                grid$pct_depth_median[here] <- round(100 * stats::median(depth_frac, na.rm = TRUE), 1)
                lowq <- suppressWarnings(as.numeric(stats::quantile(depth_frac, probs = depth_quantile_probs[1], na.rm = TRUE, names = FALSE)))
                grid$pct_depth_lowq[here] <- if (is.finite(lowq)) round(100 * lowq, 1) else NA_real_
                dmin <- suppressWarnings(min(depth_frac, na.rm = TRUE))
                grid$pct_depth_min[here] <- if (is.finite(dmin)) round(100 * dmin, 1) else NA_real_
            } else {
                grid$pct_depth_median[here] <- 0
                grid$pct_depth_lowq[here] <- 0
                grid$pct_depth_min[here] <- 0
            }
        }
    }
    grid$pct_retained <- round(100 * grid$n_features / nfeat_total, 1)
    grid$facet <- facet_label

    if (verbose){
        flog.info(sprintf("%s: counted survivors over %d PPM x %d completeness cells (%d features, %d samples). %d feature(s) carry no genome evidence in any sample; a further %d do carry it but in fewer than the %d sample(s) the prevalence rule requires, so no threshold on either axis can retain them.", facet_label, length(ppm_levels_sorted), length(comp_levels_sorted), nfeat_total, nsamp, n_no_evidence, n_lost_to_prevalence, min_samples))
        if (n_phantom > 0){
            flog.info(sprintf("%s: %d of %d rows (%.1f%%) had zero counts in every sample and were excluded from the titration - they are not features of this object. Rows like these accumulate when feature spaces are unioned by merge_JAMS_SEobj or when samples are dropped from a merged object; they are harmless to a plot but they distort any statistic computed over the feature universe.", facet_label, n_phantom, n_phantom + nfeat_total, 100 * n_phantom / (n_phantom + nfeat_total)))
        }
    }

    #A fit on a handful of samples is a fit on a handful of observations of each feature, and the
    #quantile ladders and the k-th largest statistic are both estimated from them. On the MR50
    #benchmark the chosen cell is reproducible across random subsamples from about twenty samples
    #upward and wanders badly below that, occasionally collapsing to no filtration at all.
    if ((nsamp < 20) && verbose){
        flog.warn(sprintf("%s: thresholds are being fitted on only %d samples. Below about twenty the fitted cell is not reproducible - on repeated subsamples of a benchmark cohort it wandered across a wide range and sometimes returned no filtration at all. Treat the result as indicative, set featcutoff and GenomeCompletenessCutoff explicitly, or reuse thresholds stamped from a larger cohort.", facet_label, nsamp))
    }

    #Local rate of decline: fractional feature loss relative to the looser (one-step-down) neighbour
    #in each direction, averaged. Small magnitude = a stable region (tightening barely costs features);
    #large magnitude = a steep face (each notch removes many features). This is where turning points
    #such as the ~10% completeness floor become visible as a drop-off in the rate.
    ppm_levels <- sort(unique(grid$PPM))
    comp_levels <- sort(unique(grid$Completeness))
    key <- function(p, c) paste(p, c, sep = "§")
    lut <- setNames(grid$n_features, key(grid$PPM, grid$Completeness))
    grid$decline_rate <- NA_real_
    for (i in 1:nrow(grid)){
        p <- grid$PPM[i]; c <- grid$Completeness[i]
        here <- lut[[key(p, c)]]
        rates <- c()
        pi <- match(p, ppm_levels); ci <- match(c, comp_levels)
        if (pi > 1){ prev <- lut[[key(ppm_levels[pi - 1], c)]]; if (prev > 0) rates <- c(rates, (prev - here) / prev) }
        if (ci > 1){ prev <- lut[[key(p, comp_levels[ci - 1])]]; if (prev > 0) rates <- c(rates, (prev - here) / prev) }
        if (length(rates) > 0) grid$decline_rate[i] <- mean(rates)
    }

    #Carried on the grid itself so that any consumer (notably fit_filtering_thresholds, which
    #records them as provenance) can state exactly which universe these counts were taken over
    #without having to re-derive it. The unclassified flag matters: filter_experiment does NOT
    #drop the unclassified bin, so when ignoreunclassified is TRUE the counts here exclude a row
    #that will still be present after filtration. Set after the decline-rate loop so that nothing
    #downstream can strip it by column assignment.
    attr(grid, "titration_universe") <- list(analysis = analysis, analysisname = analysisname, glomby = glomby, only_allow_CSBs = only_allow_CSBs, nfeat_total = nfeat_total, nfeat_object = nfeat_object, n_phantom = n_phantom, nsamp = nsamp, samplenames = colnames(PPMmat), min_samples = min_samples, prevalence_as_pct = prevalence_as_pct, prevalence_pct = prevalence_pct, prevalence_min_samples = prevalence_min_samples, prev_desc = prev_desc, n_no_evidence = n_no_evidence, n_lost_to_prevalence = n_lost_to_prevalence, ignoreunclassified = ignoreunclassified, have_GenomeCompleteness = have_GC, PPM_normalize_to_bases_sequenced = PPM_normalize_to_bases_sequenced)

    if (!return_plots){
        return(invisible(grid))
    }

    gdf <- grid
    gdf$PPM_f <- factor(gdf$PPM, levels = ppm_levels)
    gdf$Comp_f <- factor(gdf$Completeness, levels = comp_levels)

    tile_theme <- theme_minimal() + theme(plot.title = element_text(size = 9), axis.text = element_text(size = 7), axis.text.x = element_text(angle = 45, hjust = 1))
    subt <- paste0(facet_label, " | prevalence ", prev_desc, " | ", nfeat_total, " starting features | n = ", nsamp, " samples")

    p_count <- ggplot(gdf, aes(x = PPM_f, y = Comp_f, fill = n_features)) +
        geom_tile(colour = "grey80") +
        geom_text(aes(label = sprintf("%d\n%.0f%%", n_features, pct_retained)), size = 2.3) +
        scale_fill_gradient(low = "#f7fbff", high = "#08306b") +
        labs(x = "min PPM threshold", y = "min genome completeness (%)", fill = "features", title = paste0("Surviving feature count (count and % of starting retained)\n", subt)) +
        tile_theme

    #Decline-rate surface: dark = steep (many features lost per notch), pale = stable plateau.
    p_rate <- ggplot(gdf, aes(x = PPM_f, y = Comp_f, fill = decline_rate)) +
        geom_tile(colour = "grey80") +
        geom_text(aes(label = ifelse(is.na(decline_rate), "", sprintf("%.0f%%", 100 * decline_rate))), size = 2.3) +
        scale_fill_gradient(low = "#f7f7f7", high = "#b30000", na.value = "grey90", labels = scales::percent) +
        labs(x = "min PPM threshold", y = "min genome completeness (%)", fill = "feature loss\nper step", title = paste0("Local rate of feature loss (pale = stable plateau; dark = steep)\n", subt)) +
        tile_theme

    #Depth-retention surface: what % of sequencing depth the survivors represent. Median across
    #samples (fill and top number), with the WORST-hit sample's retention beneath (min). A cell can
    #retain few FEATURES yet high DEPTH (axed only noise) or, alarmingly, low depth (axed real
    #community). The min matters because a sample retaining, say, 15% depth is effectively gutted
    #even if the median looks comfortable.
    p_depth <- ggplot(gdf, aes(x = PPM_f, y = Comp_f, fill = pct_depth_median)) +
        geom_tile(colour = "grey80") +
        geom_text(aes(label = sprintf("%.0f%%\n(low %.0f%% / min %.0f%%)", pct_depth_median, pct_depth_lowq, pct_depth_min)), size = 2.0) +
        scale_fill_gradient(low = "#fff5eb", high = "#7f2704", limits = c(0, 100)) +
        labs(x = "min PPM threshold", y = "min genome completeness (%)", fill = "% depth\n(median)", title = paste0("Sequencing-depth retained by survivors (median across samples; worst sample in parens)\n", subt)) +
        tile_theme

    #Divergence surface: the gap between feature loss and depth loss. This is the direct read on the
    #evenness question. depth_retained - feature_retained, in percentage points:
    #  large POSITIVE (blue) = kept far more depth than features => you axed low-abundance noise (good);
    #  near ZERO (white)     = features and depth fell together => you axed a representative slice;
    #  NEGATIVE (red)        = lost more depth than features => you axed high-abundance community (bad).
    gdf$divergence <- gdf$pct_depth_median - gdf$pct_retained
    p_diverge <- ggplot(gdf, aes(x = PPM_f, y = Comp_f, fill = divergence)) +
        geom_tile(colour = "grey80") +
        geom_text(aes(label = sprintf("%+.0f", divergence)), size = 2.3) +
        scale_fill_gradient2(low = "#b2182b", mid = "#f7f7f7", high = "#2166ac", midpoint = 0) +
        labs(x = "min PPM threshold", y = "min genome completeness (%)", fill = "depth% -\nfeature%", title = paste0("Depth retained minus features retained (blue = axed mostly noise; red = axed real depth)\n", subt)) +
        tile_theme

    titrate_feature_counts_list <- list(grid = grid, plots = list(feature_count = p_count, decline_rate = p_rate, depth_retained = p_depth, count_vs_depth_divergence = p_diverge))

    if (plot_on_the_fly){
        print(titrate_feature_counts_list$plots)
        return(NULL)
    } else {
        return(titrate_feature_counts_list)
    }

}